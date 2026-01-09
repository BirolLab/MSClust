#include <algorithm>
#include <cmath>
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <iomanip>
#include <map>
#include <omp.h>
#include <set>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <argparse/argparse.hpp>

#include "btllib/bloom_filter.hpp"
#include "btllib/counting_bloom_filter.hpp"
#include "btllib/bshash.hpp"
#include <btllib/seq.hpp>
#include <btllib/seq_reader.hpp>

#include "city.h"

#include <regex>
#include <filesystem>

namespace fs = std::filesystem;

// --- Logarithmic Quality Averaging ---
double calculate_avg_phred(const std::string& qual) {
    if (qual.empty()) return 0.0;
    double sum_prob = 0.0;
    for (char c : qual) {
        int q = static_cast<int>(c) - 33;
        sum_prob += std::pow(10.0, -q / 10.0);
    }
    double avg_prob = sum_prob / qual.size();
    if (avg_prob <= 0) return 99.0; 
    return -10.0 * std::log10(avg_prob);
}

// Generate both C->T and G->A converted k-mers
std::vector<std::string> generate_converted_kmers(const std::string& kmer, bool dev) {
    std::string kmer_upper = kmer;
    std::transform(kmer_upper.begin(), kmer_upper.end(), kmer_upper.begin(), ::toupper);

    char meth_base = dev ? '1' : 'C';

    std::string ct = kmer_upper;
    std::replace(ct.begin(), ct.end(), meth_base, 'T');

    std::string ga = kmer_upper;
    std::replace(ga.begin(), ga.end(), 'G', 'A');

    return {ct, ga};
}

// Pairing logic
std::map<std::string, std::pair<std::string, std::string>> pair_fastq_files(const std::vector<std::string>& files) {
    std::map<std::string, std::pair<std::string, std::string>> paired_files;

    std::regex suffix_re("(_R?[12]|_[01]{3})\\.f(ast)?q(\\.gz)?$", std::regex::icase); // matches _1, _R1, _001, etc.

    for (const auto& file : files) {
        std::string filename = fs::path(file).filename().string();

        std::smatch match;
        if (std::regex_search(filename, match, suffix_re)) {
            std::string prefix = filename.substr(0, match.position());  // base name
            std::string suffix = match.str(1); // _1, _R2, etc.

            auto& pair = paired_files[prefix];

            if (suffix.find("1") != std::string::npos)
                pair.first = file;
            else if (suffix.find("2") != std::string::npos)
                pair.second = file;
        } else  {
            // No R1/R2 pattern → treat as single/interleaved
            paired_files[filename].first = file;  // just filename, second empty
        }
    }

    return paired_files;
}

// 1-mer entropy
float shannon_entropy(std::string_view s) {
    if (s.empty()) return 0.0f;

    int counts[256] = {};
    for (unsigned char c : s) {
        counts[c]++;
    }

    const float length = static_cast<float>(s.size());
    float entropy = 0.0f;

    for (int count : counts) {
        if (count == 0) continue;
        float p = count / length;
        entropy -= p * std::log2f(p);
    }

    return entropy;
}

// 3-mer entropy
float shannon_entropy_trimer(std::string_view s) {
    if (s.size() < 3) return 0.0f;

    std::unordered_map<std::string_view, int> counts;
    const size_t total = s.size() - 2;

    for (size_t i = 0; i < total; ++i) {
        std::string_view kmer{s.data() + i, 3};
        counts[kmer]++;
    }

    const float denom = static_cast<float>(total);
    float entropy = 0.0f;

    for (const auto& [kmer, count] : counts) {
        float p = static_cast<float>(count) / denom;
        entropy -= p * std::log2f(p);
    }

    return entropy;
}

// 2-mer entropy
float shannon_entropy_dimer(std::string_view s) {
    if (s.size() < 2) return 0.0f;

    std::unordered_map<std::string_view, int> counts;
    const size_t total = s.size() - 1;

    for (size_t i = 0; i < total; ++i) {
        std::string_view kmer{s.data() + i, 2};
        counts[kmer]++;
    }

    const float denom = static_cast<float>(total);
    float entropy = 0.0f;

    for (const auto& [kmer, count] : counts) {
        float p = static_cast<float>(count) / denom;
        entropy -= p * std::log2f(p);
    }

    return entropy;
}


std::vector<std::pair<uint64_t, bool>> get_all_methylation_kmers(
    const std::string& seq,
    unsigned k,
    const btllib::BloomFilter& clean_ct_mers,
    const btllib::BloomFilter& clean_ga_mers//, const float shannon,  const float shannon2,  const float shannon3, const std::string&  qual,int phred_threshold
) {
    
    std::vector<std::pair<uint64_t, bool>> all_kmers_hash;
    btllib::BsHashDirectional bh(seq, 1, k, "CT"); 
    btllib::BsHashDirectional bh_ga(seq, 1, k, "GA"); 
    while (bh.roll() && bh_ga.roll()) { 
        bool is_methylated = false;
        auto central_dimer = bh.center_dimer();
        if (central_dimer == "TG" || central_dimer == "CG" || central_dimer == "CA") {
            if (central_dimer == "CG") {
                is_methylated = true; // double check it exists here, if not flag kmer and hash, add entropy and kmer check
                if (clean_ct_mers.contains(bh.hashes())) { 
                    all_kmers_hash.push_back(std::make_pair(bh.hashes()[0], is_methylated));
                }
                if (clean_ga_mers.contains(bh_ga.hashes())) { 
                    all_kmers_hash.push_back(std::make_pair(bh_ga.hashes()[0], is_methylated));
                } else {
                    continue;
                }
            } 
            if (central_dimer == "TG") {
                if (clean_ct_mers.contains(bh.hashes())) { 
                    all_kmers_hash.push_back(std::make_pair(bh.hashes()[0], is_methylated));
                }
            } 
            if (central_dimer == "CA") {
                if (clean_ga_mers.contains(bh_ga.hashes())) { 
                    all_kmers_hash.push_back(std::make_pair(bh_ga.hashes()[0], is_methylated));
                }
            } 
	    }
    }
    return all_kmers_hash;
}

void export_matrices(const std::vector<std::vector<uint8_t>>& bfs1, 
                     const std::vector<std::vector<uint8_t>>& methylated_bfs1, 
                     const std::vector<size_t>& final_indices,
                     const std::vector<std::string>& sample_names,
                     bool dev) {
    
    size_t numSamples = bfs1.size();
    std::vector<size_t> cuts = {1000, 2000, 5000, 50000};
    std::vector<std::ofstream> files(cuts.size());
    
    for (size_t f = 0; f < cuts.size(); ++f) {
        size_t actual_cut = std::min(cuts[f], final_indices.size());
        files[f].open("methylation_top_" + std::to_string(actual_cut) + ".csv");
        
        // Header: Add "Sample_ID" as the first column
        files[f] << "Sample_ID";
        for (size_t j = 0; j < actual_cut; ++j) {
            files[f] << ",Site_" << final_indices[j];
        }
        files[f] << "\n";
    }

    for (size_t i = 0; i < numSamples; ++i) {
        // --- START OF YOUR NAME CLEANING LOGIC ---
        std::string raw_name = sample_names[i];
        std::string clean_name;

        if (dev) {
            size_t pos = raw_name.find("GSM");
            size_t pos2 = raw_name.find("_aligned_reads.fasta");
            if (pos != std::string::npos && pos2 != std::string::npos)
                clean_name = raw_name.substr(pos, pos2 - pos);
            else
                clean_name = raw_name;
        } else {
            size_t pos = raw_name.find_last_of("/");
            size_t pos2 = raw_name.find_first_of(".fq");
            if (pos2 == std::string::npos) pos2 = raw_name.find_first_of(".fastq");
            
            size_t start = (pos == std::string::npos) ? 0 : pos + 1;
            if (pos2 != std::string::npos && pos2 > start)
                clean_name = raw_name.substr(start, pos2 - start);
            else
                clean_name = raw_name.substr(start);
        }
        // --- END OF NAME CLEANING LOGIC ---

        // Prepare the row data
        std::vector<int> sample_row(final_indices.size());
        #pragma omp parallel for schedule(static)
        for (size_t j = 0; j < final_indices.size(); ++j) {
            size_t site_idx = final_indices[j];
            if (bfs1[i][site_idx] == 0) {
                sample_row[j] = 0;
            } else {
                sample_row[j] = (methylated_bfs1[i][site_idx] == 1) ? 1 : -1;
            }
        }

        // Write to files
        for (size_t f = 0; f < cuts.size(); ++f) {
            size_t actual_cut = std::min(cuts[f], final_indices.size());
            files[f] << clean_name; // Prepend Sample Name
            for (size_t j = 0; j < actual_cut; ++j) {
                files[f] << "," << sample_row[j];
            }
            files[f] << "\n";
        }
    }

    for (auto& f : files) f.close();
}

int main(int argc, char* argv[]) {
    argparse::ArgumentParser program("example");

    program.add_argument("-i")
        .required()
        .help("specify the first input file");
    program.add_argument("-o")
        .help("soutput_prefix")
        .default_value("_");

    program.add_argument("-j")
        .required()
        .help("specify the second input file");

    program.add_argument("-t")
        .scan<'i', int>()
        .default_value(1)
        .help("specify the number of threads (default is 1)");
    program.add_argument("-q")
        .scan<'i', int>()
        .default_value(30)
        .help("specify the min Phred (default is 30)");
    program.add_argument("-s")
        .scan<'f', float>()
        .default_value(1.1f)
        .help("Shannon Entropy (default is 1.1)");
    program.add_argument("--s2")
        .scan<'f', float>()
        .default_value(2.2f)
        .help("Dimer Shannon Entropy (default is 2.2)");
    program.add_argument("--s3")
        .scan<'f', float>()
        .default_value(3.3f)
        .help("Trimer Shannon Entropy (default is 3.3)");
    program.add_argument("-k", "--kmer")
        // unsigned
        .scan<'u', unsigned>()
        .default_value(25u)
        .help("specify the kmer size (default is 25)");
    program.add_argument("-m", "--min")
        // unsigned
        .scan<'u', unsigned>()
        .default_value(3u)
        .help("specify the min kmer");
    program.add_argument("-n", "--max")
        // unsigned
        .scan<'u', unsigned>()
        .default_value(3u)
        .help("specify the max kmer");
    // add a development flag
    program.add_argument("-d", "--dev")
        .default_value(false)
        .implicit_value(true)
        .help("enable development mode (default is false)");

        program.add_argument("-c", "--complexity")
        .default_value(false)
        .implicit_value(true)
        .help("calculate complexity (default is false)");


    try {
        program.parse_args(argc, argv);
    } catch (const std::runtime_error& err) {
        std::cerr << err.what() << std::endl;
        std::cerr << program;
        return 1;
    }

    std::string inputFile1 = program.get<std::string>("-i");
    std::string inputFile2 = program.get<std::string>("-j");
    std::string prefix = program.get<std::string>("-o");
    bool dev = program.get<bool>("-d");
    bool complexity = program.get<bool>("-c");
    int numThreads = program.get<int>("-t");
    float shannon = program.get<float>("-s");
    float shannon2 = program.get<float>("--s2");
    float shannon3 = program.get<float>("--s3");
    unsigned minKmer = program.get<unsigned>("-m");
    int phred_threshold = program.get<int>("-q");
    unsigned maxKmer = program.get<unsigned>("-n");
    unsigned k = program.get<unsigned>("-k");
omp_set_num_threads(numThreads);

    std::cout << "Input file 1: " << inputFile1 << "\n";
    std::cout << "Input file 2: " << inputFile2 << "\n";
    std::cout << "Output prefix: " << prefix << "\n";
    std::cout << "Development mode: " << (dev ? "enabled" : "disabled") << "\n";
    std::cout << "Complexity filtering: " << (complexity ? "enabled" : "disabled") << "\n";
    std::cout << "Shannon entropy threshold: " << shannon << "\n";
    std::cout << "Dimer Shannon entropy threshold: " << shannon2 << "\n";
    std::cout << "Trimer Shannon entropy threshold: " << shannon3 << "\n";
    std::cout << "Phred avg threshold: " << phred_threshold << "\n";
    std::cout << "k-mer size : " << k << "\n";
    std::cout << "Minimum k-mer occurence: " << minKmer << "\n";
    std::cout << "Maximum k-mer occurence: " << maxKmer << "\n";


    // read the first input file and store each line into a vector
    std::vector<std::string> lines1;
    std::ifstream file1(inputFile1);
    if (file1.is_open()) {
        std::string line;
        while (std::getline(file1, line)) {
            lines1.push_back(line);
        }
        file1.close();
    } else {
        std::cerr << "Unable to open file " << inputFile1 << "\n";
        return 1;
    }

    // read the second input file and store each line into a vector
    std::vector<std::string> lines2;
    std::ifstream file2(inputFile2);
    if (file2.is_open()) {
        std::string line;
        while (std::getline(file2, line)) {
            lines2.push_back(line);
        }
        file2.close();
    } else {
        std::cerr << "Unable to open file " << inputFile2 << "\n";
        return 1;
    }

    auto pairs = pair_fastq_files(lines1);

    if (complexity){
        std::cerr << "logging complexity"  << std::endl;
        int num_lines_2 = 0;

        std::ofstream cg_log("cg_complexity.tsv");
        std::ofstream cg2_log("cg2_complexity.tsv");
        std::ofstream cg3_log("cg3_complexity.tsv");
        /*methy_kmer_log << "error_rate\n";  //
        */
        for (const auto& [prefix, pair] : pairs) {
            const auto& r1_file = pair.first;
            const auto& r2_file = pair.second;

            std::cerr << num_lines_2 << std::endl;
            num_lines_2++;

            btllib::SeqReader reader(r1_file, btllib::SeqReader::Flag::SHORT_MODE);


            #pragma omp parallel
            for (const auto record : reader) {
            if (record.seq.size() < k) {
                continue;
            }
                for (size_t j = 0; j + k <= record.seq.size(); ++j) {
                    char meth_base = dev ? '1' : 'C';
                    if (record.seq[j + k / 2 - 1] == meth_base && record.seq[j + k / 2] == 'G') {
                        std::string orig_kmer = record.seq.substr(j, k);
                        auto converted_kmers = generate_converted_kmers(orig_kmer, dev);
                        float cg1 = shannon_entropy(converted_kmers[0]);
                        float cg2 = shannon_entropy_dimer(converted_kmers[0]);
                        float cg3 = shannon_entropy_trimer(converted_kmers[0]);
                        float cg1_rc = shannon_entropy(converted_kmers[1]);
                        float cg2_rc = shannon_entropy_dimer(converted_kmers[1]);
                        float cg3_rc = shannon_entropy_trimer(converted_kmers[1]);
            #pragma omp critical 
            {
                        cg_log << cg1 << std::endl;
                        cg2_log << cg2 << std::endl;
                        cg3_log << cg3 << std::endl;
                        cg_log << cg1_rc << std::endl;
                        cg2_log << cg2_rc << std::endl;
                        cg3_log << cg3_rc << std::endl;
            }
                            
                    }
                }
            }
            if (!pair.second.empty()) {
                btllib::SeqReader reader2(r2_file, btllib::SeqReader::Flag::SHORT_MODE);


                #pragma omp parallel
                for (const auto record : reader2) {
                    for (size_t j = 0; j + k <= record.seq.size(); ++j) {
                        char meth_base = dev ? '1' : 'C';
                        if (record.seq[j + k / 2 - 1] == meth_base && record.seq[j + k / 2] == 'G') {
                        std::string orig_kmer = record.seq.substr(j, k);
                        auto converted_kmers = generate_converted_kmers(orig_kmer, dev);
                        float cg1 = shannon_entropy(converted_kmers[0]);
                        float cg2 = shannon_entropy_dimer(converted_kmers[0]);
                        float cg3 = shannon_entropy_trimer(converted_kmers[0]);
                        float cg1_rc = shannon_entropy(converted_kmers[1]);
                        float cg2_rc = shannon_entropy_dimer(converted_kmers[1]);
                        float cg3_rc = shannon_entropy_trimer(converted_kmers[1]);
                #pragma omp critical 
                {
                        cg_log << cg1 << std::endl;
                        cg2_log << cg2 << std::endl;
                        cg3_log << cg3 << std::endl;
                        cg_log << cg1_rc << std::endl;
                        cg2_log << cg2_rc << std::endl;
                        cg3_log << cg3_rc << std::endl;
                }
                                
                        }
                    }
                }
            }
        }
        exit(0);
    }


    std::vector<std::vector<uint8_t>> bfs1;
    std::vector<std::vector<uint8_t>> methylated_bfs1;
    // const variable 60 mil for bf
    //const size_t bfSize = 30000000000;
    // calculate size needed for a false positive rate of 0.1
    //const int bfSize = -1 * num_elements / log(0.1);

std::cerr << "making methylated kmers dataset" << std::endl;

//std::unordered_set<uint64_t> methylated_kmers_in_dataset;
//uint64_t max_size = max_size;
uint64_t max_size = 3000000000ULL;
btllib::BloomFilter prelim_ct_mers(max_size, 1);
btllib::BloomFilter prelim_ga_mers(max_size, 1);
btllib::BloomFilter clean_ct_mers(max_size, 1);
btllib::BloomFilter clean_ga_mers(max_size, 1);
btllib::BloomFilter methylated_kmers_in_dataset(max_size, 1);
//btllib::BloomFilter all_kmers_in_dataset(max_size, 1);
btllib::CountingBloomFilter8 error_kmer_ct(max_size, 3);
btllib::CountingBloomFilter8 error_kmer_ga(max_size, 3);



int num_lines_2 = 0;

/*std::ofstream methy_kmer_log("methylated_kmers_with_phred.tsv");
methy_kmer_log << "error_rate\n";  //
*/
std::vector<std::string> sample_names;
for (const auto& [prefix, pair] : pairs) {
    const auto& r1_file = pair.first;
    const auto& r2_file = pair.second;


    sample_names.push_back(prefix);
    std::cerr << num_lines_2 << std::endl;
    std::cerr << r1_file << std::endl;
    std::cerr << r2_file << std::endl;
    num_lines_2++;

    btllib::SeqReader reader(r1_file, btllib::SeqReader::Flag::SHORT_MODE);


#pragma omp parallel
    for (const auto record : reader) {
        if (record.seq.size() < k) {
            continue;
        }
        if (calculate_avg_phred(record.qual) < 20) {
            continue;
        }

        btllib::BsHashDirectional bh(record.seq, 3, k, "CT");
        btllib::BsHashDirectional bh_ga(record.seq, 3, k, "GA");
        while(bh.roll() && bh_ga.roll()) {
            
            size_t  j = bh.get_pos();
            std::string_view kmer{record.seq.data() + j, k};

            auto central_dimer = bh.center_dimer();
            if (central_dimer == "CG") {
                bool pass_quality = true;
                double total_error_prob = 0.0;

                for (size_t m = 0; m < k; ++m) {
                    int phred = static_cast<unsigned char>(record.qual[j + m]) - 33;
                    double error_prob = std::pow(10.0, -phred / 10.0);
                    total_error_prob += error_prob;
                }

                double avg_error_prob = total_error_prob / k;
                double avg_phred_score = -10.0 * std::log10(avg_error_prob);


                if (avg_phred_score < phred_threshold) {
                    pass_quality = false;
                }
                if (!pass_quality) continue;
                std::string_view orig_kmer{record.seq.data() + j, k};
size_t num_dimers = k / 2;
size_t center_dimer = num_dimers / 2;
size_t center_pos = center_dimer * 2;

std::string_view center_from_kmer{orig_kmer.data() + center_pos, 2};
assert(center_from_kmer == central_dimer);
                if (shannon_entropy(orig_kmer) < shannon && shannon_entropy_dimer(orig_kmer) < shannon2 && shannon_entropy_trimer(orig_kmer) < shannon3) {
                    continue;
                }

                error_kmer_ct.insert(bh.hashes());
                error_kmer_ga.insert(bh_ga.hashes());
              
            }
        }
    }
    if (!pair.second.empty()) {

        btllib::SeqReader reader2(r2_file, btllib::SeqReader::Flag::SHORT_MODE);


#pragma omp parallel
        for (const auto record : reader2) {
            if (record.seq.size() < k) {
                continue;
            }
            if (calculate_avg_phred(record.qual) < 20) {
                continue;
            }

            btllib::BsHashDirectional bh(record.seq, 3, k, "CT");
            btllib::BsHashDirectional bh_ga(record.seq, 3, k, "GA");

            while(bh.roll() && bh_ga.roll()) {
                size_t  j = bh.get_pos();

            auto central_dimer = bh.center_dimer();
            if (central_dimer == "CG") {
                    bool pass_quality = true;
                    double total_error_prob = 0.0;

                    for (size_t m = 0; m < k; ++m) {
                        int phred = static_cast<unsigned char>(record.qual[j + m]) - 33;
                        double error_prob = std::pow(10.0, -phred / 10.0);
                        total_error_prob += error_prob;
                    }

                    double avg_error_prob = total_error_prob / k;
                    double avg_phred_score = -10.0 * std::log10(avg_error_prob);

                    if (avg_phred_score < phred_threshold) {
                        pass_quality = false;
                    }
                    if (!pass_quality) continue;

                    std::string_view orig_kmer{record.seq.data() + j, k};
size_t num_dimers = k / 2;
size_t center_dimer = num_dimers / 2;
size_t center_pos = center_dimer * 2;

std::string_view center_from_kmer{orig_kmer.data() + center_pos, 2};
assert(center_from_kmer == central_dimer);
                    if (shannon_entropy(orig_kmer) < shannon && shannon_entropy_dimer(orig_kmer) < shannon2 && shannon_entropy_trimer(orig_kmer) < shannon3) {
                        continue;
                    }

                    error_kmer_ct.insert(bh.hashes());
                    error_kmer_ga.insert(bh_ga.hashes());
                    
                }
            }
        }
    }
}

std::cerr << "2nd pass" << std::endl;
for (const auto& [prefix, pair] : pairs) {
    const auto& r1_file = pair.first;
    const auto& r2_file = pair.second;


    std::cerr << num_lines_2 << std::endl;
    num_lines_2++;

btllib::SeqReader reader(r1_file, btllib::SeqReader::Flag::SHORT_MODE);


#pragma omp parallel
    for (const auto record : reader) {
        if (record.seq.size() < k) {
            continue;
        }
        if (calculate_avg_phred(record.qual) < 20) {
            continue;
        }

        btllib::BsHashDirectional bh(record.seq, 3, k, "CT");
        btllib::BsHashDirectional bh_ga(record.seq, 3, k, "GA");
        while(bh.roll() && bh_ga.roll()) {
            size_t  j = bh.get_pos();

            auto central_dimer = bh.center_dimer();
            if (central_dimer == "CG") {

                bool pass_quality = true;
                double total_error_prob = 0.0;

                for (size_t m = 0; m < k; ++m) {
                    int phred = static_cast<unsigned char>(record.qual[j + m]) - 33;
                    double error_prob = std::pow(10.0, -phred / 10.0);
                    total_error_prob += error_prob;
                }

                double avg_error_prob = total_error_prob / k;
                double avg_phred_score = -10.0 * std::log10(avg_error_prob);

                if (avg_phred_score < phred_threshold) {
                    pass_quality = false;
                }
                if (!pass_quality) continue;

                std::string_view orig_kmer{record.seq.data() + j, k};
size_t num_dimers = k / 2;
size_t center_dimer = num_dimers / 2;
size_t center_pos = center_dimer * 2;

std::string_view center_from_kmer{orig_kmer.data() + center_pos, 2};
assert(center_from_kmer == central_dimer);
                if (shannon_entropy(orig_kmer) < shannon && shannon_entropy_dimer(orig_kmer) < shannon2 && shannon_entropy_trimer(orig_kmer) < shannon3) {
                    continue;
                }
                if (error_kmer_ct.contains(bh.hashes()) > minKmer && error_kmer_ct.contains(bh.hashes()) < maxKmer) {
                    prelim_ct_mers.insert(bh.hashes());

                }
                if (error_kmer_ga.contains(bh_ga.hashes()) > minKmer && error_kmer_ga.contains(bh_ga.hashes()) < maxKmer) {
                   prelim_ga_mers.insert(bh_ga.hashes());
                }
            }
        }
    }

    if (!pair.second.empty()) {
        btllib::SeqReader reader2(r2_file, btllib::SeqReader::Flag::SHORT_MODE);


    #pragma omp parallel
        for (const auto record : reader2) {
            if (record.seq.size() < k) {
                continue;
            }
            if (calculate_avg_phred(record.qual) < 20) {
                continue;
            }
            btllib::BsHashDirectional bh(record.seq, 3, k, "CT");
            btllib::BsHashDirectional bh_ga(record.seq, 3, k, "GA");
            while(bh.roll() && bh_ga.roll()) {
                size_t j = bh.get_pos();
            auto central_dimer = bh.center_dimer();
            if (central_dimer == "CG") {
                    bool pass_quality = true;
                    double total_error_prob = 0.0;

                    for (size_t m = 0; m < k; ++m) {
                        int phred = static_cast<unsigned char>(record.qual[j + m]) - 33;
                        double error_prob = std::pow(10.0, -phred / 10.0);
                        total_error_prob += error_prob;
                    }

                    double avg_error_prob = total_error_prob / k;
                    double avg_phred_score = -10.0 * std::log10(avg_error_prob);

                    if (avg_phred_score < phred_threshold) {
                        pass_quality = false;
                    }
                    if (!pass_quality) continue;
                    std::string_view orig_kmer{record.seq.data() + j, k};
                    size_t num_dimers = k / 2;
                    size_t center_dimer = num_dimers / 2;
                    size_t center_pos = center_dimer * 2;

                    std::string_view center_from_kmer{orig_kmer.data() + center_pos, 2};
                    assert(center_from_kmer == central_dimer);
                    if (shannon_entropy(orig_kmer) < shannon && shannon_entropy_dimer(orig_kmer) < shannon2 && shannon_entropy_trimer(orig_kmer) < shannon3) {
                        continue;
                    }
                    if (error_kmer_ct.contains(bh.hashes()) > minKmer && error_kmer_ct.contains(bh.hashes()) < maxKmer) {
                        prelim_ct_mers.insert(bh.hashes());
                    }
                    if (error_kmer_ga.contains(bh_ga.hashes()) > minKmer && error_kmer_ga.contains(bh_ga.hashes()) < maxKmer) {
                        prelim_ga_mers.insert(bh_ga.hashes());
                    }
                }
            }
        }
    }
}


std::cerr << "3rd pass" << std::endl;
for (const auto& [prefix, pair] : pairs) {
    const auto& r1_file = pair.first;
    const auto& r2_file = pair.second;


    std::cerr << num_lines_2 << std::endl;
    num_lines_2++;

    btllib::SeqReader reader(r1_file, btllib::SeqReader::Flag::SHORT_MODE);


#pragma omp parallel
    for (const auto record : reader) {
        if (record.seq.size() < k) {
            continue;
        }
        if (calculate_avg_phred(record.qual) < 20) {
            continue;
        }

        btllib::BsHashDirectional bh(record.seq, 3, k, "CT");
        btllib::BsHashDirectional bh_ga(record.seq, 3, k, "GA");
        while(bh.roll() && bh_ga.roll()) {
            auto central_dimer = bh.center_dimer();
            if (central_dimer == "TG") {
                if (prelim_ct_mers.contains(bh.hashes())) {
                    clean_ct_mers.insert(bh.hashes());
                    methylated_kmers_in_dataset.insert(bh.hashes());
                }
            }
            if (central_dimer == "CA") {
                if (prelim_ga_mers.contains(bh_ga.hashes())) {
                   clean_ga_mers.insert(bh_ga.hashes());
                   methylated_kmers_in_dataset.insert(bh_ga.hashes());
                }
            }
        }
    }

    if (!pair.second.empty()) {
        btllib::SeqReader reader2(r2_file, btllib::SeqReader::Flag::SHORT_MODE);


    #pragma omp parallel
        for (const auto record : reader2) {
            if (record.seq.size() < k) {
                continue;
            }
            if (calculate_avg_phred(record.qual) < 20) {
                continue;
            }
            btllib::BsHashDirectional bh(record.seq, 3, k, "CT");
            btllib::BsHashDirectional bh_ga(record.seq, 3, k, "GA");
            while(bh.roll() && bh_ga.roll()) {
                auto central_dimer = bh.center_dimer();
                if (central_dimer == "TG") {
                    if (prelim_ct_mers.contains(bh.hashes())) {
                        clean_ct_mers.insert(bh.hashes());
                        methylated_kmers_in_dataset.insert(bh.hashes());

                    }
                }
                if (central_dimer == "CA") {
                    if (prelim_ga_mers.contains(bh_ga.hashes())) {
                        clean_ga_mers.insert(bh_ga.hashes());
                        methylated_kmers_in_dataset.insert(bh_ga.hashes());
                    }
                }
            }
        }
    }
}




std::unordered_map<uint64_t, size_t> hash_to_loc_map;

/*std::atomic<size_t> all_bfSize(0);

#pragma omp parallel for schedule(static)
for (uint64_t i = 0; i < max_size * 8; ++i) {
    if (all_kmers_in_dataset.contains({i})) {  // Bloom filter says 'possibly present'
        all_bfSize.fetch_add(1, std::memory_order_relaxed);  // unique dense index
    }
}

std::cerr << "all_bfSize: " <<  all_bfSize <<std::endl;*/

std::cerr << "Mapping Bloom Filter" << std::endl;
std::atomic<size_t> bfSize(0);

#pragma omp parallel for schedule(static)
for (uint64_t i = 0; i < max_size * 8; ++i) {
    if (methylated_kmers_in_dataset.contains({i})) {  // Bloom filter says 'possibly present'
        size_t index = bfSize.fetch_add(1, std::memory_order_relaxed);  // unique dense index
        #pragma omp critical
        {
            hash_to_loc_map[i] = index;
        }
    }
}
std::unordered_map<uint64_t, std::pair<double, double>> kmer_counts;
std::mutex kmer_mutex;
std::cerr << "bfSize: " <<  bfSize <<std::endl;
std::cerr << "making bit vector" << std::endl;
int num_lines = 0;
//omp_set_num_threads(1);
for (const auto& [prefix, pair] : pairs) {
    const auto& r1_file = pair.first;
    const auto& r2_file = pair.second;


    std::cerr << "reading line " << num_lines++ << std::endl;

    /*std::vector<btllib::SeqReader::Record> records;
    for (const auto& record : btllib::SeqReader(line1, btllib::SeqReader::Flag::SHORT_MODE | btllib::SeqReader::Flag::FOLD_CASE)) {
        records.push_back(record);
    }*/

    std::vector<uint8_t> final_bf(bfSize, 0);
    std::vector<uint8_t> final_meth_bf(bfSize, 0);

    btllib::SeqReader reader(r1_file, btllib::SeqReader::Flag::SHORT_MODE);
//two pass

#pragma omp parallel
    for (const auto record : reader) {
            if (record.seq.size() < k) {
                continue;
            }
        //std::cerr << "checking get all meth" <<  std::endl;
        std::vector<std::pair<uint64_t, bool>> all_kmers;
//#pragma omp critical
//{
        all_kmers =
            get_all_methylation_kmers(record.seq, k, clean_ct_mers, clean_ga_mers);//, shannon, shannon2, shannon3, record.qual, phred_threshold);
//}

        //std::cerr << "Done checking get all meth" <<  std::endl;

        for (const auto& kmer : all_kmers) {
            //std::cerr << "Mapping kmer to hash" <<  std::endl;
            size_t idx = hash_to_loc_map[kmer.first % (max_size * 8)];
            //std::cerr << "Mapped" <<  std::endl;
            final_bf[idx] = 1;
            if (kmer.second) {
                final_meth_bf[idx] = 1;
            }

            // Track methylation/unmethylation counts
            uint64_t key = kmer.first % (max_size * 8);  // hash modulus for frequency map

            {
                std::lock_guard<std::mutex> lock(kmer_mutex);
                auto& counts = kmer_counts[key];
                if (kmer.second) {
                    counts.second += 1.0;  // methylated
                } else {
                    counts.first += 1.0;   // unmethylated
                }
            }
        }
        //std::cerr << "Inserted checking get all meth" <<  std::endl;
    }

    if (!pair.second.empty()) {
        std::cerr << "read2" << std::endl;
        btllib::SeqReader reader2(r2_file, btllib::SeqReader::Flag::SHORT_MODE);
    //two pass

    #pragma omp parallel
        for (const auto record : reader2) {
            if (record.seq.size() < k) {
                continue;
            }
        std::vector<std::pair<uint64_t, bool>> all_kmers;
//#pragma omp critical
//{
        all_kmers =
                get_all_methylation_kmers(record.seq, k, clean_ct_mers, clean_ga_mers);//, shannon, shannon2, shannon3, record.qual, phred_threshold);
//}

            for (const auto& kmer : all_kmers) {
                size_t idx = hash_to_loc_map[kmer.first % (max_size * 8)];
                final_bf[idx] = 1;
                if (kmer.second) {
                    final_meth_bf[idx] = 1;
                }

                // Track methylation/unmethylation counts
                uint64_t key = kmer.first % (max_size * 8);  // hash modulus for frequency map

                {
                    std::lock_guard<std::mutex> lock(kmer_mutex);
                    auto& counts = kmer_counts[key];
                    if (kmer.second) {
                        counts.second += 1.0;  // methylated
                    } else {
                        counts.first += 1.0;   // unmethylated
                    }
                }
            }
        }
    }

    bfs1.emplace_back(std::move(final_bf));
    methylated_bfs1.emplace_back(std::move(final_meth_bf));
}

double total_occurrences = 0.0;

// Step 1: Sum total counts across all k-mers
for (const auto& [_, counts] : kmer_counts) {
    total_occurrences += counts.first + counts.second;
}

// Step 2: Replace each entry with TF-IDF
for (auto& [key, counts] : kmer_counts) {
    double total = counts.first + counts.second;
    if (total == 0.0) continue;

    double tf_unmeth = counts.first / total;
    double tf_meth = counts.second / total;

    // Use inverse TF as weights (to penalize frequent terms)
    double inv_tf_unmeth = (tf_unmeth > 0.0) ? 1.0 / tf_unmeth : 0.0;
    double inv_tf_meth   = (tf_meth   > 0.0) ? 1.0 / tf_meth   : 0.0;

    // Scale so the larger one is 1.0, and the other is relative to it
    double max_val = std::max(inv_tf_unmeth, inv_tf_meth);
    if (max_val > 0.0) {
        counts.first  = inv_tf_unmeth / max_val;
        counts.second = inv_tf_meth   / max_val;
    } else {
        counts.first = 0.0;
        counts.second = 0.0;
    }
}
// compute the Jaccard similarity between the two sets of k-mers in bfs1 against itself
// increase intersection only when both are 1

/*for (int i = 0; i < bfs1.size(); ++i) {
    for (int j = i + 1; j < bfs1.size(); ++j) {
        int intersection = 0;
        int union_size = 0;
        for (int k = 0; k < bfSize; ++k) {
            if (bfs1[i][k] == 1 && bfs1[j][k] == 1) {
                ++intersection;
            }
            if (bfs1[i][k] == 1 || bfs1[j][k] == 1) {
                ++union_size;
            }
        }
        double jaccard = static_cast<double>(intersection) / union_size;
        std::cout << "Jaccard similarity between " << i << " and " << j << ": " << jaccard << "\n";
    }
}*/

std::cerr << "calculating stats" << std::endl;
std::cerr << "using tf inverse" << std::endl;
// --- keep your includes and existing code ---

// After you finish filling bfs1 and methylated_bfs1 vectors and bfSize is set:


std::ofstream output_hamming(prefix + "hamming.tsv");
std::ofstream output_cosine(prefix + "cosine.tsv");
std::ofstream output_pearson(prefix + "pearson.tsv");

if (!output_hamming.is_open() || !output_cosine.is_open() || !output_pearson.is_open()) {
    std::cerr << "Unable to open output file(s)\n";
    return 1;
}

size_t numSamples = bfs1.size();
std::vector<std::atomic<uint32_t>> site_counts(bfSize);
for (size_t k = 0; k < bfSize; ++k) site_counts[k].store(0);

#pragma omp parallel for schedule(dynamic)
for (size_t i = 0; i < numSamples; ++i) {
    for (size_t k = 0; k < bfSize; ++k) {
        if (bfs1[i][k] == 1) {
            site_counts[k].fetch_add(1, std::memory_order_relaxed);
        }
    }
}

// Identify top 20% sites (Sorting is typically fast enough on main thread)
std::vector<size_t> indices(bfSize);
std::iota(indices.begin(), indices.end(), 0);
std::sort(indices.begin(), indices.end(), [&](size_t a, size_t b) {
    return site_counts[a].load() > site_counts[b].load();
});

size_t top20_limit = static_cast<size_t>(bfSize * 0.20);
std::vector<size_t> top_observed_indices(indices.begin(), indices.begin() + top20_limit);


struct SiteInfo {
    size_t index;
    double entropy;
};

std::vector<SiteInfo> site_results(top20_limit);

#pragma omp parallel for schedule(dynamic)
for (size_t idx = 0; idx < top20_limit; ++idx) {
    size_t k = top_observed_indices[idx];
    uint32_t count0 = 0, count1 = 0;

    for (size_t i = 0; i < numSamples; ++i) {
        if (bfs1[i][k] == 1) {
            (methylated_bfs1[i][k] == 1) ? count1++ : count0++;
        }
    }

    double total = count0 + count1;
    double entropy = 0.0;
    if (total > 0) {
        double p0 = count0 / total;
        double p1 = count1 / total;
        if (p0 > 0) entropy -= p0 * std::log2(p0);
        if (p1 > 0) entropy -= p1 * std::log2(p1);
    }
    site_results[idx] = {k, entropy};
}

// Filter to top 50,000 sites by entropy
std::sort(site_results.begin(), site_results.end(), [](const SiteInfo& a, const SiteInfo& b) {
    return a.entropy > b.entropy;
});

size_t final_count = std::min((size_t)50000, site_results.size());
std::vector<size_t> final_indices(final_count);
for (size_t i = 0; i < final_count; ++i) final_indices[i] = site_results[i].index;

export_matrices(bfs1, methylated_bfs1, final_indices, sample_names, dev);


    return 0;
}
