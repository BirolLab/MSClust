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

float shannon_entropy(const std::string & s) {
    int counts[256] = {};
    for (unsigned char c : s) counts[c]++;
    float entropy = 0.0f;
    float length = static_cast<float>(s.size());
    for (int count : counts) {
        if (count == 0) continue;
        float p = count / length;
        entropy -= p * std::log2f(p);
    }
    return entropy;
}

float shannon_entropy_trimer(const std::string& s) {
    if (s.size() < 3) return 0.0f;

    std::unordered_map<std::string, int> counts;
    int total = 0;

    for (size_t i = 0; i + 2 < s.size(); ++i) {
        std::string kmer = s.substr(i, 3);
        counts[kmer]++;
        total++;
    }

    float entropy = 0.0f;
    for (const auto& [kmer, count] : counts) {
        float p = static_cast<float>(count) / total;
        entropy -= p * std::log2f(p);
    }

    return entropy;
}


float shannon_entropy_dimer(const std::string& s) {
    if (s.size() < 2) return 0.0f;

    std::unordered_map<std::string, int> counts;
    int total = 0;

    for (size_t i = 0; i + 1 < s.size(); ++i) {
        std::string kmer = s.substr(i, 2);
        counts[kmer]++;
        total++;
    }

    float entropy = 0.0f;
    for (const auto& [kmer, count] : counts) {
        float p = static_cast<float>(count) / total;
        entropy -= p * std::log2f(p);
    }

    return entropy;
}


/*std::vector<std::pair<uint64_t, bool>> get_all_methylation_kmers(const std::string& seq, unsigned k, btllib::BloomFilter& methylated_kmers_in_dataset, bool dev) {
    std::vector<std::pair<uint64_t, bool>> all_kmers_hash;
    for (size_t i = 0; i + k <= seq.size(); ++i) {
        bool is_methylated = false;
        // locate the k-mer where the center base is a CpG site encoded as a 1 or a CpG site encoded as a C and the next base is a G
        char meth_base = dev ? '1' : 'C';
        if (seq[i + k / 2] == meth_base || (seq[i + k / 2] == 'T' && seq[i + k / 2 + 1] == 'G')) {
            if (seq[i + k / 2] == meth_base) {
                is_methylated = true;
            } else {
                // check if the kmer is in the dataset
                std::string kmer = seq.substr(i, k);
                // covert all 1 to C
                std::replace(kmer.begin(), kmer.end(), meth_base, 'T');
                // upper case the kmer
                std::transform(kmer.begin(), kmer.end(), kmer.begin(), ::toupper);
                std::string reverse_kmer = btllib::get_reverse_complement(kmer);

                uint64_t hash_fwd = CityHash64WithSeed(kmer.c_str(), k,0);
                uint64_t hash_rev = CityHash64WithSeed(reverse_kmer.c_str(), k,0);

                // Pick smaller (canonical form)
                std::vector<uint64_t> hashes = {std::min(hash_fwd, hash_rev)};
                if (!methylated_kmers_in_dataset.contains(hashes)) {
                    continue;
                }
            }
            std::string kmer = seq.substr(i, k);
            // covert all 1 to T
            std::replace(kmer.begin(), kmer.end(), meth_base, 'T');
            // upper case the kmer
            std::transform(kmer.begin(), kmer.end(), kmer.begin(), ::toupper);
            std::string reverse_kmer = btllib::get_reverse_complement(kmer);

            uint64_t hash_fwd = CityHash64WithSeed(kmer.c_str(), k,0);
            uint64_t hash_rev = CityHash64WithSeed(reverse_kmer.c_str(), k,0);

            // Pick smaller (canonical form)
            std::vector<uint64_t> hashes = {std::min(hash_fwd, hash_rev)};
            if (!methylated_kmers_in_dataset.contains(hashes)) {
                continue;
            }
            all_kmers_hash.push_back(std::make_pair(hashes[0], is_methylated));            
        }
    }
    return all_kmers_hash;
}*/


std::vector<std::pair<uint64_t, bool>> get_all_methylation_kmers(
    const std::string& seq,
    unsigned k,
    btllib::BloomFilter& methylated_kmers_in_dataset,
    bool dev
) {
    dev = false;
    if (dev) {
        k = 1;
    }
    std::vector<std::pair<uint64_t, bool>> all_kmers_hash;

    const size_t seq_len = seq.size();
    if (seq_len < k) return all_kmers_hash;

    for (size_t i = 0; i + k <= seq_len; ++i) {
        const size_t mid = i + (k - 1) / 2;
        //if (mid + 1 >= seq_len) continue;

        // Only look at the two bases in the middle
        char base1 = seq[mid];
        char base2 = seq[mid + 1];

        bool is_methylated = false;
        std::vector<std::string> converted_kmers;

        if (base1 == 'C' && base2 == 'G') {        // CG
            is_methylated = true;
        } else if (base1 == 'T' && base2 == 'G') { // TG
            is_methylated = false;
        } else if (base1 == 'C' && base2 == 'A') { // CA
            is_methylated = false;
        } else {
            continue; // skip all other central dimers
        }

        // Check for N/n in the k-mer window
        bool has_N = false;
        for (size_t j = i; j < i + k; ++j) {
            char c = seq[j];
            if (c == 'N' || c == 'n') {
                has_N = true;
                break;
            }
        }
        if (has_N) continue;


        // Now that we know it's valid, generate the full k-mer
        std::string kmer = seq.substr(i, k);
        btllib::BsHashDirectional bh(kmer, 1, k, "CT"); 
        btllib::BsHashDirectional bh_ga(kmer, 1, k, "GA"); 

        // CG: do both C→T and G→A conversions
        /*if (base1 == 'C' && base2 == 'G') {
            std::string ct_kmer = kmer;
            std::replace(ct_kmer.begin(), ct_kmer.end(), dev ? '1' : 'C', 'T');
            std::transform(ct_kmer.begin(), ct_kmer.end(), ct_kmer.begin(), ::toupper);
            converted_kmers.push_back(ct_kmer);

            std::string ga_kmer = kmer;
            std::replace(ga_kmer.begin(), ga_kmer.end(), 'G', 'A');
            std::transform(ga_kmer.begin(), ga_kmer.end(), ga_kmer.begin(), ::toupper);
            converted_kmers.push_back(ga_kmer);
        }
        // TG: only C→T
        else if (base1 == 'T' && base2 == 'G') {
            std::string ct_kmer = kmer;
            std::replace(ct_kmer.begin(), ct_kmer.end(), dev ? '1' : 'C', 'T');
            std::transform(ct_kmer.begin(), ct_kmer.end(), ct_kmer.begin(), ::toupper);
            converted_kmers.push_back(ct_kmer);
        }
        // CA: only G→A
        else if (base1 == 'C' && base2 == 'A') {
            std::string ga_kmer = kmer;
            std::replace(ga_kmer.begin(), ga_kmer.end(), 'G', 'A');
            std::transform(ga_kmer.begin(), ga_kmer.end(), ga_kmer.begin(), ::toupper);
            converted_kmers.push_back(ga_kmer);
        }

        for (const auto& ck : converted_kmers) {
            std::string rev_ck = btllib::get_reverse_complement(ck);
            uint64_t hash_fwd = CityHash64WithSeed(ck.c_str(), k, 0);
            uint64_t hash_rev = CityHash64WithSeed(rev_ck.c_str(), k, 0);
            uint64_t canonical_hash = std::min(hash_fwd, hash_rev);

            if (!methylated_kmers_in_dataset.contains({canonical_hash})) continue;

            all_kmers_hash.emplace_back(canonical_hash, is_methylated);
        }*/
        if (base1 == 'C' && base2 == 'G') {
            bh.roll();
            bh_ga.roll();
            if (methylated_kmers_in_dataset.contains(bh.hashes())) { 
                all_kmers_hash.push_back(std::make_pair(bh.hashes()[0], is_methylated));
            } else if (methylated_kmers_in_dataset.contains(bh_ga.hashes())) { 
                all_kmers_hash.push_back(std::make_pair(bh_ga.hashes()[0], is_methylated));
            } else {
                continue;
            }
        }
        // TG: only C→T
        else if (base1 == 'T' && base2 == 'G') {
            bh.roll();
            if (!methylated_kmers_in_dataset.contains(bh.hashes())) { 
                continue; 
            } else {
                all_kmers_hash.push_back(std::make_pair(bh.hashes()[0], is_methylated));
            }
        }
        // CA: only G→A
        else if (base1 == 'C' && base2 == 'A') {
            bh_ga.roll();
            if (!methylated_kmers_in_dataset.contains(bh_ga.hashes())) { 
                continue; 
            } else {
                all_kmers_hash.push_back(std::make_pair(bh_ga.hashes()[0], is_methylated));
            }
        }
    }

    return all_kmers_hash;
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
/*for (const auto& [prefix, pair] : pairs) {
    std::cout << pair.first << std::endl;
    std::cout << pair.second << std::endl;
}

    exit(0);*/


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
btllib::BloomFilter methylated_kmers_in_dataset(max_size, 1);
//btllib::BloomFilter all_kmers_in_dataset(max_size, 1);
btllib::CountingBloomFilter8 error_kmer(max_size, 3);



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
    num_lines_2++;

    btllib::SeqReader reader(r1_file, btllib::SeqReader::Flag::SHORT_MODE);


#pragma omp parallel
    for (const auto record : reader) {
        for (size_t j = 0; j + k <= record.seq.size(); ++j) {
            char meth_base = dev ? '1' : 'C';
            if (record.seq[j + k / 2 - 1] == meth_base && record.seq[j + k / 2] == 'G') {
                // Check for N/n in the k-mer window
                bool has_N = false;
                for (size_t z = j; z< j + k; ++z) {
                    char c = record.seq[z];
                    if (c == 'N' || c == 'n') {
                        has_N = true;
                        break;
                    }
                }
                if (has_N) continue;


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
                /*std::string kmer = record.seq.substr(j, k);
                std::replace(kmer.begin(), kmer.end(), meth_base, 'T');
                std::transform(kmer.begin(), kmer.end(), kmer.begin(), ::toupper);
                if (shannon_entropy(kmer) < shannon && shannon_entropy_dimer(kmer) < shannon2 && shannon_entropy_trimer(kmer) < shannon3) {
                    continue;
                }
                std::string reverse_kmer = btllib::get_reverse_complement(kmer);
                std::vector<uint64_t> hashes;
                for (unsigned seed = 0; seed < 3; seed++)
                {
                    uint64_t hash_fwd = CityHash64WithSeed(kmer.c_str(), k, seed);
                    uint64_t hash_rev = CityHash64WithSeed(reverse_kmer.c_str(), k, seed);
                    hashes.push_back(std::min(hash_fwd, hash_rev));
                }*/
                std::string orig_kmer = record.seq.substr(j, k);
                if (shannon_entropy(orig_kmer) < shannon && shannon_entropy_dimer(orig_kmer) < shannon2 && shannon_entropy_trimer(orig_kmer) < shannon3) {
                    continue;
                }
                /*auto converted_kmers = generate_converted_kmers(orig_kmer, dev);

                std::vector<uint64_t> hashes;
                for (auto& ck : converted_kmers) {
                    std::string rev_ck = btllib::get_reverse_complement(ck);

                    for (unsigned seed = 0; seed < 3; seed++) {
                        uint64_t hash_fwd = CityHash64WithSeed(ck.c_str(), k, seed);
                        uint64_t hash_rev = CityHash64WithSeed(rev_ck.c_str(), k, seed);
                        hashes.push_back(std::min(hash_fwd, hash_rev)); // canonical
                    }
                }

                std::vector<uint64_t> ct_hashes(hashes.begin(), hashes.begin() + 3);
                std::vector<uint64_t> ga_hashes;
                if (hashes.size() > 3) {
                    ga_hashes.assign(hashes.begin() + 3, hashes.end());
                }


                // Pick smaller (canonical form)
                
                //all_kmers_in_dataset.insert(hashes);
                error_kmer.insert(ct_hashes);
                error_kmer.insert(ga_hashes);*/
                btllib::BsHashDirectional bh(orig_kmer, 3, k, "CT");
                btllib::BsHashDirectional bh_ga(orig_kmer, 3, k, "GA");
                bh.roll();
                bh_ga.roll();
                error_kmer.insert(bh.hashes());
                error_kmer.insert(bh_ga.hashes());
                
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
                    // Check for N/n in the k-mer window
                    bool has_N = false;
                for (size_t z = j; z< j + k; ++z) {
                    char c = record.seq[z];
                    if (c == 'N' || c == 'n') {
                        has_N = true;
                        break;
                    }
                }
                if (has_N) continue;

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
                    /*std::string kmer = record.seq.substr(j, k);
                    std::replace(kmer.begin(), kmer.end(), meth_base, 'T');
                    std::transform(kmer.begin(), kmer.end(), kmer.begin(), ::toupper);
                    if (shannon_entropy(kmer) < shannon && shannon_entropy_dimer(kmer) < shannon2 && shannon_entropy_trimer(kmer) < shannon3) {
                        continue;
                    }
                    std::string reverse_kmer = btllib::get_reverse_complement(kmer);
                    std::vector<uint64_t> hashes;
                    for (unsigned seed = 0; seed < 3; seed++)
                    {
                        uint64_t hash_fwd = CityHash64WithSeed(kmer.c_str(), k, seed);
                        uint64_t hash_rev = CityHash64WithSeed(reverse_kmer.c_str(), k, seed);
                        hashes.push_back(std::min(hash_fwd, hash_rev));
                    }*/
                    std::string orig_kmer = record.seq.substr(j, k);
                    if (shannon_entropy(orig_kmer) < shannon && shannon_entropy_dimer(orig_kmer) < shannon2 && shannon_entropy_trimer(orig_kmer) < shannon3) {
                        continue;
                    }
                    /*auto converted_kmers = generate_converted_kmers(orig_kmer, dev);

                    std::vector<uint64_t> hashes;
                    for (auto& ck : converted_kmers) {
                        std::string rev_ck = btllib::get_reverse_complement(ck);

                        for (unsigned seed = 0; seed < 3; seed++) {
                            uint64_t hash_fwd = CityHash64WithSeed(ck.c_str(), k, seed);
                            uint64_t hash_rev = CityHash64WithSeed(rev_ck.c_str(), k, seed);
                            hashes.push_back(std::min(hash_fwd, hash_rev)); // canonical
                        }
                    }
                    std::vector<uint64_t> ct_hashes(hashes.begin(), hashes.begin() + 3);
                    std::vector<uint64_t> ga_hashes;
                    if (hashes.size() > 3) {
                        ga_hashes.assign(hashes.begin() + 3, hashes.end());
                    }


                    // Pick smaller (canonical form)
                    
                    //all_kmers_in_dataset.insert(hashes);
                    error_kmer.insert(ct_hashes);
                    error_kmer.insert(ga_hashes);*/
                    btllib::BsHashDirectional bh(orig_kmer, 3, k, "CT");
                    btllib::BsHashDirectional bh_ga(orig_kmer, 3, k, "GA");
                    bh.roll();
                    bh_ga.roll();
                    error_kmer.insert(bh.hashes());
                    error_kmer.insert(bh_ga.hashes());
                    
                }
            }
        }
    }
}


for (const auto& [prefix, pair] : pairs) {
    const auto& r1_file = pair.first;
    const auto& r2_file = pair.second;


    std::cerr << num_lines_2 << std::endl;
    num_lines_2++;

btllib::SeqReader reader(r1_file, btllib::SeqReader::Flag::SHORT_MODE);


#pragma omp parallel
    for (const auto record : reader) {
        for (size_t j = 0; j + k <= record.seq.size(); ++j) {
            char meth_base = dev ? '1' : 'C';
            if (record.seq[j + k / 2 - 1] == meth_base && record.seq[j + k / 2] == 'G') {
                // Check for N/n in the k-mer window
                bool has_N = false;
                for (size_t z = j; z< j + k; ++z) {
                    char c = record.seq[z];
                    if (c == 'N' || c == 'n') {
                        has_N = true;
                        break;
                    }
                }
                if (has_N) continue;

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
                /*std::string kmer = record.seq.substr(j, k);
                std::replace(kmer.begin(), kmer.end(), meth_base, 'T');
                std::transform(kmer.begin(), kmer.end(), kmer.begin(), ::toupper);
                if (shannon_entropy(kmer) < shannon && shannon_entropy_dimer(kmer) < shannon2 && shannon_entropy_trimer(kmer) < shannon3) {
                    continue;
                }
                std::string reverse_kmer = btllib::get_reverse_complement(kmer);
                std::vector<uint64_t> hashes;
                for (unsigned seed = 0; seed < 3; seed++)
                {
                    uint64_t hash_fwd = CityHash64WithSeed(kmer.c_str(), k, seed);
                    uint64_t hash_rev = CityHash64WithSeed(reverse_kmer.c_str(), k, seed);
                    hashes.push_back(std::min(hash_fwd, hash_rev));
                }*/
                std::string orig_kmer = record.seq.substr(j, k);
                if (shannon_entropy(orig_kmer) < shannon && shannon_entropy_dimer(orig_kmer) < shannon2 && shannon_entropy_trimer(orig_kmer) < shannon3) {
                    continue;
                }
                /*auto converted_kmers = generate_converted_kmers(orig_kmer, dev);

                std::vector<uint64_t> hashes;
                for (auto& ck : converted_kmers) {
                    std::string rev_ck = btllib::get_reverse_complement(ck);

                    for (unsigned seed = 0; seed < 3; seed++) {
                        uint64_t hash_fwd = CityHash64WithSeed(ck.c_str(), k, seed);
                        uint64_t hash_rev = CityHash64WithSeed(rev_ck.c_str(), k, seed);
                        hashes.push_back(std::min(hash_fwd, hash_rev)); // canonical
                    }
                }
                std::vector<uint64_t> ct_hashes(hashes.begin(), hashes.begin() + 3);
                std::vector<uint64_t> ga_hashes;
                if (hashes.size() > 3) {
                    ga_hashes.assign(hashes.begin() + 3, hashes.end());
                }


                if (error_kmer.contains(ct_hashes) > minKmer && error_kmer.contains(ct_hashes) < maxKmer) {
                    methylated_kmers_in_dataset.insert(ct_hashes);
                    // optional logging code
                }



                if (error_kmer.contains(ga_hashes) > minKmer && error_kmer.contains(ga_hashes) < maxKmer) {
                    methylated_kmers_in_dataset.insert(ga_hashes);
                    // optional logging code
                }*/
                btllib::BsHashDirectional bh(orig_kmer, 3, k, "CT");
                btllib::BsHashDirectional bh_ga(orig_kmer, 3, k, "GA");
                bh.roll();
                bh_ga.roll();
                if (error_kmer.contains(bh.hashes()) > minKmer && error_kmer.contains(bh.hashes()) < maxKmer) {
                    methylated_kmers_in_dataset.insert(bh.hashes());
                    // optional logging code
                }
                if (error_kmer.contains(bh_ga.hashes()) > minKmer && error_kmer.contains(bh_ga.hashes()) < maxKmer) {
                    methylated_kmers_in_dataset.insert(bh_ga.hashes());
                    // optional logging code
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
                    // Check for N/n in the k-mer window
                    bool has_N = false;
                for (size_t z = j; z< j + k; ++z) {
                    char c = record.seq[z];
                    if (c == 'N' || c == 'n') {
                        has_N = true;
                        break;
                    }
                }
                if (has_N) continue;

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
                    /*std::string kmer = record.seq.substr(j, k);
                    std::replace(kmer.begin(), kmer.end(), meth_base, 'T');
                    std::transform(kmer.begin(), kmer.end(), kmer.begin(), ::toupper);
                    if (shannon_entropy(kmer) < shannon && shannon_entropy_dimer(kmer) < shannon2 && shannon_entropy_trimer(kmer) < shannon3) {
                        continue;
                    }
                    std::string reverse_kmer = btllib::get_reverse_complement(kmer);
                    std::vector<uint64_t> hashes;
                    for (unsigned seed = 0; seed < 3; seed++)
                    {
                        uint64_t hash_fwd = CityHash64WithSeed(kmer.c_str(), k, seed);
                        uint64_t hash_rev = CityHash64WithSeed(reverse_kmer.c_str(), k, seed);
                        hashes.push_back(std::min(hash_fwd, hash_rev));
                    }

                    if (error_kmer.contains(hashes) > minKmer && error_kmer.contains(hashes) < maxKmer) {
                        methylated_kmers_in_dataset.insert(hashes);
                        double error_sum = 0.0;
                        for (size_t q = j; q < j + k; ++q) {
                            int phred_score = record.qual[q] - 33; // ASCII to Phred
                            error_sum += std::pow(10.0, -phred_score / 10.0);
                        }
                        double avg_error_rate = error_sum / k * 100;

                        methy_kmer_log << avg_error_rate << '\n';
                    }*/
                    std::string orig_kmer = record.seq.substr(j, k);
                    if (shannon_entropy(orig_kmer) < shannon && shannon_entropy_dimer(orig_kmer) < shannon2 && shannon_entropy_trimer(orig_kmer) < shannon3) {
                        continue;
                    }
                    /*auto converted_kmers = generate_converted_kmers(orig_kmer, dev);

                    std::vector<uint64_t> hashes;
                    for (auto& ck : converted_kmers) {
                        std::string rev_ck = btllib::get_reverse_complement(ck);

                        for (unsigned seed = 0; seed < 3; seed++) {
                            uint64_t hash_fwd = CityHash64WithSeed(ck.c_str(), k, seed);
                            uint64_t hash_rev = CityHash64WithSeed(rev_ck.c_str(), k, seed);
                            hashes.push_back(std::min(hash_fwd, hash_rev)); // canonical
                        }
                    }


                    std::vector<uint64_t> ct_hashes(hashes.begin(), hashes.begin() + 3);
                    std::vector<uint64_t> ga_hashes;
                    if (hashes.size() > 3) {
                        ga_hashes.assign(hashes.begin() + 3, hashes.end());
                    }

                    if (error_kmer.contains(ct_hashes) > minKmer && error_kmer.contains(ct_hashes) < maxKmer) {
                        methylated_kmers_in_dataset.insert(ct_hashes);
                        // optional logging code
                    }



                    if (error_kmer.contains(ga_hashes) > minKmer && error_kmer.contains(ga_hashes) < maxKmer) {
                        methylated_kmers_in_dataset.insert(ga_hashes);
                        // optional logging code
                    }*/
                    btllib::BsHashDirectional bh(orig_kmer, 3, k, "CT");
                    btllib::BsHashDirectional bh_ga(orig_kmer, 3, k, "GA");
                    bh.roll();
                    bh_ga.roll();
                    if (error_kmer.contains(bh.hashes()) > minKmer && error_kmer.contains(bh.hashes()) < maxKmer) {
                        methylated_kmers_in_dataset.insert(bh.hashes());
                        // optional logging code
                    }
                    if (error_kmer.contains(bh_ga.hashes()) > minKmer && error_kmer.contains(bh_ga.hashes()) < maxKmer) {
                        methylated_kmers_in_dataset.insert(bh_ga.hashes());
                        // optional logging code
                    }
                }
            }
        }
    }
}


std::unordered_map<uint64_t, size_t> hash_to_loc_map;


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
std::cerr << "making bloom filter" << std::endl;
int num_lines = 0;

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

        std::vector<std::pair<uint64_t, bool>> all_kmers =
            get_all_methylation_kmers(record.seq, k, methylated_kmers_in_dataset, dev);

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
    if (!pair.second.empty()) {
    btllib::SeqReader reader2(r2_file, btllib::SeqReader::Flag::SHORT_MODE);
//two pass

    #pragma omp parallel
        for (const auto record : reader2) {

            std::vector<std::pair<uint64_t, bool>> all_kmers =
                get_all_methylation_kmers(record.seq, k, methylated_kmers_in_dataset, dev);

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

std::cerr << "calculating jaccard" << std::endl;
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

#pragma omp parallel for schedule(dynamic)
for (size_t i = 0; i < bfs1.size(); ++i) {
    for (size_t j = i + 1; j < bfs1.size(); ++j) {
        double intersection = 0;
        double methylated_intersection = 0;
        double dot = 0.0;
        double sumAi2_cos = 0.0, sumAj2_cos = 0.0;
        double sumAiAj = 0.0, sumAi2 = 0.0, sumAj2 = 0.0;
        double shared_sites = 0;
        double sumA = 0.0, sumB = 0.0;

        for (size_t k = 0; k < bfSize; ++k) {
            if (bfs1[i][k] == 1 && bfs1[j][k] == 1) {
                ++intersection;
                uint8_t A = methylated_bfs1[i][k];
                uint8_t B = methylated_bfs1[j][k];

                if (A == B) {
                    auto it = kmer_counts.find(k);
                    if (it != kmer_counts.end()) {
                        if (A == 0)
                            methylated_intersection += it->second.first;   // unmethylated TF-IDF
                        else
                            methylated_intersection += it->second.second;  // methylated TF-IDF
                    }
                }

                // For Cosine
                dot += A * B;
                sumAi2_cos += A * A;
                sumAj2_cos += B * B;

                // For Pearson
                sumA += A;
                sumB += B;
                ++shared_sites;
            }
        }

        double cosine_sim = (sumAi2_cos > 0 && sumAj2_cos > 0) ? dot / (sqrt(sumAi2_cos) * sqrt(sumAj2_cos)) : 0.0;
        double pearson_sim = 0.0;
        double jaccard = intersection > 0 ? static_cast<double>(methylated_intersection) / intersection : 0.0;

        if (shared_sites > 0) {
            double meanA = sumA / shared_sites;
            double meanB = sumB / shared_sites;

            for (size_t k = 0; k < bfSize; ++k) {
                if (bfs1[i][k] == 1 && bfs1[j][k] == 1) {
                    double A = methylated_bfs1[i][k];
                    double B = methylated_bfs1[j][k];
                    double dA = A - meanA;
                    double dB = B - meanB;
                    sumAiAj += dA * dB;
                    sumAi2 += dA * dA;
                    sumAj2 += dB * dB;
                }
            }

            pearson_sim = (sumAi2 > 0 && sumAj2 > 0) ? sumAiAj / (sqrt(sumAi2) * sqrt(sumAj2)) : 0.0;
        }

        /*size_t pos = 0, pos2 = 0, pos3 = 0, pos4 = 0;
        if (dev) {
            pos = sample_names[i].find("GSM");
            pos2 = sample_names[i].find("_aligned_reads.fasta");
            pos3 = sample_names[j].find("GSM");
            pos4 = sample_names[j].find("_aligned_reads.fasta");
        } else {
            pos = sample_names[i].find_last_of("/");
            pos3 = sample_names[j].find_last_of("/");
            pos2 = sample_names[i].find_first_of(".fq");
            pos4 = sample_names[j].find_first_of(".fq");
            if (pos2 == std::string::npos) pos2 = sample_names[i].find_first_of(".fastq");
            if (pos4 == std::string::npos) pos4 = sample_names[j].find_first_of(".fastq");
        }*/

        std::string name1 = sample_names[i];
        std::string name2 = sample_names[j];
        #pragma omp critical
        {
            output_hamming << name1 << "\t" << name2 << "\t" << jaccard << "\n";
            output_cosine << name1 << "\t" << name2 << "\t" << cosine_sim << "\n";
            output_pearson << name1 << "\t" << name2 << "\t" << pearson_sim << "\n";
        }
    }
}


output_hamming.close();
output_cosine.close();
output_pearson.close();





    return 0;
}