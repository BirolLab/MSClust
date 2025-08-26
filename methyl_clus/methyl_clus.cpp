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
#include "btllib/nthash.hpp"
#include <btllib/seq.hpp>
#include <btllib/seq_reader.hpp>

#include <boost/container_hash/hash.hpp>

#include "city.h"

size_t hash_combine( size_t lhs, size_t rhs ) {
  lhs ^= rhs + 0x9e3779b9 + (lhs << 6) + (lhs >> 2);
  return lhs;
}

uint64_t combine_triple(uint64_t h1, uint64_t h2, uint64_t h3) {
    std::size_t seed = 0;
    seed = hash_combine(seed, h1);
    seed = hash_combine(seed, h2);
    seed = hash_combine(seed, h3);
    return static_cast<uint64_t>(seed);
}

std::vector<std::pair<uint64_t, bool>> get_all_methylation_kmers(const std::string& seq, unsigned k, btllib::BloomFilter& methylated_kmers_in_dataset, bool dev) {
    btllib::NtHash itr(seq, 1, k);
    itr.roll();
    std::vector<std::pair<uint64_t, bool>> all_kmers_hash;
    for (size_t i = 0; i < seq.size() - k + 1; ++i) {
        bool is_methylated = false;
        // locate the k-mer where the center base is a CpG site encoded as a 1 or a CpG site encoded as a C and the next base is a G
        char meth_base = dev ? '1' : 'C';
        if (seq[i + k / 2] == meth_base || (seq[i + k / 2] == 'T' && seq[i + k / 2 + 1] == 'G')) {
            if (seq[i + k / 2] == meth_base) {
                is_methylated = true;
            } else {
                // check if the kmer is in the dataset
                //itr.sub({k/2 - 1}, {'T'});

                // Pick smaller (canonical form)
                if (!methylated_kmers_in_dataset.contains(itr.hashes())) {
                    itr.roll();
                    continue;
                }
            }
            itr.sub({k/2 - 1}, {'T'});
            all_kmers_hash.push_back(std::make_pair(itr.hashes()[0], is_methylated));
            if (is_methylated) {
                itr.sub({k/2 - 1}, {'C'}); 
            }           
        }
        itr.roll();
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
    program.add_argument("-k", "--kmer")
        // unsigned
        .scan<'u', unsigned>()
        .default_value(25u)
        .help("specify the kmer size (default is 25)");
    // add a development flag
    program.add_argument("-d", "--dev")
        .default_value(false)
        .implicit_value(true)
        .help("enable development mode (default is false)");

    program.add_argument("-b", "--debug")
        .default_value(false)
        .implicit_value(true)
        .help("enable debug message (default is false)");


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
    bool debug = program.get<bool>("-b");
    int numThreads = program.get<int>("-t");
    unsigned k = program.get<unsigned>("-k");
omp_set_num_threads(numThreads);

    std::cout << "Number of threads: " << numThreads << "\n";

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


    // const variable 60 mil for bf
    //const size_t bfSize = 3000000000;
    // calculate size needed for a false positive rate of 0.1
    //const int bfSize = -1 * num_elements / log(0.1);

std::cerr << "making methylated kmers dataset" << std::endl;

std::unordered_map<uint64_t, size_t> hash_to_loc_map;
std::atomic<size_t> bfSize(0);
std::unordered_map<uint64_t, size_t> trip_hash_to_loc_map;
std::atomic<size_t> trip_bfSize(0);

std::cerr << "bf size: " << bfSize << std::endl;
std::cerr << "trip bf size: " << trip_bfSize << std::endl;

//std::unordered_set<uint64_t> methylated_kmers_in_dataset;
btllib::BloomFilter methylated_kmers_in_dataset(30000000000ULL, 1);
btllib::BloomFilter triple_methylated_kmers_in_dataset(30000000000ULL, 1);




int num_lines_2 = 0;

for (const auto& line1 : lines1) {
    std::cerr << num_lines_2 << std::endl;
    num_lines_2++;

btllib::SeqReader reader(line1, btllib::SeqReader::Flag::SHORT_MODE);
//two pass

#pragma omp parallel
    for (const auto record : reader) {
        btllib::NtHash itr(record.seq, 1, k);
        itr.roll();
        //std::vector<uint64_t> list_of_meth_states;
        for (size_t j = 0; j + k <= record.seq.size(); ++j) {
            
            char meth_base = dev ? '1' : 'C';
            if (record.seq[j + k / 2] == meth_base && record.seq[j + k / 2 + 1] == 'G') {
                itr.sub({k/2 - 1}, {'T'});
                //list_of_meth_states.emplace_back(hashes[0]);
                methylated_kmers_in_dataset.insert(itr.hashes());
                itr.sub({k/2 - 1}, {'C'});
            }
            itr.roll();
            // insert into triple

        }
    }
}


num_lines_2 = 0;

for (const auto& line1 : lines1) {
    std::cerr << num_lines_2 << std::endl;
    num_lines_2++;

btllib::SeqReader reader(line1, btllib::SeqReader::Flag::SHORT_MODE);
//two pass

#pragma omp parallel
    for (const auto record : reader) {
        btllib::NtHash itr(record.seq, 1, k);
        itr.roll();

        std::vector<uint64_t> list_of_meth_states;
        for (size_t j = 0; j + k <= record.seq.size(); ++j) {
            
            char meth_base = dev ? '1' : 'C';
            if ((record.seq[j + k / 2] == meth_base || record.seq[j + k / 2] == 'T' )&& record.seq[j + k / 2 + 1] == 'G') {
                itr.sub({k/2 - 1}, {'T'});
                if (!methylated_kmers_in_dataset.contains(itr.hashes())) {
                    if (record.seq[j + k / 2] == meth_base ) {
                        itr.sub({k/2 - 1}, {'C'});
                    }
                    continue;
                }
                list_of_meth_states.emplace_back(itr.hashes()[0]);
                if (record.seq[j + k / 2] == meth_base ) {
                    itr.sub({k/2 - 1}, {'C'});
                }
            }
            itr.roll();
            // insert into triple

        }
        for (size_t i = 0; i + 2 < list_of_meth_states.size(); ++i) {
    uint64_t combined_hash = combine_triple(
        list_of_meth_states[i],
        list_of_meth_states[i + 1],
        list_of_meth_states[i + 2]
    );
    triple_methylated_kmers_in_dataset.insert({combined_hash});
}
    }
}



#pragma omp parallel for schedule(static)
for (uint64_t i = 0; i < 30000000000ULL * 8; ++i) {
    if (methylated_kmers_in_dataset.contains({i})) {  // Bloom filter says 'possibly present'
        size_t index = bfSize.fetch_add(1, std::memory_order_relaxed);  // unique dense index
        #pragma omp critical
        {
            hash_to_loc_map[i] = index;
        }
    }
}


#pragma omp parallel for schedule(static)
for (uint64_t i = 0; i < 30000000000ULL * 8; ++i) {
    if (triple_methylated_kmers_in_dataset.contains({i})) {  // Bloom filter says 'possibly present'
        size_t index = trip_bfSize.fetch_add(1, std::memory_order_relaxed);  // unique dense index
        #pragma omp critical
        {
            trip_hash_to_loc_map[i] = index;
        }
    }
}




/*std::vector<std::vector<uint8_t>> bfs1;
std::vector<std::vector<uint8_t>> methylated_bfs1;*/
std::vector<std::vector<uint8_t>> consolidated_single_bfs1;

std::vector<std::vector<uint8_t>> trip_bfs1;
std::vector<std::vector<uint8_t>> trip_methylated_bfs1;


std::cerr << "making bloom filter" << std::endl;
int num_lines = 0;

for (const auto& line1 : lines1) {
    std::cerr << "[DEBUG] Reading line " << num_lines++ << std::endl;

    /*std::vector<btllib::SeqReader::Record> records;
    try {
        if (debug) std::cerr << "[DEBUG] Opening SeqReader for: " << line1 << std::endl;
        for (const auto& record : btllib::SeqReader(line1, btllib::SeqReader::Flag::SHORT_MODE | btllib::SeqReader::Flag::FOLD_CASE)) {
            records.push_back(record);
        }
        if (debug) std::cerr << "[DEBUG] Total records read: " << records.size() << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "[ERROR] SeqReader failed: " << e.what() << std::endl;
        throw;
    }*/

    /*std::vector<uint8_t> final_bf(bfSize, 0);
    std::vector<uint8_t> final_meth_bf(bfSize, 0);*/
    std::vector<uint8_t> final_consolidated_bf(bfSize * 2, 0);
    std::vector<uint8_t> final_trip_bf(trip_bfSize, 0);
    std::vector<uint8_t> final_trip_meth_bf(
        trip_bfSize * 8, 0 
    );

    if (debug) std::cerr << "[DEBUG] Starting OpenMP parallel loop" << std::endl;
btllib::SeqReader reader(line1, btllib::SeqReader::Flag::SHORT_MODE);
//two pass

#pragma omp parallel
    for (const auto record : reader) {


        std::vector<std::pair<uint64_t, bool>> all_kmers =
            get_all_methylation_kmers(record.seq, k, methylated_kmers_in_dataset, dev);

        for (const auto& kmer : all_kmers) {
            size_t idx = hash_to_loc_map[kmer.first % (30000000000ULL * 8)];
            /*final_bf[idx] = 1;
            if (kmer.second) {
                final_meth_bf[idx] = 1;
            }*/
            final_consolidated_bf[idx * 2] = 1;
            if (kmer.second) {
                final_consolidated_bf[idx * 2 + 1] = 1;
            }
        }
        if (debug) std::cerr << "[DEBUG] Finished regular hash" << std::endl;

        for (size_t j = 0; j + 2 < all_kmers.size(); ++j) {
            uint64_t combined_hash = combine_triple(
                all_kmers[j].first,
                all_kmers[j + 1].first,
                all_kmers[j + 2].first
            );

            size_t idx = trip_hash_to_loc_map[combined_hash % (30000000000ULL * 8)];
            if (debug) {
                std::cerr << "combined_hash : " << combined_hash << std::endl;
                std::cerr << "combined_hash mod : " << combined_hash % (30000000000ULL * 8) << std::endl;
                std::cerr << "idx : " << idx << std::endl;
                std::cerr << "final_trip_bf size : " << trip_bfSize << std::endl;
                if (trip_hash_to_loc_map.count(combined_hash % (30000000000ULL * 8)) == 0) {
                    std::cerr << "Combined hash not found" << std::endl;
                } else {
                    std::cerr << "Combined hash found" << std::endl;
                }
            }
            final_trip_bf[idx] = 1;

            uint8_t pattern = (all_kmers[j].second   ? 4 : 0) |
                              (all_kmers[j + 1].second ? 2 : 0) |
                              (all_kmers[j + 2].second ? 1 : 0);

            final_trip_meth_bf[idx * 8 + pattern] = 1;
            if (debug) std::cerr << "[DEBUG] Finished triple hash" << std::endl;
        }
    }
    if (debug) std::cerr << "[DEBUG] Finished processing line" << std::endl;

    /*bfs1.emplace_back(std::move(final_bf));
    methylated_bfs1.emplace_back(std::move(final_meth_bf));*/
    consolidated_single_bfs1.emplace_back(std::move(final_consolidated_bf));
    trip_bfs1.emplace_back(std::move(final_trip_bf));
    trip_methylated_bfs1.emplace_back(std::move(final_trip_meth_bf));
    if (debug) std::cerr << "[DEBUG] Data pushed to BFS vectors" << std::endl;
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
for (size_t i = 0; i < consolidated_single_bfs1.size(); ++i) {
    for (size_t j = i + 1; j < consolidated_single_bfs1.size(); ++j) {
        int intersection = 0;
        int methylated_intersection = 0;
        double dot = 0.0;
        double sumAi2_cos = 0.0, sumAj2_cos = 0.0;
        double sumAiAj = 0.0, sumAi2 = 0.0, sumAj2 = 0.0;
        double shared_sites = 0;
        double sumA = 0.0, sumB = 0.0;

        for (size_t k = 0; k < bfSize; ++k) {
            if (consolidated_single_bfs1[i][k*2] == 1 && consolidated_single_bfs1[j][k*2] == 1) {
                ++intersection;
                uint8_t A = consolidated_single_bfs1[i][k * 2 + 1];
                uint8_t B = consolidated_single_bfs1[j][k * 2 + 1];

                if (A == B) ++methylated_intersection;

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
                if (consolidated_single_bfs1[i][k * 2] == 1 &&consolidated_single_bfs1[j][k * 2] == 1) {
                    double A = consolidated_single_bfs1[i][k * 2 + 1];
                    double B = consolidated_single_bfs1[j][k * 2 + 1];
                    double dA = A - meanA;
                    double dB = B - meanB;
                    sumAiAj += dA * dB;
                    sumAi2 += dA * dA;
                    sumAj2 += dB * dB;
                }
            }

            pearson_sim = (sumAi2 > 0 && sumAj2 > 0) ? sumAiAj / (sqrt(sumAi2) * sqrt(sumAj2)) : 0.0;
        }

        size_t pos = 0, pos2 = 0, pos3 = 0, pos4 = 0;
        if (dev) {
            pos = lines1[i].find("GSM");
            pos2 = lines1[i].find("_aligned_reads.fasta");
            pos3 = lines1[j].find("GSM");
            pos4 = lines1[j].find("_aligned_reads.fasta");
        } else {
            pos = lines1[i].find_last_of("/");
            pos3 = lines1[j].find_last_of("/");
            pos2 = lines1[i].find_first_of(".fq");
            pos4 = lines1[j].find_first_of(".fq");
            if (pos2 == std::string::npos) pos2 = lines1[i].find_first_of(".fastq");
            if (pos4 == std::string::npos) pos4 = lines1[j].find_first_of(".fastq");
        }

        std::string name1 = lines1[i].substr(pos, pos2 - pos);
        std::string name2 = lines1[j].substr(pos3, pos4 - pos3);

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


// --- Keep your existing bfs1 / methylated_bfs1 calculation loop above ---

std::cerr << "calculating triple jaccard" << std::endl;

std::ofstream output_trip_hamming(prefix + "trip_hamming.tsv");
std::ofstream output_trip_cosine(prefix + "trip_cosine.tsv");
std::ofstream output_trip_pearson(prefix + "trip_pearson.tsv");

if (!output_trip_hamming.is_open() || !output_trip_cosine.is_open() || !output_trip_pearson.is_open()) {
    std::cerr << "Unable to open triple output file(s)\n";
    return 1;
}

#pragma omp parallel for schedule(dynamic)
for (size_t i = 0; i < trip_bfs1.size(); ++i) {
    for (size_t j = i + 1; j < trip_bfs1.size(); ++j) {
        int intersection = 0;
        int methylated_intersection = 0;
        double dot = 0.0;
        double sumAi2_cos = 0.0, sumAj2_cos = 0.0;
        double sumAiAj = 0.0, sumAi2 = 0.0, sumAj2 = 0.0;
        double shared_sites = 0;
        double sumA = 0.0, sumB = 0.0;

        for (size_t k = 0; k < trip_bfSize; ++k) {
            if (trip_bfs1[i][k] == 1 && trip_bfs1[j][k] == 1) {
                ++intersection;

                // Compare across all 8 methylation patterns
                int match_count = 0;
                for (int p = 0; p < 8; ++p) {
                    uint8_t A = trip_methylated_bfs1[i][k * 8 +p];
                    uint8_t B = trip_methylated_bfs1[j][k * 8 +p];
                    if (A == B) ++match_count;

                    // Cosine
                    dot += A * B;
                    sumAi2_cos += A * A;
                    sumAj2_cos += B * B;

                    // Pearson sums
                    sumA += A;
                    sumB += B;
                    ++shared_sites;
                }

                methylated_intersection += match_count;
            }
        }

        double cosine_sim = (sumAi2_cos > 0 && sumAj2_cos > 0) ? dot / (sqrt(sumAi2_cos) * sqrt(sumAj2_cos)) : 0.0;
        double pearson_sim = 0.0;
        double jaccard = intersection > 0 ? static_cast<double>(methylated_intersection) / (intersection * 8) : 0.0;

        if (shared_sites > 0) {
            double meanA = sumA / shared_sites;
            double meanB = sumB / shared_sites;

            for (size_t k = 0; k < trip_bfSize; ++k) {
                if (trip_bfs1[i][k] == 1 && trip_bfs1[j][k] == 1) {
                    for (int p = 0; p < 8; ++p) {
                        double A = trip_methylated_bfs1[i][k * 8 + p];
                        double B = trip_methylated_bfs1[j][k * 8 + p];
                        double dA = A - meanA;
                        double dB = B - meanB;
                        sumAiAj += dA * dB;
                        sumAi2 += dA * dA;
                        sumAj2 += dB * dB;
                    }
                }
            }

            pearson_sim = (sumAi2 > 0 && sumAj2 > 0) ? sumAiAj / (sqrt(sumAi2) * sqrt(sumAj2)) : 0.0;
        }

        size_t pos = 0, pos2 = 0, pos3 = 0, pos4 = 0;
        if (dev) {
            pos = lines1[i].find("GSM");
            pos2 = lines1[i].find("_aligned_reads.fasta");
            pos3 = lines1[j].find("GSM");
            pos4 = lines1[j].find("_aligned_reads.fasta");
        } else {
            pos = lines1[i].find_last_of("/");
            pos3 = lines1[j].find_last_of("/");
            pos2 = lines1[i].find_first_of(".fq");
            pos4 = lines1[j].find_first_of(".fq");
            if (pos2 == std::string::npos) pos2 = lines1[i].find_first_of(".fastq");
            if (pos4 == std::string::npos) pos4 = lines1[j].find_first_of(".fastq");
        }

        std::string name1 = lines1[i].substr(pos, pos2 - pos);
        std::string name2 = lines1[j].substr(pos3, pos4 - pos3);

        #pragma omp critical
        {
            output_trip_hamming << name1 << "\t" << name2 << "\t" << jaccard << "\n";
            output_trip_cosine << name1 << "\t" << name2 << "\t" << cosine_sim << "\n";
            output_trip_pearson << name1 << "\t" << name2 << "\t" << pearson_sim << "\n";
        }
    }
}

output_trip_hamming.close();
output_trip_cosine.close();
output_trip_pearson.close();



    return 0;
}