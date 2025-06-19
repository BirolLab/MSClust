#include <algorithm>
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
#include <btllib/seq.hpp>
#include <btllib/seq_reader.hpp>

#include "city.h"

std::vector<std::pair<uint64_t, bool>> get_all_methylation_kmers(const std::string& seq, unsigned k, std::unordered_set<uint64_t>& methylated_kmers_in_dataset, bool dev) {
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
                std::string kmer = seq.substr(i, k);
                // covert all 1 to C
                std::replace(kmer.begin(), kmer.end(), meth_base, 'T');
                // upper case the kmer
                std::transform(kmer.begin(), kmer.end(), kmer.begin(), ::toupper);
                if (methylated_kmers_in_dataset.find(CityHash64(kmer.c_str(), k)) == methylated_kmers_in_dataset.end()) {
                    continue;
                }
            }
            std::string kmer = seq.substr(i, k);
            // covert all 1 to T
            std::replace(kmer.begin(), kmer.end(), meth_base, 'T');
            // upper case the kmer
            std::transform(kmer.begin(), kmer.end(), kmer.begin(), ::toupper);
            all_kmers_hash.push_back(std::make_pair(CityHash64(kmer.c_str(), k), is_methylated));            
        }
    }
    return all_kmers_hash;
}

int main(int argc, char* argv[]) {
    argparse::ArgumentParser program("example");

    program.add_argument("-i")
        .required()
        .help("specify the first input file");

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


    try {
        program.parse_args(argc, argv);
    } catch (const std::runtime_error& err) {
        std::cerr << err.what() << std::endl;
        std::cerr << program;
        return 1;
    }

    std::string inputFile1 = program.get<std::string>("-i");
    std::string inputFile2 = program.get<std::string>("-j");
    bool dev = program.get<bool>("-d");
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

    std::vector<std::vector<uint8_t>> bfs1;
    std::vector<std::vector<uint8_t>> methylated_bfs1;
    // const variable 60 mil for bf
    //const size_t bfSize = 3000000000;
    // calculate size needed for a false positive rate of 0.1
    //const int bfSize = -1 * num_elements / log(0.1);

std::cerr << "making methylated kmers dataset" << std::endl;

std::unordered_set<uint64_t> methylated_kmers_in_dataset;




for (const auto& line1 : lines1) {
    // First, collect all records into a vector to allow indexing
    std::vector<btllib::SeqReader::Record> records;
    for (const auto& record : btllib::SeqReader(line1, btllib::SeqReader::Flag::SHORT_MODE | btllib::SeqReader::Flag::FOLD_CASE)) {
        records.push_back(record);
    }

    // Parallelize over the indexed records
    std::vector<std::unordered_set<uint64_t>> thread_local_sets(omp_get_max_threads());

    #pragma omp parallel for schedule(dynamic)
    for (size_t i = 0; i < records.size(); ++i) {
        const auto& record = records[i];
        int tid = omp_get_thread_num();
        auto& local_set = thread_local_sets[tid];

        for (size_t j = 0; j + k <= record.seq.size(); ++j) {
            char meth_base = dev ? '1' : 'C';
            if (record.seq[j + k / 2] == meth_base) {
                std::string kmer = record.seq.substr(j, k);
                std::replace(kmer.begin(), kmer.end(), meth_base, 'T');
                std::transform(kmer.begin(), kmer.end(), kmer.begin(), ::toupper);
                local_set.insert(CityHash64(kmer.c_str(), k));
            }
        }
    }

    // Merge thread-local sets into the global set
    for (const auto& local_set : thread_local_sets) {
        methylated_kmers_in_dataset.insert(local_set.begin(), local_set.end());
    }
}

std::unordered_map<uint64_t, size_t> hash_to_loc_map;

size_t bfSize = 0;
for (const auto& hash : methylated_kmers_in_dataset) {
    hash_to_loc_map[hash] = bfSize++;
}

std::cerr << "making bloom filter" << std::endl;
int num_lines = 0;

for (const auto& line1 : lines1) {
    std::cerr << "reading line " << num_lines++ << std::endl;

    std::vector<btllib::SeqReader::Record> records;
    for (const auto& record : btllib::SeqReader(line1, btllib::SeqReader::Flag::SHORT_MODE | btllib::SeqReader::Flag::FOLD_CASE)) {
        records.push_back(record);
    }

    std::vector<uint8_t> final_bf(bfSize, 0);
    std::vector<uint8_t> final_meth_bf(bfSize, 0);

    #pragma omp parallel for
    for (size_t i = 0; i < records.size(); ++i) {
        auto& record = records[i];
        std::vector<std::pair<uint64_t, bool>> all_kmers =
            get_all_methylation_kmers(record.seq, k, methylated_kmers_in_dataset, dev);

        for (const auto& kmer : all_kmers) {
            size_t idx = hash_to_loc_map[kmer.first];
            final_bf[idx] = 1;
            if (kmer.second) {
                final_meth_bf[idx] = 1;
            }
        }
    }

    bfs1.emplace_back(std::move(final_bf));
    methylated_bfs1.emplace_back(std::move(final_meth_bf));
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

// make an output file
    std::ofstream output("output.txt");
    if (!output.is_open()) {
        std::cerr << "Unable to open output file\n";
        return 1;
    }

    // calculate the Jaccard similarity between the two sets of k-mers in bfs1
    // against itself

#pragma omp parallel for
    for (size_t i = 0; i < bfs1.size(); ++i) {
        for (size_t j = i + 1; j < bfs1.size(); ++j) {
            // intersect both bfs1 and methylated_bfs1 and then divide the intersection of methylated_bfs1 by the intersection of bfs1
            int intersection = 0;
            // no union size needed since we are only interested in the intersection of methylated_bfs1
            for (size_t k = 0; k < bfSize; ++k) {
                if (bfs1[i][k] == 1 && bfs1[j][k] == 1) {
                    ++intersection;
                }
            }
            int methylated_intersection = 0;
            for (size_t k = 0; k < bfSize; ++k) {
                if ((bfs1[i][k] == 1 && bfs1[j][k] == 1 && methylated_bfs1[i][k] == methylated_bfs1[j][k] ) ) {
                    ++methylated_intersection;
                }
            }
            double jaccard = static_cast<double>(methylated_intersection) / intersection;
#pragma omp critical
            {
                // output the Jaccard similarity to the output file using names of the file
                // file name is in the format of /projects/ABySS_research/grant_prep/R21_single-cell-methylation/data/simulated_genomes/sim_haplotype/GSM5652179_Aorta-Endothel-Z00000422.hg38/11/_aligned_reads.fasta
                // strip /projects/ABySS_research/grant_prep/R21_single-cell-methylation/data/simulated_genomes/sim_haplotype/ and /_aligned_reads.fasta
                
                // identify postition of GSM in string
                size_t pos = 0;
                size_t pos2 = 0;
                size_t pos3 = 0;
                size_t pos4 = 0;
                if (dev) {
                    pos = lines1[i].find("GSM");
                    // identify the postion of _aligned_reads.fasta
                    pos2 = lines1[i].find("_aligned_reads.fasta");

                    pos3 = lines1[j].find("GSM");
                    pos4 = lines1[j].find("_aligned_reads.fasta");
                }
                else {
                    // get the base name of the full path
                    pos = lines1[i].find_last_of("/");
                    pos3 = lines1[j].find_last_of("/");
                    // find fastq or fq in the name and use that as anchor pos
                    pos2 = lines1[i].find_first_of(".fq");
                    pos4 = lines1[j].find_first_of(".fq");
                    if (pos2 == std::string::npos) {
                        pos2 = lines1[i].find_first_of(".fastq");
                    }
                    if (pos4 == std::string::npos) {
                        pos4 = lines1[j].find_first_of(".fastq");
                    }


                }


                    std::string name1 = lines1[i].substr(pos, pos2 - pos);
                    std::string name2 = lines1[j].substr(pos3, pos4 - pos3 );
                output << name1 << "\t" << name2 << "\t" << jaccard << std::endl;
            }
        }
    }
    





    

    return 0;
}