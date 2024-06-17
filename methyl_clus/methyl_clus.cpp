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

std::vector<std::pair<uint64_t, bool>> get_all_methylation_kmers(const std::string& seq, unsigned k) {
    std::vector<std::pair<uint64_t, bool>> all_kmers_hash;
    for (size_t i = 0; i < seq.size() - k + 1; ++i) {
        bool is_methylated = false;
        // locate the k-mer where the center base is a CpG site encoded as a 1 or a CpG site encoded as a C and the next base is a G
        if (seq[i + k / 2] == '1' || (seq[i + k / 2] == 'C' && seq[i + k / 2 + 1] == 'G')) {
            if (seq[i + k / 2] == '1') {
                is_methylated = true;
            }
            std::string kmer = seq.substr(i, k);
            // covert all 1 to C
            std::replace(kmer.begin(), kmer.end(), '1', 'C');
            // upper case the kmer
            std::transform(kmer.begin(), kmer.end(), kmer.begin(), ::toupper);
            all_kmers_hash.push_back(std::make_pair(CityHash64(kmer.c_str(), k), is_methylated));            
        }

        
        
        
        
        
        /*if (seq[i + k / 2] == '1') {
            std::string kmer = seq.substr(i, k);
            // upper case the kmer
            std::transform(kmer.begin(), kmer.end(), kmer.begin(), ::toupper);
            all_kmers_hash.push_back(CityHash64(kmer.c_str(), k));            
        }*/
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

    try {
        program.parse_args(argc, argv);
    } catch (const std::runtime_error& err) {
        std::cerr << err.what() << std::endl;
        std::cerr << program;
        return 1;
    }

    std::string inputFile1 = program.get<std::string>("-i");
    std::string inputFile2 = program.get<std::string>("-j");
    int numThreads = program.get<int>("-t");
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

    std::vector<std::vector<bool>> bfs1;
    std::vector<std::vector<bool>> methylated_bfs1;
    // const variable 60 mil for bf
    const size_t bfSize = 6000000000;
    // calculate size needed for a false positive rate of 0.1
    //const int bfSize = -1 * num_elements / log(0.1);

    std::cerr << "making bloom filter" << std::endl;
    int num_lines = 0;

    for (const auto& line1 : lines1) {
        std::cerr << "reading line " << num_lines++ << std::endl;
        bfs1.emplace_back(bfSize, 0);
        methylated_bfs1.emplace_back(bfSize, 0);
#pragma omp parallel
        for (const auto& record : btllib::SeqReader(
                line1, btllib::SeqReader::Flag::LONG_MODE)) {
            std::vector<std::pair<uint64_t, bool>> all_kmers = get_all_methylation_kmers(record.seq, 25);
            for (const auto& kmer : all_kmers) {
                // create vector of a single value equal to kmer
                //std::vector<uint64_t> kmer_vec = {kmer};
                size_t idx = kmer.first % bfSize;
                bfs1.back()[idx] = 1;
                if (kmer.second) {
                    methylated_bfs1.back()[idx] = 1;
                }
            }
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
                size_t pos = lines1[i].find("GSM");
                // identify the postion of _aligned_reads.fasta
                size_t pos2 = lines1[i].find("_aligned_reads.fasta");

                size_t pos3 = lines1[j].find("GSM");
                size_t pos4 = lines1[j].find("_aligned_reads.fasta");


                std::string name1 = lines1[i].substr(pos, pos2 - pos);
                std::string name2 = lines1[j].substr(pos3, pos4 - pos3);
                output << name1 << "\t" << name2 << "\t" << jaccard << std::endl;
            }
        }
    }
    





    

    return 0;
}