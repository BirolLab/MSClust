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

#include <cmath>
#include <boost/math/distributions/normal.hpp>

#include "city.h"
#include <numeric>
#include <boost/math/distributions/students_t.hpp>

#include <iostream>
#include <cmath>

// Function to calculate the Welch t-statistic
double welchTTest(double mean1, double mean2, double stdErr1, double stdErr2) {
    return (mean1 - mean2) / std::sqrt((stdErr1 * stdErr1) + (stdErr2 * stdErr2));
}

// Function to calculate the degrees of freedom for the Welch t-test
double welchDegreesOfFreedom(double stdErr1, double stdErr2, int n1, int n2) {
    //std::cerr << "n1: " << n1 << std::endl;
    //std::cerr << "n2: " << n2 << std::endl;
    if (stdErr1 <= 0 || stdErr2 <= 0) {
        std::cerr << "Standard error cannot be zero or negative" << std::endl;
        return std::numeric_limits<double>::quiet_NaN();
    }
    if (n1 <= 1 || n2 <= 1) {
        std::cerr << "Sample sizes must be greater than 1" << std::endl;
        return std::numeric_limits<double>::quiet_NaN();
    }
    double se1Squared = stdErr1 * stdErr1;
    double se2Squared = stdErr2 * stdErr2;
    double numerator = (se1Squared + se2Squared) * (se1Squared + se2Squared);
    double denominator = (se1Squared * se1Squared) / (n1 - 1) + (se2Squared * se2Squared) / (n2 - 1);
    //std::cerr << "numerator: " << numerator << std::endl;
    //std::cerr << "denominator: " << denominator << std::endl;
    //std::cerr << se1Squared << std::endl;
    //std::cerr << se2Squared << std::endl;
    return numerator / denominator;
}

// Function to calculate the p-value from the t-statistic and degrees of freedom
double calculatePValue(double tStatistic, double df, bool isLowerTail = false) {
  // Use the appropriate tail based on the flag
  boost::math::students_t dist(df);
  double pValue;
  if (isLowerTail) {
    pValue = boost::math::cdf(dist, tStatistic);
  } else {
    pValue = boost::math::cdf(boost::math::complement(dist, tStatistic));
  }
  return pValue;
}

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
        .scan<'i', size_t>()
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
    size_t numThreads = program.get<size_t>("-t");
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
    //const size_t bfSize = -1 * num_elements / log(0.1);

    std::cerr << "making bloom filter" << std::endl;
    size_t num_lines = 0;

    std::unordered_map<size_t, std::string> group_to_name;
    std::unordered_map<std::string, size_t> name_to_group;
    std::string old_name = "";

    for (const auto& line1 : lines1) {
        num_lines++;

        // identify postition of GSM in string
        size_t pos = line1.find("GSM");
        // identify the postion of _aligned_reads.fasta
        size_t pos2 = line1.find("-Z000");
        std::string name1 = line1.substr(pos, pos2 - pos);
        size_t pos3 = name1.find("_");
        name1 = name1.substr(pos3 + 1);

        

        if (old_name != name1) {
            bfs1.emplace_back(bfSize, 0);
            methylated_bfs1.emplace_back(bfSize, 0);
            old_name = name1;
        }

        
        
        // if the name is not in the map, add it
        if (name_to_group.find(name1) == name_to_group.end()) {
            std::cerr << "name1: " << name1 << std::endl;
            size_t group = bfs1.size();
            std::cerr << "group: " << group << std::endl;
            group_to_name[group] = name1;
            name_to_group[name1] = group;
        }
        btllib::SeqReader reader(
                line1, btllib::SeqReader::Flag::LONG_MODE);
#pragma omp parallel
        for (const auto& record : reader) {
            auto all_kmers = get_all_methylation_kmers(record.seq, 33);
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


std::cerr << "calculating jaccard" << std::endl;

// make an output file
    std::ofstream output("classification.txt");
    if (!output.is_open()) {
        std::cerr << "Unable to open output file\n";
        return 1;
    }
    std::vector<std::ofstream> output_fasta;
    // open 20 output files
    // name the output files after the group name
    for (size_t i = 0; i < bfs1.size(); ++i) {
        std::string output_file_name = group_to_name[i + 1] + ".fasta";
        output_fasta.emplace_back(output_file_name);
        if (!output_fasta[i].is_open()) {
            std::cerr << "Unable to open output file\n";
            return 1;
        }
    }

struct result {
    size_t intersection;
    size_t num_samples;
    float ratio;
    
    // stuff needed for vector of results

};

struct output_entry {
    std::string expected_name;
    std::string predicted_name;
    std::string welch_result;
    float p_value;
    size_t intersection;
    size_t num_samples;
    float ratio;
    size_t intersection2;
    size_t num_samples2;
    float ratio2;
    std::string next_predicted_name;
    double stdErr1;
    double stdErr2;
};
output << "expected_name\tpredicted_name\twelch_result\tp_value\tintersection\tnum_samples\tratio\tintersection2\tnum_samples2\tratio2\t2nd_predicted_name\tstdErr1\tstdErr2\n";

    for (const auto& line2 : lines2) {
        std::cerr << "reading from second file" << std::endl;
         btllib::SeqReader reader(
                line2, btllib::SeqReader::Flag::LONG_MODE);
        
#pragma omp parallel
        for (const auto& record : reader) {
            auto all_kmers = get_all_methylation_kmers(record.seq, 33);
            std::vector<result> intersection(bfs1.size(), {0, 0, 0});
            //std::cerr << "calculating intersection" << std::endl;
            for (const auto& kmer : all_kmers) {
                size_t idx = kmer.first % bfSize;
                for (size_t i = 0; i < bfs1.size(); ++i) {
                    
                    if (bfs1[i][idx] == 1) {
                        intersection[i].num_samples++;
                        if (kmer.second == methylated_bfs1[i][idx]) {
                            intersection[i].intersection++;
                        }
                    }
                }
            }

            for (size_t i = 0; i < bfs1.size(); ++i) {
                if (intersection[i].num_samples == 0) {
                    continue;
                }
                intersection[i].ratio = static_cast<float>(intersection[i].intersection) / intersection[i].num_samples;
            }
            // set ratio to 0 if sample is <=10
            for (size_t i = 0; i < bfs1.size(); ++i) {
                if (intersection[i].num_samples <= 10) {
                    intersection[i].ratio = 0;
                }
            }

            //find index of largest ratio and second largest ratio
            size_t largest_group = 0;
            size_t second_largest_group = 0;
            for (size_t i = 0; i < bfs1.size(); ++i) {
                if (intersection[i].ratio > intersection[largest_group].ratio) {
                    second_largest_group = largest_group;
                    largest_group = i;
                } else if (intersection[i].ratio > intersection[second_largest_group].ratio) {
                    second_largest_group = i;
                }
            }
            // if the largest group has a ratio of 1, check and make sure that of the other groups that has a ratio of 1, it has the largest sample size
            for (size_t i = 0; i < bfs1.size(); ++i) {
                if (intersection[i].ratio == 1 && intersection[i].num_samples > intersection[largest_group].num_samples) {
                    largest_group = i;
                }
            }
            if (largest_group == second_largest_group) {
                // try to find the second largest group
                // set the second largest group to the first group that is not the largest group
                second_largest_group = 0;
                for (size_t i = 0; i < bfs1.size(); ++i) {
                    if (i != largest_group) {
                        second_largest_group = i;
                        break;
                    }
                }
                // find the group with the second largest ratio
                for (size_t i = 0; i < bfs1.size(); ++i) {
                    if (i != largest_group && intersection[i].ratio > intersection[second_largest_group].ratio) {
                        second_largest_group = i;
                    }
                }
            }
            // if the second largest group has a ratio of 1, check and make sure that of the other groups that has a ratio of 1, it has the largest sample size that is not the largest group
            for (size_t i = 0; i < bfs1.size(); ++i) {
                if (intersection[i].ratio == 1 && intersection[i].num_samples > intersection[second_largest_group].num_samples && i != largest_group) {
                    second_largest_group = i;
                }
            }

            output_entry entry;
            size_t pos3 = record.id.find("GSM");
            size_t pos4 = record.id.find("-Z000");
            std::string name2 = record.id.substr(pos3, pos4 - pos3);
            size_t pos5 = name2.find("-");
            name2 = name2.substr(pos5 + 1);
            entry.expected_name = name2;
            entry.predicted_name = group_to_name[largest_group + 1];
            entry.next_predicted_name = group_to_name[second_largest_group + 1];
            // calculate the p value for the difference in ratio using welch's t-test
            //std::cerr << "calculating welch's t-test" << std::endl;
            double mean1 = intersection[largest_group].ratio;
            double mean2 = intersection[second_largest_group].ratio;
            // if the combination (AND) of means1 and 2 are 0 and 1
            // then the standard error is 0
            if ((mean1 == 0 && mean2 == 1) || (mean1 == 1 && mean2 == 0) || (mean1 == 0 && mean2 == 0) || (mean1 == 1 && mean2 == 1) || intersection[largest_group].num_samples <= 1 || intersection[second_largest_group].num_samples <= 1 || mean1 == 0 || mean2 == 0 || mean1 == 1 || mean2 == 1) {
                entry.stdErr1 = 0;
                entry.stdErr2 = 0;
                entry.welch_result = "cannot calculate";
                entry.p_value = 1;
                entry.intersection = intersection[largest_group].intersection;
                entry.num_samples = intersection[largest_group].num_samples;
                entry.ratio = intersection[largest_group].ratio;
                entry.intersection2 = intersection[second_largest_group].intersection;
                entry.num_samples2 = intersection[second_largest_group].num_samples;
                entry.ratio2 = intersection[second_largest_group].ratio;
#pragma omp critical
                {
                    // output the results to the output file in tsv format
                    output << entry.expected_name << "\t" << entry.predicted_name << "\t" << entry.welch_result << "\t" << entry.p_value << "\t" << entry.intersection << "\t" << entry.num_samples << 
                    "\t" << entry.ratio << "\t" << entry.intersection2 << "\t" << entry.num_samples2 << "\t" << entry.ratio2 << "\t" << 
                    entry.next_predicted_name <<   "\t" << entry.stdErr1 << "\t" << entry.stdErr2 << std::endl;
                }
                continue;
            }


            double stdErr1 = std::sqrt(mean1 * (1 - mean1) / intersection[largest_group].num_samples);
            double stdErr2 = std::sqrt(mean2 * (1 - mean2) / intersection[second_largest_group].num_samples);
            double tStatistic = welchTTest(mean1, mean2, stdErr1, stdErr2);
            double df = welchDegreesOfFreedom(stdErr1, stdErr2, intersection[largest_group].num_samples, intersection[second_largest_group].num_samples);
            if (std::isnan(df)) {
                entry.stdErr1 = 0;
                entry.stdErr2 = 0;
                entry.welch_result = "cannot calculate";
                entry.p_value = 1;
                entry.intersection = intersection[largest_group].intersection;
                entry.num_samples = intersection[largest_group].num_samples;
                entry.ratio = intersection[largest_group].ratio;
                entry.intersection2 = intersection[second_largest_group].intersection;
                entry.num_samples2 = intersection[second_largest_group].num_samples;
                entry.ratio2 = intersection[second_largest_group].ratio;
#pragma omp critical
                {
                    // output the results to the output file in tsv format
                    output << entry.expected_name << "\t" << entry.predicted_name << "\t" << entry.welch_result << "\t" << entry.p_value << "\t" << entry.intersection << "\t" << entry.num_samples << 
                    "\t" << entry.ratio << "\t" << entry.intersection2 << "\t" << entry.num_samples2 << "\t" << entry.ratio2 << "\t" << 
                    entry.next_predicted_name <<   "\t" << entry.stdErr1 << "\t" << entry.stdErr2 << std::endl;
                }
                continue;
            }
            double pValue = calculatePValue(tStatistic, df);
            


            entry.welch_result = (pValue < 0.05) ? "significant" : "not significant";
            entry.p_value = pValue;
            entry.intersection = intersection[largest_group].intersection;
            entry.num_samples = intersection[largest_group].num_samples;
            entry.ratio = intersection[largest_group].ratio;
            entry.intersection2 = intersection[second_largest_group].intersection;
            entry.num_samples2 = intersection[second_largest_group].num_samples;
            entry.ratio2 = intersection[second_largest_group].ratio;
           

            
#pragma omp critical
            {
                // output the results to the output file in tsv format
                output << entry.expected_name << "\t" << entry.predicted_name << "\t" << entry.welch_result << "\t" << entry.p_value << "\t" << entry.intersection << "\t" << entry.num_samples << 
                "\t" << entry.ratio << "\t" << entry.intersection2 << "\t" << entry.num_samples2 << "\t" << entry.ratio2 << "\t" << 
                entry.next_predicted_name <<   "\t" << entry.stdErr1 << "\t" << entry.stdErr2 << std::endl;
                if (entry.welch_result == "significant") {
                    output_fasta[largest_group] << ">" << record.id << "\n" << record.seq << "\n";
                }

            }
        }

    }

    // calculate the Jaccard similarity between the two sets of k-mers in bfs1
    // against itself

/*#pragma omp parallel for
    for (size_t i = 0; i < bfs1.size(); ++i) {
        for (size_t j = i + 1; j < bfs1.size(); ++j) {
            int intersection = 0;
            int union_size = 0;
            for (size_t k = 0; k < bfSize; ++k) {
                if (bfs1[i][k] == 1 && bfs1[j][k] == 1) {
                    ++intersection;
                }
                if (bfs1[i][k] == 1 || bfs1[j][k] == 1) {
                    ++union_size;
                }
            }

            double jaccard = static_cast<double>(intersection) / union_size;
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
                output << name1 << " " << name2 << " " << jaccard << "\n";
            }
        }
    }
    
*/




    

    return 0;
}