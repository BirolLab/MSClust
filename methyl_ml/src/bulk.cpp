#include <argparse/argparse.hpp>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <kseq.h>
#include <npy.hpp>
#include <omp.h>
#include <regex>
#include <string>
#include <unistd.h>
#include <vector>

#include "utils.hpp"

KSEQ_INIT(int, read)

class ProgramArguments
{
public:
  unsigned kmer_length;
  unsigned num_threads;
  std::string reads_path;
  std::string vectors_path;
  std::string out_path;

  ProgramArguments(int argc, char* argv[])
  {
    argparse::ArgumentParser parser("onehot");
    parser.add_argument("-r").help("path to bulk reads file").required();
    parser.add_argument("-v").help("path to vectors file").required();
    parser.add_argument("-k")
      .help("k-mer length")
      .scan<'u', unsigned>()
      .default_value(26U);
    parser.add_argument("-t")
      .help("number of threads")
      .scan<'u', unsigned>()
      .default_value(1U);
    parser.add_argument("-o").help("output path").required();
    try {
      parser.parse_args(argc, argv);
    } catch (const std::exception& err) {
      std::cerr << err.what() << std::endl;
      std::cerr << parser;
      std::exit(1);
    }
    reads_path = parser.get("-r");
    vectors_path = parser.get("-v");
    kmer_length = parser.get<unsigned>("-k");
    num_threads = parser.get<unsigned>("-t");
    out_path = parser.get("-o");
  }
};

std::vector<std::pair<size_t, int>>
get_sparse_vector(char* seq,
                  size_t seq_len,
                  size_t kmer_length,
                  size_t num_features)
{
  std::vector<std::pair<size_t, int>> vec;
  for (size_t i = 0; i < seq_len - kmer_length; i++) {
    if ((seq[i] == '1' || seq[i] == 'C' || seq[i] == 'c') &&
        (seq[i + 1] == 'G' || seq[i + 1] == 'g')) {
      size_t k_i = std::max(0L, (ssize_t)i - (ssize_t)kmer_length / 2);
      const char backup = seq[i];
      seq[i] = 'C';
      size_t feature = hash(seq, k_i, kmer_length) % num_features;
      seq[i] = backup;
      vec.emplace_back(std::make_pair(feature, seq[i] == '1' ? 1 : -1));
    }
  }
  return vec;
}

double
norm(const int* vec, size_t num_features)
{
  double n = 0;
  for (size_t i = 0; i < num_features; i++) {
    n += vec[i] * vec[i];
  }
  return n;
}

std::string
get_read_label(const std::string& read_name)
{
  std::regex rgx("GSM\\d+-(.+?)-Z0+");
  std::smatch matches;
  std::regex_search(read_name, matches, rgx);
  return matches[1].str();
}

int
main(int argc, char* argv[])
{
  const auto args = ProgramArguments(argc, argv);
  omp_set_num_threads(args.num_threads);
  std::cout << "loading vectors... " << std::flush;
  const auto vectors = npy::read_npy<int>(args.vectors_path);
  std::cout << "done" << std::endl;
  size_t num_samples = vectors.shape[0], num_features = vectors.shape[1];
  std::cout << num_samples << "x" << num_features << std::endl;
  std::cout << "calculating vector norms... " << std::flush;
  double* vec_norms = new double[num_samples];
#pragma omp parallel for
  for (size_t i = 0; i < num_samples; i++) {
    vec_norms[i] = norm(vectors.data.data() + (i * num_features), num_features);
  }
  std::cout << "done" << std::endl;
  std::cout << "calculating distances..." << std::endl;
  unsigned num_reads = 0;
  std::vector<double> distances;
  FILE* fp = fopen(args.reads_path.data(), "r");
  kseq_t* seq = kseq_init(fileno(fp));
  std::ofstream labels(args.out_path + ".txt");
  while (kseq_read(seq) >= 0) {
    const auto vec =
      get_sparse_vector(seq->seq.s, seq->seq.l, args.kmer_length, num_features);
    const auto norm_a = vec.size();
    for (size_t i = 0; i < num_samples; i++) {
      double a_dot_b = 0;
      for (const auto& item : vec) {
        a_dot_b += item.second * vectors.data[i * num_features + item.first];
      }
      distances.push_back(1 - a_dot_b / (norm_a * vec_norms[i]));
    }
    std::string read_name(seq->name.s, seq->name.l);
    labels << get_read_label(read_name) << std::endl;
    ++num_reads;
    if (num_reads % 1000 == 0) {
      std::cout << "  processed " << num_reads << " reads" << std::endl;
      break;
    }
  }
  kseq_destroy(seq);
  fclose(fp);
  std::cout << "done processing " << num_reads << " reads" << std::endl;
  std::cout << "saving distances... " << std::flush;
  npy::npy_data_ptr<double> dists_ptr;
  dists_ptr.data_ptr = distances.data();
  dists_ptr.shape = { num_reads, num_samples };
  dists_ptr.fortran_order = false;
  npy::write_npy(args.out_path + ".npy", dists_ptr);
  std::cout << "done" << std::endl;
  return EXIT_SUCCESS;
}