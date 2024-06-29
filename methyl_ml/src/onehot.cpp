#include <argparse/argparse.hpp>
#include <filesystem>
#include <fstream>
#include <glob/glob.hpp>
#include <indicators/progress_bar.hpp>
#include <indicators/terminal_size.hpp>
#include <iostream>
#include <kseq.h>
#include <npy.hpp>
#include <omp.h>
#include <string>
#include <vector>

#include "utils.hpp"

KSEQ_INIT(int, read)

class ProgramArguments {
public:
  std::vector<std::string> read_paths;
  unsigned num_features;
  unsigned kmer_length;
  unsigned num_threads;
  std::string out_path;

  ProgramArguments(int argc, char *argv[]) {
    argparse::ArgumentParser parser("onehot");
    parser.add_argument("-d")
        .help("feature vector size")
        .required()
        .scan<'u', unsigned>();
    parser.add_argument("-k")
        .help("k-mer length")
        .scan<'u', unsigned>()
        .default_value(26U);
    parser.add_argument("-t")
        .help("number of threads")
        .scan<'u', unsigned>()
        .default_value(1U);
    parser.add_argument("-o").help("output path").required();
    parser.add_argument("samples").help("paths to sample folders").remaining();
    try {
      parser.parse_args(argc, argv);
    } catch (const std::exception &err) {
      std::cerr << err.what() << std::endl;
      std::cerr << parser;
      std::exit(1);
    }
    num_features = parser.get<unsigned>("-d");
    kmer_length = parser.get<unsigned>("-k");
    num_threads = parser.get<unsigned>("-t");
    out_path = parser.get("-o");
    const auto samples = parser.get<std::vector<std::string>>("samples");
    for (const auto &sample : glob::glob(samples)) {
      const auto cells = std::filesystem::recursive_directory_iterator(sample);
      for (const auto &cell : cells) {
        const auto reads = cell / std::filesystem::path("_aligned_reads.fasta");
        if (cell.is_directory() && std::filesystem::exists(reads)) {
          read_paths.emplace_back(reads.u8string());
        }
      }
    }
  }
};

class Dataset {
private:
  std::vector<int> data;
  std::vector<std::string> labels;
  size_t num_cells, num_features;

public:
  Dataset(size_t num_cells, size_t num_features)
      : data(std::vector<int>(num_cells * num_features, 0)),
        labels(std::vector<std::string>(num_cells, "")), num_cells(num_cells),
        num_features(num_features) {}

  void set(size_t cell, std::string label) { labels[cell] = label; }

  void set(size_t cell, size_t feature, int value) {
    data[cell * num_features + feature] = value;
  }

  void save(const std::string &path) {
    std::filesystem::path out_dir(path);
    std::filesystem::path vec_file("vectors.npy");
    std::filesystem::path lbl_file("labels.txt");
    npy::npy_data_ptr<int> np_arr;
    np_arr.data_ptr = data.data();
    np_arr.shape = {num_cells, num_features};
    np_arr.fortran_order = false;
    npy::write_npy(out_dir / vec_file, np_arr);
    std::ofstream labels_file(out_dir / lbl_file);
    for (const auto &label : labels) {
      labels_file << label << std::endl;
    }
  }

  size_t get_num_features() { return num_features; }
};

inline void process_file(const std::string &path, Dataset &data,
                         unsigned cell_number, unsigned kmer_length) {
  const auto num_features = data.get_num_features();
  FILE *fp = fopen(path.data(), "r");
  kseq_t *seq = kseq_init(fileno(fp));
  while (kseq_read(seq) >= 0) {
    for (size_t i = 0; i < seq->seq.l - kmer_length; i++) {
      const auto seq_s = seq->seq.s;
      if ((seq_s[i] == '1' || seq_s[i] == 'C' || seq_s[i] == 'c') &&
          (seq_s[i + 1] == 'G' || seq_s[i + 1] == 'g')) {
        size_t k_i = std::max(0L, (ssize_t)i - (ssize_t)kmer_length / 2);
        const char backup = seq_s[i];
        seq_s[i] = 'C';
        size_t feature = hash(seq_s, k_i, kmer_length) % num_features;
        seq_s[i] = backup;
        data.set(cell_number, feature, seq_s[i] == '1' ? 1 : -1);
      }
    }
  }
  kseq_destroy(seq);
  fclose(fp);
}

int main(int argc, char *argv[]) {
  const auto args = ProgramArguments(argc, argv);
  omp_set_num_threads(args.num_threads);
  std::cout << "allocating... " << std::flush;
  Dataset data(args.read_paths.size(), args.num_features);
  std::cout << "done" << std::endl;
  indicators::ProgressBar bar{
      indicators::option::BarWidth{indicators::terminal_width() - 30},
      indicators::option::Start{"["},
      indicators::option::Fill{"="},
      indicators::option::Lead{">"},
      indicators::option::Remainder{" "},
      indicators::option::End{"]"},
      indicators::option::ShowElapsedTime{true},
      indicators::option::ShowRemainingTime{true},
  };
  unsigned files_done = 0;
#pragma omp parallel for
  for (size_t cell = 0; cell < args.read_paths.size(); cell++) {
    const auto file = std::filesystem::path(args.read_paths[cell]);
    const auto label = file.parent_path().parent_path().filename();
    data.set(cell, label);
    process_file(file, data, cell, args.kmer_length);
#pragma omp critical
    {
      ++files_done;
      const auto r = (double)files_done / (double)args.read_paths.size();
      const auto prefix = std::to_string(files_done) + "/" +
                          std::to_string(args.read_paths.size()) + " ";
      bar.set_option(indicators::option::PrefixText{prefix});
      bar.set_progress(r * 100.0);
    }
  }
  std::cout << "saving... " << std::flush;
  data.save(args.out_path);
  std::cout << "done" << std::endl;
  return EXIT_SUCCESS;
}