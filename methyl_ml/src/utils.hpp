#define XXH_INLINE_ALL

#include <algorithm>
#include <string>
#include <vector>
#include <xxhash.h>

inline XXH64_hash_t
hash(const char* seq, size_t i, unsigned k)
{
  std::string kmer(seq + i, seq + i + k);
  std::transform(kmer.begin(), kmer.end(), kmer.begin(), ::toupper);
  return XXH64(kmer.data(), kmer.size(), 42);
}