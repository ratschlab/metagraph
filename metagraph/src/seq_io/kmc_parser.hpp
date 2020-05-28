#ifndef __KMC_KMERS__
#define __KMC_KMERS__

#include <functional>
#include <string>


namespace mtg {
namespace seq_io {

// Read k-mers from KMC database.
// This retrieves only canonical k-mers from canonical
// KMC databases (those constructed by KMC without flag '-b' passed).
//
// - If |call_both_from_canonical| = true and the KMC database is canonical
//   (stores only canonical k-mers with counts), call k-mers with their
//   reverse-complement (non-canonical) ones.
//   This flag does not change the behavior for reading from non-canonical
//   KMC databases (those consructed by KMC with flag '-b' passed).
//
// - |min_count| -- minimum k-mer abundance (including the value passed)
// - |max_count| -- maximum k-mer abundance (excluding the value passed)
void read_kmers(const std::string &kmc_base_filename,
                const std::function<void(std::string_view, uint64_t count)> &callback,
                bool call_both_from_canonical,
                uint64_t min_count = 1,
                uint64_t max_count = -1);

void read_kmers(const std::string &kmc_filename,
                const std::function<void(std::string_view)> &callback,
                bool call_both_from_canonical,
                uint64_t min_count = 1,
                uint64_t max_count = -1);

} // namespace seq_io
} // namespace mtg

#endif // __KMC_KMERS__
