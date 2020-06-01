#include "kmc_parser.hpp"

#include <kmc_file.h>

#include "common/utils/string_utils.hpp"


namespace mtg {
namespace seq_io {

const auto kFileSuffixes = { ".kmc_suf", ".kmc_pre" };

void read_kmers(const std::string &kmc_filename,
                const std::function<void(std::string_view)> &callback,
                bool call_both_from_canonical,
                uint64_t min_count,
                uint64_t max_count) {
    read_kmers(kmc_filename,
               [&](std::string_view sequence, uint64_t) {
                   callback(sequence);
               },
               call_both_from_canonical,
               min_count,
               max_count);
}

void read_kmers(const std::string &kmc_filename,
                const std::function<void(std::string_view, uint64_t)> &callback,
                bool call_both_from_canonical,
                uint64_t min_count,
                uint64_t max_count) {
    if (min_count >= max_count)
        return;

    std::string kmc_base_filename = kmc_filename;
    for (const auto &suffix : kFileSuffixes) {
        kmc_base_filename = utils::remove_suffix(kmc_base_filename, suffix);
    }

    CKMCFile kmc_database;
    if (!kmc_database.OpenForListing(kmc_base_filename))
        throw std::runtime_error("Error: Can't open KMC database " + kmc_base_filename);

    kmc_database.SetMinCount(min_count);
    kmc_database.SetMaxCount(max_count - 1);

    size_t k = kmc_database.KmerLength();
    CKmerAPI kmer(k);
    std::string kmer_str(k, '\0');
    uint64 count;

    while (kmc_database.ReadNextKmer(kmer, count)) {
        kmer.to_string(kmer_str.data());
        callback(kmer_str, count);
        if (call_both_from_canonical && kmc_database.GetBothStrands()) {
            kmer.reverse();
            kmer.to_string(kmer_str.data());
            callback(kmer_str, count);
        }
    }
    kmc_database.Close();
}

} // namespace seq_io
} // namespace mtg
