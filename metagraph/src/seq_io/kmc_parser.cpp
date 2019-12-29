#include "kmc_parser.hpp"

#include <kmc_file.h>

#include "common/utils/string_utils.hpp"


namespace kmc {

const auto kFileSuffixes = { ".kmc_suf", ".kmc_pre" };

void read_kmers(const std::string &kmc_filename,
                const std::function<void(std::string&&)> &callback,
                bool call_both_from_canonical,
                uint64_t min_count,
                uint64_t max_count) {
    read_kmers(kmc_filename,
               [&](std::string&& sequence, uint64_t) {
                   callback(std::move(sequence));
               },
               call_both_from_canonical,
               min_count,
               max_count);
}

void read_kmers(const std::string &kmc_filename,
                const std::function<void(std::string&&, uint64_t)> &callback,
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
        throw std::runtime_error(
                std::string("Error: Can't open KMC database ") + kmc_base_filename);

    kmc_database.SetMinCount(min_count);
    kmc_database.SetMaxCount(max_count - 1);

    CKmerAPI kmer(kmc_database.KmerLength());
    uint64 count;

    while (kmc_database.ReadNextKmer(kmer, count)) {
        callback(kmer.to_string(), count);
        if (call_both_from_canonical && kmc_database.GetBothStrands()) {
            kmer.reverse();
            callback(kmer.to_string(), count);
        }
    }
    kmc_database.Close();
}

} // namespace kmc
