#include "kmc_parser.hpp"

#include <kmc_file.h>

#include "utils.hpp"


namespace kmc {

const auto kFileSuffixes = { ".kmc_suf", ".kmc_pre" };


void read_kmers(const std::string &kmc_filename,
                const std::function<void(std::string&&)> &callback,
                uint64_t min_count) {
    std::string kmc_base_filename = kmc_filename;
    for (const auto &suffix : kFileSuffixes) {
        kmc_base_filename = utils::remove_suffix(kmc_base_filename, suffix);
    }

    CKMCFile kmc_database;
    if (!kmc_database.OpenForListing(kmc_base_filename))
        throw std::runtime_error(
                std::string("Error: Can't open KMC database ") + kmc_base_filename);

    kmc_database.SetMinCount(min_count);

    CKmerAPI kmer(kmc_database.KmerLength());
    uint64 count;

    while (kmc_database.ReadNextKmer(kmer, count)) {
        callback(kmer.to_string());
        if (kmc_database.GetBothStrands()) {
            kmer.reverse();
            callback(kmer.to_string());
        }
    }
    kmc_database.Close();
}

} // namespace kmc
