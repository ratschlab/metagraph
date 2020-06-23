#include "formats.hpp"

#include <algorithm>


namespace mtg {
namespace seq_io {

std::string file_format(const std::string &fname) {
    size_t dotind = fname.rfind(".");
    if (dotind == std::string::npos)
        return "";

    std::string ext = fname.substr(dotind);

    if (ext == ".kmc_pre" || ext == ".kmc_suf")
        return "KMC";

    if (ext == ".gz" || ext == ".bgz") {
        size_t nextind = fname.substr(0, dotind - 1).rfind(".");
        if (nextind == std::string::npos)
            return "";

        ext = fname.substr(nextind, dotind - nextind);
    }

    std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);

    if (ext == ".vcf") {
        return "VCF";

    } else if ((ext == ".fq") || (ext == ".fastq")) {
        return "FASTQ";

    } else {
        return "FASTA";
    }
}

} // namespace seq_io
} // namespace mtg
