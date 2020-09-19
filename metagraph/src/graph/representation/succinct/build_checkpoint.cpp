#include "build_checkpoint.hpp"

#include "common/logger.hpp"

namespace mtg {
namespace graph {
namespace boss {

BuildCheckpoint::BuildCheckpoint(const std::filesystem::path &output_prefix,
                                 uint32_t phase)
    : enabled_(!output_prefix.empty()),
      phase_(phase),
      checkpoint_(0),
      checkpoint_file_(output_prefix.string() + ".checkpoint") {
    if (output_prefix.empty()) {
        return;
    }
    if (std::filesystem::exists(checkpoint_file_)) {
        std::ifstream f(checkpoint_file_);
        f >> checkpoint_;
        if (checkpoint_ > 0) {
            f >> kmer_dir_;
            if (!std::filesystem::exists(kmer_dir_)) {
                common::logger->error("Found checkpoint {}, but the k-mer directory {} "
                        "doesn't exist. Remove the checkpoint to start the computation "
                        "from scratch", checkpoint_file_, kmer_dir_);
                std::exit(1);
            }
            common::logger->info(
                    "Found an interrupted computation, at phase {}, kmer directory "
                    "{}. Will attempt to continue - if this is not intended, please "
                    "remove the checkpoint file {}",
                    checkpoint_, kmer_dir_, checkpoint_file_);
        }
    }
}

void BuildCheckpoint::set_checkpoint(uint32_t checkpoint) {
    checkpoint_ = checkpoint;
    if (!enabled_)
        return;

    std::ofstream f(checkpoint_file_);
    f << checkpoint_ << std::endl;
    if (checkpoint_ > 0) {
        f << kmer_dir_;
    }
    f.close();
}

void BuildCheckpoint::remove_checkpoint() const {
    std::filesystem::remove(checkpoint_file_);
    std::filesystem::remove_all(kmer_dir_);
}

} // namespace boss
} // namespace graph
} // namespace mtg
