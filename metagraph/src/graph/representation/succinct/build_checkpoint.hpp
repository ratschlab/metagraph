#pragma once
#include <cstdint>
#include <filesystem>
#include <fstream>

namespace mtg {
namespace graph {
namespace boss {

/**
 * Stores checkpointing and phase information for resuming disk-based building of succinct
 * graphs.
 */
class BuildCheckpoint {
  public:
    BuildCheckpoint(const std::filesystem::path &output_prefix, uint32_t phase);

    uint32_t phase() const { return phase_; }

    uint32_t checkpoint() const { return checkpoint_; }
    void set_checkpoint(uint32_t checkpoint);

    const std::filesystem::path& tmp_dir() const { return checkpoint_file_; }

    const std::filesystem::path& kmer_dir() const { return kmer_dir_; }
    void set_kmer_dir(const std::filesystem::path &kmer_dir) { kmer_dir_ = kmer_dir; }

    void remove_checkpoint() const;

  private:
    void store() const;

    bool enabled_;
    uint32_t phase_;
    uint32_t checkpoint_;
    std::filesystem::path checkpoint_file_;
    std::filesystem::path kmer_dir_;
};

} // namespace boss
} // namespace graph
} // namespace mtg
