#pragma once
#include <cstdint>
#include <filesystem>
#include <fstream>

namespace mtg {
namespace graph {
namespace boss {

/**
 * Stores checkpointing information for resuming disk-based building of succinct graphs.
 */
class BuildCheckpoint {
  public:
    BuildCheckpoint(bool enabled, const std::filesystem::path &output_prefix, uint32_t phase);

    uint32_t phase() const { return phase_; }
    uint32_t checkpoint() const { return checkpoint_; }
    std::filesystem::path tmp_dir() const { return checkpoint_file_; }
    std::filesystem::path kmer_dir() const { return kmer_dir_; }

    void set_kmer_dir(const std::filesystem::path &kmer_dir) { kmer_dir_ = kmer_dir; }

    void set_checkpoint(uint32_t checkpoint) { checkpoint_ = checkpoint; }

    void done() { std::filesystem::remove(checkpoint_file_); }

    void store() const;

  private:
    bool enabled_;
    uint32_t phase_;
    uint32_t checkpoint_;
    std::filesystem::path checkpoint_file_;
    std::filesystem::path kmer_dir_;
};

} // namespace boss
} // namespace graph
} // namespace mtg
