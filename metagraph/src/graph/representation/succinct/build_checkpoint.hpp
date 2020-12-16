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
 * In order to assist building large graphs, which can take days or even weeks to
 * construct, metagraph checkpoints the computation at certain stages in order to allow
 * continuing an interrupted computation (e.g. bc of running out of memory, or disk
 * failure, etc.).
 * Currently, metagraph checkpoints after the following operations:
 * 1. Collecting k-mers (de-duping, and sorting into chunks)
 * 2. Computing reverse complements (if built with `--canonical`)
 * 3. Splitting into 16 chunks (for computing dummy k-mers in parallel)
 * 4. Computing dummy-1 source and dummy sink k-mers
 * 5. Concatenate dummy sink and real k-mers
 * 6. Generate dummy-2+ source k-mers
 * Checkpoints 1 and 6 are exposed to the client as phase 1 and 2, so that users can
 * break down graph construction into stages. The following three phases are currently
 * supported:
 * Phase 1: kmer collection. This phase stops after having collected, sorted and de-duped
 * k-mers from all the input files. If you have a large amount of input files, you can
 * shard them into n chunks, and process each chunk on a different machine (doesn't need
 * much RAM or CPU), making sure you specify the same (empty) temporary directory
 * using --disk-swap.
 * Phase 2: kmer generation. This phase stops after having generated reverse complements
 * (if --canonical is present) and dummy source and sink k-mers. This needs to run on a
 * single machine, but doesn't need much RAM.
 * Phase 3: graph generation. This is the final phase and consists of processing the
 * k-mers and generating the BOSS table. This phase requires enough RAM to fit the BOSS
 * table, which is 0.6 * kmer_count bytes. If using `--count-kmers`, additional memory
 * to store the count for each k-mer is needed. For `--count-width 16` (2 bytes), the
 * memory used by the BOSS table is 2.6 * kmer_count.
 */
class BuildCheckpoint {
  public:
    BuildCheckpoint(const std::filesystem::path &output_prefix, uint32_t phase);

    uint32_t phase() const { return phase_; }

    uint32_t checkpoint() const { return checkpoint_; }
    void set_checkpoint(uint32_t checkpoint);

    static uint32_t checkpoint_for_phase(uint32_t phase) {
        return phase <= 1 ? phase : 6;
    }

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
