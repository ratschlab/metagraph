#pragma once

namespace mtg {
namespace kmer {

/**
 * What type of data structure to use in the #KmerCollector for k-mer storage.
 */
enum class ContainerType {
    VECTOR,
    /**
     * Uses several vectors that are written to disk and then merged, as defined
     * in #SortedSetDisk
     */
    VECTOR_DISK
};

} // namespace kmer
} // namespace mtg
