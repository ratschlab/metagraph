#pragma once

namespace mg {
namespace kmer {

/**
 * What type of data structure to use in the #KmerCollector for k-mer storage.
 * Note: for k-mer counting, only VECTOR is supported.
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
} // namespace mg
