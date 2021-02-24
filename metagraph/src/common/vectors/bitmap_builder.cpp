#include "bitmap_builder.hpp"

#include "common/elias_fano.hpp"
#include "common/utils/file_utils.hpp"


const uint64_t kStreamBufferSize = 100'000;


bitmap_builder_set_disk::bitmap_builder_set_disk(uint64_t size,
                                                 size_t num_threads,
                                                 uint64_t buffer_size_in_bytes,
                                                 const std::string &swap_dir)
      : size_(size),
        swap_dir_(swap_dir),
        chunks_tmp_dir_(utils::create_temp_dir(swap_dir, "bitmap_chunks")),
        set_bit_positions_(num_threads, buffer_size_in_bytes / sizeof(uint64_t),
                           chunks_tmp_dir_, -1, 16) {}

bitmap_builder_set_disk::~bitmap_builder_set_disk() {
    utils::remove_temp_dir(chunks_tmp_dir_);
}

bitmap_builder_set_disk::InitializationData
bitmap_builder_set_disk::get_initialization_data() {
    if (merged_)
        throw std::runtime_error("ERROR: Trying to flush the bitmap twice");

    // create temp dir for the merged result
    auto tmp_dir = utils::create_temp_dir(swap_dir_, "bitmap");

    mtg::common::EliasFanoEncoderBuffered<uint64_t> encoder(tmp_dir/"merged",
                                                            kStreamBufferSize);
    auto &merged = set_bit_positions_.data();
    uint64_t num_set_bits = 0;
    for (auto &it = merged.begin(); it != merged.end(); ++it) {
        assert(*it < size_ && "Indexes cannot be greater than bitmap's size");
        encoder.add(*it);
        num_set_bits++;
    }
    encoder.finish();

    set_bit_positions_.clear();

    auto call_indexes = [&](auto callback) {
        mtg::common::EliasFanoDecoder<uint64_t> decoder(tmp_dir/"merged", true);
        while (std::optional<uint64_t> next = decoder.next()) {
            callback(next.value());
        }
        utils::remove_temp_dir(tmp_dir);
    };

    merged_ = true;

    return { size_, num_set_bits, call_indexes };
}
