#include "bitmap_builder.hpp"

#include "common/elias_fano.hpp"


const uint64_t kStreamBufferSize = 100'000;


bitmap_builder_set_disk::InitializationData
bitmap_builder_set_disk::get_initialization_data() const {
    if (merged_)
        throw std::runtime_error("ERROR: Trying to flush the bitmap twice");

    mtg::common::EliasFanoEncoderBuffered<uint64_t> encoder(tmp_dir_/"merged",
                                                            kStreamBufferSize);
    auto &merged = const_cast<bitmap_builder_set_disk*>(this)->set_bit_positions_.data();
    uint64_t num_set_bits = 0;
    for (auto &it = merged.begin(); it != merged.end(); ++it) {
        assert(*it < size_ && "Indexes cannot be greater than bitmap's size");
        encoder.add(*it);
        num_set_bits++;
    }
    encoder.finish();

    const_cast<bitmap_builder_set_disk*>(this)->set_bit_positions_.clear();

    auto call_indexes = [&](auto callback) {
        mtg::common::EliasFanoDecoder<uint64_t> decoder(tmp_dir_/"merged", true);
        while (std::optional<uint64_t> next = decoder.next()) {
            callback(next.value());
        }
    };

    const_cast<bitmap_builder_set_disk*>(this)->merged_ = true;

    return { size_, num_set_bits, call_indexes };
}
