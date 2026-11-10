#include "bitmap_builder.hpp"

#include <iostream>
#include <tuple>

#include "common/elias_fano/elias_fano.hpp"
#include "common/utils/file_utils.hpp"


const uint64_t kStreamBufferSize = 100'000;


bitmap_builder_set_disk::bitmap_builder_set_disk(uint64_t size,
                                                 size_t num_threads,
                                                 uint64_t buffer_size,
                                                 const std::string &swap_dir)
      : size_(size),
        num_threads_(num_threads),
        buffer_size_(buffer_size),
        swap_dir_(swap_dir),
        in_memory_(num_threads) {}

bitmap_builder_set_disk::~bitmap_builder_set_disk() {
    if (!on_disk_)
        return;

    try {
        on_disk_->clear();
        utils::remove_temp_dir(tmp_dir_);

    } catch (const std::exception &e) {
        std::cerr << "ERROR: Failed to destruct bitmap_builder_set_disk: "
                  << e.what() << std::endl;
    } catch (...) {
        std::cerr << "ERROR: Failed to destruct bitmap_builder_set_disk" << std::endl;
    }
}

void bitmap_builder_set_disk::add_ones(const uint64_t *begin, const uint64_t *end) {
    {
        std::shared_lock<std::shared_mutex> lock(mutex_);
        if (on_disk_) {
            on_disk_->insert(begin, end);
            return;
        }
    }

    std::unique_lock<std::shared_mutex> lock(mutex_);
    if (!on_disk_) {
        const uint64_t batch_size = end - begin;
        if (num_in_memory_ + batch_size > buffer_size_) {
            // drop the duplicates before deciding to go to disk
            num_in_memory_ = in_memory_.data().size();
        }
        // go to disk once the distinct positions fill half of the buffer,
        // so that the duplicates are not dropped again on every insert
        if (num_in_memory_ <= buffer_size_ / 2
                && num_in_memory_ + batch_size <= buffer_size_) {
            in_memory_.insert(begin, end);
            num_in_memory_ += batch_size;
            return;
        }
        spill();
    }
    on_disk_->insert(begin, end);
}

void bitmap_builder_set_disk::spill() {
    assert(!on_disk_);

    tmp_dir_ = utils::create_temp_dir(swap_dir_, "bitmap");
    on_disk_ = std::make_unique<mtg::common::SortedSetDisk<uint64_t>>(
            num_threads_, buffer_size_, tmp_dir_, -1, 16);

    auto &positions = in_memory_.data();
    on_disk_->insert(positions.begin(), positions.end());
    in_memory_.clear();
    num_in_memory_ = 0;
}

bitmap_builder_set_disk::InitializationData
bitmap_builder_set_disk::get_initialization_data() {
    if (merged_)
        throw std::runtime_error("ERROR: Trying to flush the bitmap twice");

    merged_ = true;

    if (!on_disk_) {
        // everything fitted in the buffer
        auto &positions = in_memory_.data();
        auto call_positions = [&positions, size = size_](auto callback) {
            std::ignore = size;
            for (uint64_t pos : positions) {
                assert(pos < size && "Indexes cannot be greater than bitmap's size");
                callback(pos);
            }
        };
        return { size_, positions.size(), call_positions };
    }

    mtg::elias_fano::EliasFanoEncoderBuffered<uint64_t> encoder(tmp_dir_/"merged",
                                                                kStreamBufferSize);
    auto &merged = on_disk_->data();
    uint64_t num_set_bits = 0;
    for (auto &it = merged.begin(); it != merged.end(); ++it) {
        assert(*it < size_ && "Indexes cannot be greater than bitmap's size");
        encoder.add(*it);
        num_set_bits++;
    }
    encoder.finish();

    on_disk_->clear();

    auto call_indexes = [&](auto callback) {
        mtg::elias_fano::EliasFanoDecoder<uint64_t> decoder(tmp_dir_/"merged", true);
        while (std::optional<uint64_t> next = decoder.next()) {
            callback(next.value());
        }
    };

    return { size_, num_set_bits, call_indexes };
}
