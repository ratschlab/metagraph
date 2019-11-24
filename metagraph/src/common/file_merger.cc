#include "file_merger.hh"

namespace mg {
namespace common {
static void merge_file(const std::vector<std::string> &sources,
                       const std::string &out_file,
                       ChunkedWaitQueue<T> *result) {
    // start merging disk chunks by using a heap to store the current element
    // from each chunk
    std::vector<std::fstream> chunk_files(chunk_count);

    auto comp = [](const auto &a, const auto &b) { return a.first > b.first; };

    std::priority_queue<std::pair<T, uint32_t>, std::vector<std::pair<T, uint32_t>>, decltype(comp)>
            merge_heap(comp);
    for (uint32_t i = 0; i < chunk_count; ++i) {
        chunk_files[i].open(sources[i], std::ios::in | std::ios::binary);
        T data_item;
        if (chunk_files[i].good()) {
            chunk_files[i].read(reinterpret_cast<char *>(&data_item), sizeof(data_item));
            merge_heap.push({ data_item, i });
        } else {
            throw std::runtime_error("Unable to open chunk file " + sources[i]);
        }
    }
    uint64_t totalSize = 0;

    // initialized to suppress maybe-uninitialized warnings in GCC
    T last_written = {};
    std::fstream sorted_file;
    if (out_file != "") {
        sorted_file = std::fstream(out_file, std::ios::binary | std::ios::out);
        if (!sorted_file) {
            std::cerr << "Error: Could not create file " << out_file << std::endl;
            std::exit(EXIT_FAILURE);
        }
    }
    bool has_written = false;
    while (!merge_heap.empty()) {
        std::pair<T, uint32_t> smallest = merge_heap.top();
        merge_heap.pop();
        if (!has_written || smallest.first != last_written) {
            has_written = true;
            merge_queue->push_front(smallest.first);
            if (out_file != "") {
                // TODO(ddanciu) - consider writing asynchronously if this proves to
                // be a bottleneck (a simple produce/consumer queue should work)
                if (!sorted_file.write(reinterpret_cast<char *>(&smallest.first), sizeof(T))) {
                    std::cerr << "Error: Writing of merged data to " + out_file
                                    + " failed."
                              << std::endl;
                    std::exit(EXIT_FAILURE);
                }
            }

            last_written = smallest.first;
            totalSize++;
        }
        if (chunk_files[smallest.second].good()) {
            T data_item;
            if (chunk_files[smallest.second].read(reinterpret_cast<char *>(&data_item),
                                                  sizeof(data_item))) {
                merge_heap.push({ data_item, smallest.second });
            }
        }
    }
    merge_queue->shutdown();
    sorted_file.close();

    for (uint32_t i = 0; i < chunk_count; ++i) {
        std::filesystem::remove(sources[i]);
    }
}
} // namespace common
} // namespace mg
