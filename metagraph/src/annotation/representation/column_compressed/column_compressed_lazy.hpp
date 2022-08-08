#ifndef __COLUMN_COMPRESSED_LAZY_HPP__
#define __COLUMN_COMPRESSED_LAZY_HPP__


#include "annotation/representation/column_compressed/annotate_column_compressed.hpp"

namespace mtg {
namespace annot {

template <typename Label = std::string>
class ColumnCompressedLazy
      : public AnnotationCategory<typename MultiLabelEncoded<Label>::Index,
                                  typename MultiLabelEncoded<Label>::VLabels> {
  public:
    using Index = typename MultiLabelEncoded<Label>::Index;
    using VLabels = typename MultiLabelEncoded<Label>::VLabels;
    using LabelIndexCallback = typename MultiLabelEncoded<Label>::LabelIndexCallback;

    ColumnCompressedLazy(size_t num_objects, const std::vector<std::string> &files)
          : num_objects_(num_objects), files_(files) {}

    void call_label_indices(const LabelIndexCallback &callback, size_t num_threads = 1) const override final {
        ColumnCompressed<Label>::merge_load(files_, [&](size_t, const Label &label, auto&& bitmap) {
            if (bitmap->size() != num_objects_) {
                common::logger->error("Label {} has incorrect number of rows: {} != {}",
                                      bitmap->size(), num_objects_);
                exit(1);
            }
            callback(label, *bitmap);
        }, num_threads);
    }

    const std::string& get_file(size_t i) const { return files_[i]; }

    VLabels get(Index) const override final { throw std::runtime_error("get not implemented"); }
    void set(Index, const VLabels &) override final { throw std::runtime_error("set not implemented"); }

    void serialize(const std::string &) const override final {
        common::logger->warn("Nothing to serialize");
    }

    bool load(const std::string &) override final {
        common::logger->warn("Nothing to load");
        return true;
    }

  private:
    size_t num_objects_;
    const std::vector<std::string> &files_;
};

} // namespace annot
} // namespace mtg

#endif // __COLUMN_COMPRESSED_LAZY_HPP__
