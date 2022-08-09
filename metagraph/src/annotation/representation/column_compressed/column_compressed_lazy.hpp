#ifndef __COLUMN_COMPRESSED_LAZY_HPP__
#define __COLUMN_COMPRESSED_LAZY_HPP__

#include <tsl/hopscotch_map.h>

#include "annotation/representation/column_compressed/annotate_column_compressed.hpp"

namespace mtg {
namespace annot {

template <typename Label = std::string>
class ColumnCompressedLazy
      : public Annotation<typename MultiLabelEncoded<Label>::Index,
                          typename MultiLabelEncoded<Label>::VLabels> {
  public:
    using Index = typename MultiLabelEncoded<Label>::Index;
    using VLabels = typename MultiLabelEncoded<Label>::VLabels;
    using LabelIndexCallback = typename MultiLabelEncoded<Label>::LabelIndexCallback;

    ColumnCompressedLazy(size_t num_objects, const std::vector<std::string> &files)
          : num_objects_(num_objects), files_(files) {}

    uint64_t num_objects() const override final { return num_objects_; }

    void call_label_objects(const VLabels &labels,
                            const LabelIndexCallback &callback,
                            size_t num_threads = 1) const override final {
        tsl::hopscotch_map<Label, size_t> label_set;
        for (size_t i = 0; i < labels.size(); ++i) {
            label_set[labels[i]] = i;
        }

        ColumnCompressed<Label>::merge_load(files_, [&](size_t, const Label &label, auto&& bitmap) {
            if (bitmap->size() != num_objects_) {
                common::logger->error("Label {} has incorrect number of rows: {} != {}",
                                      bitmap->size(), num_objects_);
                exit(1);
            }

            auto find = label_set.find(label);
            if (find != label_set.end())
                callback(find->second, *bitmap);

        }, num_threads);
    }

    const std::string& get_file(size_t i) const { return files_[i]; }

  private:
    size_t num_objects_;
    const std::vector<std::string> &files_;
};

} // namespace annot
} // namespace mtg

#endif // __COLUMN_COMPRESSED_LAZY_HPP__
