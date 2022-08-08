#ifndef __COLUMN_COMPRESSED_LAZY_HPP__
#define __COLUMN_COMPRESSED_LAZY_HPP__


#include "annotation/representation/column_compressed/annotate_column_compressed.hpp"

namespace mtg {
namespace annot {

template <typename Label = std::string>
class ColumnCompressedLazy : public MultiLabelEncoded<Label> {
  public:
    typedef ColumnCompressed<Label> BaseAnnotator;
    typedef typename BaseAnnotator::ColumnCallback ColumnCallback;
    using Index = typename MultiLabelEncoded<Label>::Index;
    using VLabels = typename MultiLabelEncoded<Label>::VLabels;

    ColumnCompressedLazy(size_t num_objects, const std::vector<std::string> &files)
          : num_objects_(num_objects), files_(files) {}

    uint64_t num_objects() const override final { return num_objects_; }

    void call_columns(const ColumnCallback &callback, size_t num_threads = 1) const {
        BaseAnnotator::merge_load(files_, callback, num_threads);
    }

    const std::string& get_file(size_t i) const { return files_[i]; }

    void set(Index, const VLabels &) override final { throw std::runtime_error("set not implemented"); }
    void serialize(const std::string &) const override final { throw std::runtime_error("serialize not implemented"); }
    bool load(const std::string &) override final { throw std::runtime_error("load not implemented"); }
    void add_labels(const std::vector<Index> &, const VLabels &) override final { throw std::runtime_error("add_labels not implemented"); }
    bool has_label(Index, const Label &) const override final { throw std::runtime_error("has_label not implemented"); }
    bool has_labels(Index, const VLabels &) const override final { throw std::runtime_error("has_labels not implemented"); }
    void insert_rows(const std::vector<Index> &) override final { throw std::runtime_error("insert_rows not implemented"); }
    void rename_labels(const tsl::hopscotch_map<Label, Label> &) override final { throw std::runtime_error("rename_labels not implemented"); }
    size_t num_labels() const override final { throw std::runtime_error("num_labels not implemented"); }
    uint64_t num_relations() const override final { throw std::runtime_error("num_relations not implemented"); }
    void call_objects(const Label &, std::function<void(Index)>) const override final { throw std::runtime_error("call_objects not implemented"); }
    std::string file_extension() const override final { throw std::runtime_error("file_extension not implemented"); }
    const binmat::BinaryMatrix& get_matrix() const override final { throw std::runtime_error("get_matrix not implemented"); }

  private:
    size_t num_objects_;
    const std::vector<std::string> &files_;
};

} // namespace annot
} // namespace mtg

#endif // __COLUMN_COMPRESSED_LAZY_HPP__
