#ifndef __SEQ_INDEXED_HPP__
#define __SEQ_INDEXED_HPP__

#include "annotation/representation/base/annotation.hpp"
#include "annotation/binary_matrix/base/binary_matrix.hpp"
#include "annotation/int_matrix/base/int_matrix.hpp"

namespace mtg::annot {

template <typename Label = std::string>
class SeqIndexedAnnotator : public MultiLabelAnnotation<Label> {
  public:
    typedef MultiLabelAnnotation<Label> Annotation;
    typedef typename Annotation::Index Index;
    typedef typename Annotation::VLabels VLabels;
    typedef matrix::MultiIntMatrix::RowTuples RowTuples;

    typedef Vector<matrix::MultiIntMatrix::Tuple> Ids;
    typedef Vector<std::pair<matrix::BinaryMatrix::Column, Ids>> SeqIds;

    SeqIndexedAnnotator(std::shared_ptr<const Annotation> base_anno,
                        std::vector<std::shared_ptr<const Annotation>>&& seq_indexes)
          : base_anno_(base_anno), seq_indexes_(std::move(seq_indexes)) {
        if (base_anno_) {
            this->label_encoder_ = base_anno_->get_label_encoder();
            assert(std::all_of(seq_indexes_.begin(), seq_indexes_.end(),
                               [](const auto &i) -> bool { return i.get(); }));
            assert(std::all_of(seq_indexes_.begin(), seq_indexes_.end(), [this](const auto &i) {
                return i->get_label_encoder().get_labels() == this->label_encoder_.get_labels();
            }));
        }
    }

    std::vector<SeqIds> get_seq_ids(const std::vector<RowTuples> &row_tuples) const;

    const matrix::BinaryMatrix& get_matrix() const override { return base_anno_->get_matrix(); }
    uint64_t num_objects() const override { return base_anno_->num_objects(); }
    uint64_t num_relations() const override { return base_anno_->num_relations(); }

    std::string file_extension() const override {
        throw std::runtime_error("No extension");
    }
    void serialize(const std::string &filename) const override {
        std::ignore = filename;
        throw std::runtime_error("Serializing SeqIndexedAnnotator not supported");
    }

    bool load(const std::string &filename) override {
        std::ignore = filename;
        throw std::runtime_error("Loading SeqIndexedAnnotator not supported");
    }
    void add_labels(const std::vector<Index> &indices,
                    const VLabels &labels) override {
        throw std::runtime_error("Not implemented for static");
    }

    void insert_rows(const std::vector<Index> &rows) override {
        throw std::runtime_error("Not implemented for static");
    }

  private:
    std::shared_ptr<const Annotation> base_anno_;
    std::vector<std::shared_ptr<const Annotation>> seq_indexes_;
};

} // namespace mtg::annot

#endif // __SEQ_INDEXED_HPP__