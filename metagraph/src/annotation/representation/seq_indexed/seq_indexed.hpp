#ifndef __SEQ_INDEXED_HPP__
#define __SEQ_INDEXED_HPP__

#include "annotation/representation/annotation_matrix/annotation_matrix.hpp"
#include "annotation/int_matrix/base/int_matrix.hpp"

namespace mtg::annot {

template <class BinaryMatrix, typename Label = std::string>
class SeqIndexedAnnotator : public StaticBinRelAnnotator<BinaryMatrix, Label> {
  public:
    typedef StaticBinRelAnnotator<BinaryMatrix, Label> annotation_type;
    typedef matrix::MultiIntMatrix::RowTuples RowTuples;
    typedef Vector<matrix::MultiIntMatrix::Tuple> Ids;
    typedef Vector<std::pair<matrix::BinaryMatrix::Column, Ids>> SeqIds;

    SeqIndexedAnnotator(const SeqIndexedAnnotator&) = delete;
    SeqIndexedAnnotator(SeqIndexedAnnotator&&) = default;

    SeqIndexedAnnotator(std::unique_ptr<BinaryMatrix> matrix,
                        std::vector<std::unique_ptr<matrix::BinaryMatrix>>&& seq_indexes,
                        const LabelEncoder<Label> &label_encoder)
          : annotation_type(std::move(matrix), label_encoder),
            seq_indexes_(std::move(seq_indexes)) {}

    void serialize(const std::string &filename) const override {
        std::ignore = filename;
        throw std::runtime_error("Serializing SeqIndexedAnnotator not supported");
    }

    bool load(const std::string &filename) override {
        std::ignore = filename;
        throw std::runtime_error("Loading SeqIndexedAnnotator not supported");
    }

    std::vector<SeqIds> get_seq_ids(const std::vector<RowTuples> &row_tuples) const;

  private:
    std::vector<std::unique_ptr<matrix::BinaryMatrix>> seq_indexes_;
};

} // namespace mtg::annot

#endif // __SEQ_INDEXED_HPP__