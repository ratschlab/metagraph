#ifndef __ANNOTATION_HPP__
#define __ANNOTATION_HPP__

#include <cassert>
#include <cstdint>
#include <string>
#include <vector>
#include <memory>
#include <functional>

#include <tsl/hopscotch_map.h>

#include "common/vector.hpp"
#include "annotation/binary_matrix/base/binary_matrix.hpp"


namespace mtg {
namespace annot {


template <typename IndexType, typename LabelType>
class AnnotationCategory {
  public:
    virtual ~AnnotationCategory() {}

    virtual LabelType get(IndexType i) const = 0;
    virtual void set(IndexType i, const LabelType &label) = 0;

    virtual void serialize(const std::string &filename) const = 0;
    virtual bool load(const std::string &filename) { return merge_load({ filename }); }
    // TODO: remove merge_load and merge by annotation converters/mergers?
    virtual bool merge_load(const std::vector<std::string> &filenames) = 0;
};


// Graph annotation
// An annotated graph is a graph with labeled nodes/edges.
// A labeled node/edge is an node/edge carrying certain labels.
template <typename IndexType, typename LabelType>
class MultiLabelAnnotation
      : public AnnotationCategory<IndexType, std::vector<LabelType>> {

  public:
    typedef IndexType Index;
    typedef LabelType Label;
    typedef std::vector<Label> VLabels;

    virtual ~MultiLabelAnnotation() {}

    /******************* General functionality *******************/

    virtual void add_labels(const std::vector<Index> &indices,
                            const VLabels &labels) = 0;

    virtual bool has_label(Index i, const Label &label) const = 0;
    virtual bool has_labels(Index i, const VLabels &labels) const = 0;

    virtual void insert_rows(const std::vector<Index> &rows) = 0;

    // For each pair (L, L') in the dictionary, replaces label |L| with |L'|
    // and merges all relations (*, L') with matching labels L', if supported.
    virtual void rename_labels(const tsl::hopscotch_map<Label, Label> &dict) = 0;

    /************************* Properties *************************/

    virtual uint64_t num_objects() const = 0;
    virtual size_t num_labels() const = 0;
    virtual uint64_t num_relations() const = 0;
    virtual const std::vector<Label>& get_all_labels() const = 0;

    virtual bool label_exists(const Label &label) const = 0;

    virtual void call_objects(const Label &label,
                              std::function<void(Index)> callback) const = 0;

    virtual std::string file_extension() const = 0;
};


// Basic dictionary mapping abstract objects to integer indexes: [0, 1, 2, ...]
template <typename Label = std::string>
class LabelEncoder {
  public:
    /**
     * If the label passed does not exist, insert
     * that label and return its code.
     */
    size_t insert_and_encode(const Label &label);

    /**
     * Return the code of the label passed.
     * Throw exception if it does not exist.
     */
    size_t encode(const Label &label) const { return encode_label_.at(label); }

    /**
     * Check if the label has been added to encoder.
     */
    bool label_exists(const Label &label) const { return encode_label_.count(label); }

    /**
     * Throws an exception if a bad code is passed.
     */
    const Label& decode(size_t code) const { return decode_label_.at(code); }

    const std::vector<Label>& get_labels() const { return decode_label_; }

    void merge(const LabelEncoder<Label> &other);

    size_t size() const { return decode_label_.size(); }

    void clear() { encode_label_.clear(); decode_label_.clear(); }

    bool load(std::istream &instream);
    void serialize(std::ostream &outstream) const;

  private:
    tsl::hopscotch_map<Label, uint64_t> encode_label_;
    std::vector<Label> decode_label_;
};


template <typename LabelType>
class MultiLabelEncoded : public MultiLabelAnnotation<uint64_t, LabelType> {
  public:
    using Index = typename MultiLabelAnnotation<uint64_t, LabelType>::Index;
    using Label = typename MultiLabelAnnotation<uint64_t, LabelType>::Label;
    using VLabels = typename MultiLabelAnnotation<uint64_t, LabelType>::VLabels;

    virtual ~MultiLabelEncoded() {}

    /******************* General functionality *******************/

    virtual VLabels get(Index i) const override final;

    virtual inline const LabelEncoder<Label>& get_label_encoder() const final {
        return label_encoder_;
    }

    // For each pair (L, L') in the dictionary, replaces label |L| with |L'|
    // and merges all relations (*, L') with matching labels L', if supported.
    virtual void rename_labels(const tsl::hopscotch_map<Label, Label> &dict) override;

    /************************* Properties *************************/

    virtual inline size_t num_labels() const override {
        assert(label_encoder_.size() == get_matrix().num_columns());
        return label_encoder_.size();
    }

    virtual const std::vector<Label>& get_all_labels() const override final {
        return label_encoder_.get_labels();
    }

    virtual inline bool label_exists(const Label &label) const override final {
        return label_encoder_.label_exists(label);
    }

    // TODO: return a shared_ptr to const bitmap
    virtual void call_objects(const Label &label,
                              std::function<void(Index)> callback) const override;

    virtual const binmat::BinaryMatrix& get_matrix() const = 0;

    /*********************** Special queries **********************/

    /**
     * Return all labels for which counts are greater than or equal to |min_count|.
     * Stop counting if count is greater than |count_cap|.
     */
    virtual std::vector<std::pair<uint64_t /* label_code */, size_t /* count */>>
    count_labels(const std::vector<std::pair<Index, size_t>> &index_counts,
                 size_t min_count = 1,
                 size_t count_cap = std::numeric_limits<size_t>::max()) const;

  protected:
    LabelEncoder<Label> label_encoder_;
};

} // namespace annot
} // namespace mtg

#endif // __ANNOTATION_HPP__
