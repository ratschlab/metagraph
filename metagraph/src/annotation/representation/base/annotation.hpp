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

template <typename Label = std::string>
class LabelEncoder;


// Graph annotation
// An annotated graph is a graph with labeled nodes/edges.
// A labeled node/edge is an node/edge carrying certain labels.
template <typename LabelType>
class MultiLabelAnnotation {
  public:
    typedef uint64_t Index;
    typedef LabelType Label;
    typedef std::vector<Label> VLabels;

    virtual ~MultiLabelAnnotation() {}

    /******************* General functionality *******************/

    virtual VLabels get_labels(Index i) const final;

    virtual void add_labels(const std::vector<Index> &indices,
                            const VLabels &labels) = 0;
    // For each label and index 'indices[i]' add count 'counts[i]'
    // Note: must be thread-safe
    virtual void add_label_counts(const std::vector<Index> &indices,
                                  const VLabels &labels,
                                  const std::vector<uint64_t> &counts);
    // for each label and index 'i' add numeric attribute 'coord'
    virtual void add_label_coord(Index i, const VLabels &labels, uint64_t coord);
    // for each label and index 'i' add numeric attribute 'coord'
    virtual void add_label_coords(const std::vector<std::pair<Index, uint64_t>> &coords,
                                  const VLabels &labels);

    virtual void insert_rows(const std::vector<Index> &rows) = 0;

    // For each pair (L, L') in the dictionary, replaces label |L| with |L'|
    // and merges all relations (*, L') with matching labels L', if supported.
    virtual void rename_labels(const tsl::hopscotch_map<Label, Label> &dict);

    virtual inline const LabelEncoder<Label>& get_label_encoder() const final {
        return label_encoder_;
    }

    virtual const matrix::BinaryMatrix& get_matrix() const = 0;

    // TODO: return a shared_ptr to const bitmap
    virtual void call_objects(const Label &label,
                              std::function<void(Index)> callback) const;

    virtual void serialize(const std::string &filename) const = 0;
    virtual bool load(const std::string &filename) = 0;

    /************************* Properties *************************/

    virtual uint64_t num_objects() const = 0;
    virtual inline size_t num_labels() const {
        assert(label_encoder_.size() == get_matrix().num_columns());
        return label_encoder_.size();
    }
    virtual uint64_t num_relations() const = 0;

    virtual std::string file_extension() const = 0;

  protected:
    LabelEncoder<Label> label_encoder_;
};


// Basic dictionary mapping abstract objects to integer indexes: [0, 1, 2, ...]
template <typename Label>
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

} // namespace annot
} // namespace mtg

#endif // __ANNOTATION_HPP__
