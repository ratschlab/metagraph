#ifndef __ANNOTATE_HPP__
#define __ANNOTATE_HPP__

#include <cstdint>
#include <string>
#include <vector>
#include <unordered_map>


namespace annotate {

// class GenomeAnnotation {
//   public:
//     typedef uint64_t Index;
//     enum Label {
//         OUTGOING_EDGE_LABELS = 0,
//         SEQUENCE_NAME,
//     };

//     virtual const Annotation& get(Index i, Label j) const = 0;
//     virtual void set(Index i, Label j, const Annotation &label) = 0;
//     virtual void merge(const GenomeAnnotation &annotation) = 0;

//     virtual bool load(const std::string &filename) = 0;
//     virtual void serialize(const std::string &filename) const = 0;
// };


template <typename Index, typename LabelType>
class AnnotationCategory {
  public:
    virtual ~AnnotationCategory() {}

    virtual LabelType get(Index i) const = 0;
    virtual void add(Index i, const LabelType &label) = 0;
    virtual void set(Index i, const LabelType &label) = 0;

    virtual void serialize(const std::string &filename) const = 0;
    virtual bool load(const std::string &filename) { return merge_load({ filename }); }
    virtual bool merge_load(const std::vector<std::string> &filenames) = 0;
};


// Graph annotation
// An annotated graph is a graph with labeled edges.
// A labeled edge is an edge carrying certain labels.
template <typename IndexType, typename LabelType>
class MultiLabelAnnotation
      : public AnnotationCategory<IndexType, std::vector<LabelType>> {

  public:
    typedef IndexType Index;
    typedef LabelType Label;
    typedef std::vector<Label> VLabels;

    virtual ~MultiLabelAnnotation() {}

    /***************** Inherited member functions ****************/

    virtual VLabels get(Index i) const override final {
        return get_labels(i);
    }
    virtual void set(Index i, const VLabels &labels) override final {
        set_labels(i, labels);
    }
    virtual void add(Index i, const VLabels &labels) override final {
        add_labels(i, labels);
    }

    /******************* General functionality *******************/

    virtual void set_labels(Index i, const VLabels &labels) = 0;
    virtual VLabels get_labels(Index i) const = 0;

    virtual void add_label(Index i, const Label &label) = 0;
    virtual void add_labels(Index i, const VLabels &labels) = 0;
    virtual void add_labels(const std::vector<Index> &indices,
                            const VLabels &labels) = 0;

    virtual bool has_label(Index i, const Label &label) const = 0;
    virtual bool has_labels(Index i, const VLabels &labels) const = 0;

    virtual void serialize(const std::string &filename) const override = 0;
    virtual bool merge_load(const std::vector<std::string> &filenames) override = 0;

    virtual void insert_rows(const std::vector<Index> &rows) = 0;

    /*********************** Special queries **********************/

    // Get labels that occur at least in |presence_ratio| rows.
    // If |presence_ratio| = 0, return all occurring labels.
    virtual VLabels get_labels(const std::vector<Index> &indices,
                               double presence_ratio) const = 0;

    // Count all labels collected from the given rows
    // and return top |num_top| with the their counts.
    virtual std::vector<std::pair<Label, size_t>>
    get_top_labels(const std::vector<Index> &indices,
                   size_t num_top = static_cast<size_t>(-1)) const = 0;

    /************************* Properties *************************/

    virtual uint64_t num_objects() const = 0;
    virtual size_t num_labels() const = 0;
    virtual uint64_t num_relations() const = 0;
};


// A dictionary to encode annotation labels
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
    size_t encode(const Label &label) const;

    /**
     * Throws an exception if a bad code is passed.
     */
    const Label& decode(size_t code) const { return decode_label_.at(code); }

    size_t size() const { return decode_label_.size(); }

    bool load(std::istream &instream);
    void serialize(std::ostream &outstream) const;

    void clear() { encode_label_.clear(); decode_label_.clear(); }

  private:
    std::unordered_map<Label, uint64_t> encode_label_;
    std::vector<Label> decode_label_;
};


template <typename IndexType, typename LabelType>
class MultiLabelEncoded
      : public MultiLabelAnnotation<IndexType, LabelType> {

  public:
    using Index = typename MultiLabelAnnotation<IndexType, LabelType>::Index;
    using Label = typename MultiLabelAnnotation<IndexType, LabelType>::Label;
    using VLabels = typename MultiLabelAnnotation<IndexType, LabelType>::VLabels;

    virtual ~MultiLabelEncoded() {}

    /******************* General functionality *******************/

    virtual void set_labels(Index i, const VLabels &labels) override = 0;
    virtual VLabels get_labels(Index i) const override = 0;

    virtual void add_label(Index i, const Label &label) override = 0;
    virtual void add_labels(Index i, const VLabels &labels) override = 0;
    virtual void add_labels(const std::vector<Index> &indices,
                            const VLabels &labels) override = 0;

    virtual bool has_label(Index i, const Label &label) const override = 0;
    virtual bool has_labels(Index i, const VLabels &labels) const override = 0;

    virtual void serialize(const std::string &filename) const override = 0;
    virtual bool merge_load(const std::vector<std::string> &filenames) override = 0;

    virtual void insert_rows(const std::vector<Index> &rows) override = 0;

    /*********************** Special queries **********************/

    // Get labels that occur at least in |presence_ratio| rows.
    // If |presence_ratio| = 0, return all occurring labels.
    virtual VLabels get_labels(const std::vector<Index> &indices,
                               double presence_ratio) const override = 0;

    // Count all labels collected from the given rows
    // and return top |num_top| with the their counts.
    virtual std::vector<std::pair<Label, size_t>>
    get_top_labels(const std::vector<Index> &indices,
                   size_t num_top = static_cast<size_t>(-1)) const override final;

    /************************* Properties *************************/

    virtual uint64_t num_objects() const override = 0;
    virtual size_t num_labels() const override = 0;
    virtual uint64_t num_relations() const override = 0;

  protected:
    virtual std::vector<uint64_t>
    count_labels(const std::vector<Index> &indices) const = 0;

    LabelEncoder<Label> label_encoder_;
};

} // namespace annotate

#endif // __ANNOTATE_HPP__
