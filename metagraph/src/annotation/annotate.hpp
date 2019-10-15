#ifndef __ANNOTATE_HPP__
#define __ANNOTATE_HPP__

#include <cstdint>
#include <string>
#include <vector>
#include <unordered_map>
#include <memory>
#include <functional>


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

    virtual void insert_rows(const std::vector<Index> &rows) = 0;

    // For each pair (L, L') in the dictionary, replaces label |L| with |L'|
    // and merges all relations (*, L') with matching labels L', if supported.
    virtual void rename_labels(const std::unordered_map<Label, Label> &dict) = 0;

    /*********************** Special queries **********************/

    virtual void call_objects(const Label &label,
                              std::function<void(Index)> callback) const = 0;

    /************************* Properties *************************/

    virtual uint64_t num_objects() const = 0;
    virtual size_t num_labels() const = 0;
    virtual uint64_t num_relations() const = 0;
    virtual const VLabels& get_all_labels() const = 0;

    virtual bool label_exists(const Label &label) const = 0;

    virtual std::string file_extension() const = 0;
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
     * Check if the label has been added to encoder.
     */
    bool label_exists(const Label &label) const { return encode_label_.count(label); }

    /**
     * Throws an exception if a bad code is passed.
     */
    const Label& decode(size_t code) const { return decode_label_.at(code); }

    const std::vector<Label>& get_labels() const { return decode_label_; }

    size_t size() const { return decode_label_.size(); }

    bool load(std::istream &instream);
    void serialize(std::ostream &outstream) const;

    void merge(const LabelEncoder<Label> &other);

    void clear() { encode_label_.clear(); decode_label_.clear(); }

  private:
    std::unordered_map<Label, uint64_t> encode_label_;
    std::vector<Label> decode_label_;
};


class IterateRows {
  public:
    virtual ~IterateRows() {}
    virtual std::vector<uint64_t> next_row() = 0;
};


template <typename IndexType, typename LabelType>
class MultiLabelEncoded
      : public MultiLabelAnnotation<IndexType, LabelType> {
  public:
    using Index = typename MultiLabelAnnotation<IndexType, LabelType>::Index;
    using Label = typename MultiLabelAnnotation<IndexType, LabelType>::Label;
    using VLabels = typename MultiLabelAnnotation<IndexType, LabelType>::VLabels;

    virtual ~MultiLabelEncoded() {}

    virtual VLabels get_labels(Index i) const override final;

    virtual std::unique_ptr<IterateRows> iterator() const;
    virtual std::vector<uint64_t> get_label_codes(Index i) const = 0;

    virtual const LabelEncoder<Label>& get_label_encoder() const final { return label_encoder_; }

    /******************* General functionality *******************/

    // For each pair (L, L') in the dictionary, replaces label |L| with |L'|
    // and merges all relations (*, L') with matching labels L', if supported.
    virtual void rename_labels(const std::unordered_map<Label, Label> &dict) override;

    // For each Index in indices, call row_callback on the vector of its
    // corresponding label indices. Terminate early if terminate returns true.
    virtual void call_rows(const std::vector<Index> &indices,
                           const std::function<void(std::vector<uint64_t>&&)> &row_callback,
                           const std::function<bool()> &terminate = []() { return false; }) const;

    virtual bool label_exists(const Label &label) const override final {
        return label_encoder_.label_exists(label);
    }

    virtual const VLabels& get_all_labels() const override final {
        return label_encoder_.get_labels();
    }

    virtual std::string file_extension() const override = 0;

  protected:
    LabelEncoder<Label> label_encoder_;
};

} // namespace annotate

#endif // __ANNOTATE_HPP__
