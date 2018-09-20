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


// Terminology
//
// An annotated graph is a graph with labeled edges.
//
// A labeled edge is an edge carrying a corresponding label.
//
// If the edge labels represent inclusion of the edges into
// some edge categories, each of these categories is called
// an edge color, and the edge labels are called the edge colorings.
template <typename IndexType, typename ColorType>
class MultiColorAnnotation
      : public AnnotationCategory<IndexType, std::vector<ColorType>> {

  public:
    typedef IndexType Index;
    typedef ColorType Color;
    typedef std::vector<ColorType> Coloring;

    virtual ~MultiColorAnnotation() {}

    /***************** Inherited member functions ****************/

    virtual std::vector<Color> get(Index i) const override final {
        return get_coloring(i);
    }
    virtual void set(Index i, const std::vector<Color> &label) override final {
        set_coloring(i, label);
    }
    virtual void add(Index i, const std::vector<Color> &label) override final {
        add_colors(i, label);
    }

    /******************* General functionality *******************/

    virtual void set_coloring(Index i, const Coloring &coloring) = 0;
    virtual Coloring get_coloring(Index i) const = 0;

    virtual void add_color(Index i, const Color &color) = 0;
    virtual void add_colors(Index i, const Coloring &coloring) = 0;
    virtual void add_colors(const std::vector<Index> &indices,
                            const Coloring &coloring) = 0;

    virtual bool has_color(Index i, const Color &color) const = 0;
    virtual bool has_colors(Index i, const Coloring &coloring) const = 0;

    virtual void serialize(const std::string &filename) const = 0;
    virtual bool merge_load(const std::vector<std::string> &filenames) = 0;

    virtual void insert_rows(const std::vector<Index> &rows) = 0;

    /*********************** Special queries **********************/

    // Get colors that occur at least in |discovery_ratio| colorings.
    // If |discovery_ratio| = 0, return the union of colorings.
    virtual Coloring aggregate_colors(const std::vector<Index> &indices,
                                      double discovery_ratio = 1) const = 0;

    // Count all colors collected from extracted colorings
    // and return top |num_top| with the counts computed.
    virtual std::vector<std::pair<Color, size_t>>
    get_most_frequent_colors(const std::vector<Index> &indices,
                             size_t num_top = static_cast<size_t>(-1)) const = 0;

    /************************* Properties *************************/

    virtual size_t num_colors() const = 0;
    virtual double sparsity() const = 0;
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

    void serialize(std::ostream &outstream) const;
    bool load(std::istream &instream);

  private:
    std::unordered_map<Label, uint32_t> encode_label_;
    std::vector<Label> decode_label_;
};

} // namespace annotate

#endif // __ANNOTATE_HPP__
