#ifndef __ANNOTATE_HPP__
#define __ANNOTATE_HPP__

#include <cstdint>
#include <string>
#include <vector>
#include <unordered_map>


namespace annotate {

typedef uint64_t Index;


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


template <typename LabelType>
class AnnotationCategory {
  public:
    virtual ~AnnotationCategory() {}
    virtual LabelType get(Index i) const = 0;
    virtual void set_label(Index i, const LabelType &label) = 0;

    virtual void add_labels(const std::string &sequence, const LabelType &labels) = 0;

    virtual bool has_label(Index i, const LabelType &label) const = 0;

    virtual std::vector<std::string> get_label_names() const;

    virtual bool load(const std::string &filename) = 0;
    virtual void serialize(const std::string &filename) const = 0;
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
template <typename Index, typename Color>
class MultiColorAnnotation {
  public:
    typedef std::vector<Color> Coloring;

    virtual ~MultiColorAnnotation() {}

    // General functionality
    virtual Coloring get_coloring(Index i) const = 0;
    virtual std::vector<Coloring>
    get_colorings(const std::vector<Index> &indices) const = 0;

    virtual void add_color(Index i, const Color &color) = 0;
    virtual void add_colors(Index i, const Coloring &coloring) = 0;
    virtual void add_colors(const std::vector<Index> &indices,
                            const Coloring &coloring) = 0;

    virtual bool is_colored(Index i, const Color &color) const = 0;

    virtual bool load(const std::string &filename) = 0;
    virtual void serialize(const std::string &filename) const = 0;

    // Special queries

    // Compute the union of colorings excluding colors
    // observed in less than |filtering_ratio| colorings.
    virtual Coloring aggregate_colors(const std::vector<Index> &indices,
                                      double filtering_ratio) const = 0;

    // Count all colors collected from extracted colorings
    // and return top |num_top| with the counts computed.
    virtual std::unordered_map<Color, size_t>
    get_most_frequent_colors(const std::vector<Index> &indices,
                             size_t num_top) const = 0;
};


/*
class EdgeWiseMatrix : public UncompressedMatrix {
  public:
};
*/


/*
class WaveletTrie : public GenomeAnnotation {
  public:
};


class EdgeCompressed : public GenomeAnnotation {
  public:
};

*/


} // namespace annotate

#endif // __ANNOTATE_HPP__
