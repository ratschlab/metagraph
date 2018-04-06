#ifndef __ANNOTATE_HPP__
#define __ANNOTATE_HPP__

#include <cstdint>
#include <string>
#include <vector>

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
