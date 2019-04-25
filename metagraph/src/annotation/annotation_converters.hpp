#ifndef __ANNOTATION_CONVERTERS_HPP__
#define __ANNOTATION_CONVERTERS_HPP__

#include <memory>
#include <vector>


namespace annotate {

template <typename Label>
class RowCompressed;

template <typename IndexType, typename LabelType>
class MultiLabelEncoded;


template <class StaticAnnotation, typename Label>
typename std::unique_ptr<StaticAnnotation>
convert(RowCompressed<Label>&& annotation);


template <typename Label>
class ColumnCompressed;

template <class StaticAnnotation, typename Label>
std::unique_ptr<StaticAnnotation> convert(ColumnCompressed<Label>&& annotation);

template <class StaticAnnotation>
std::unique_ptr<StaticAnnotation> convert(const std::string &filename);

template <class StaticAnnotation, typename Label>
typename std::unique_ptr<StaticAnnotation>
convert_to_simple_BRWT(ColumnCompressed<Label>&& annotation,
                       size_t grouping_arity = 2,
                       size_t num_threads = 1);

template <class StaticAnnotation, typename Label>
typename std::unique_ptr<StaticAnnotation>
convert_to_greedy_BRWT(ColumnCompressed<Label>&& annotation,
                       size_t num_threads = 1);

template <class StaticAnnotation>
void relax_BRWT(StaticAnnotation *annotation,
                size_t relax_max_arity,
                size_t num_threads = 1);


template <class ToAnnotation, typename Label = std::string>
void merge(const std::vector<const MultiLabelEncoded<uint64_t, Label>*> &annotators,
           const std::vector<std::string> &filenames,
           const std::string &outfile);

template <class ToAnnotation, class Annotation>
void merge(const std::vector<std::shared_ptr<Annotation>> &annotators,
           const std::vector<std::string> &filenames,
           const std::string &outfile) {
    std::vector<const MultiLabelEncoded<uint64_t, typename Annotation::Label>*> anno_vector;
    for (const auto &annotator : annotators) {
        anno_vector.push_back(annotator.get());
    }
    merge<ToAnnotation>(anno_vector, filenames, outfile);
}

template <class ToAnnotation, class Annotation>
void merge(const std::vector<std::unique_ptr<Annotation>> &annotators,
           const std::vector<std::string> &filenames,
           const std::string &outfile) {
    std::vector<const MultiLabelEncoded<uint64_t, typename Annotation::Label>*> anno_vector;
    for (const auto &annotator : annotators) {
        anno_vector.push_back(annotator.get());
    }
    merge<ToAnnotation>(anno_vector, filenames, outfile);
}

} // namespace annotate

#endif // __ANNOTATION_CONVERTERS_HPP__
