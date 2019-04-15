#ifndef __ANNOTATION_CONVERTERS_HPP__
#define __ANNOTATION_CONVERTERS_HPP__

#include <memory>
#include <vector>
#include <functional>

#include "annotate.hpp"


namespace annotate {

typedef LabelEncoder<std::string> LEncoder;


template <typename Label>
class RowCompressed;

template <class StaticAnnotation, typename Label>
typename std::unique_ptr<StaticAnnotation>
convert(RowCompressed<Label>&& annotation);

template <class FromAnnotation, class ToAnnotation, typename Label, bool sparse = false>
typename std::unique_ptr<ToAnnotation>
convert(const std::string &filename);


template <class FromAnnotation, class ToAnnotation, typename Label, bool sparse = false>
uint64_t merge(const std::vector<std::string> &filenames, const std::string &outfile);

template <class FromAnnotation, class ToAnnotation, typename Label, bool sparse = false>
uint64_t merge(const std::vector<const FromAnnotation*> &annotators, const std::string &outfile);


template <typename Label>
class ColumnCompressed;

template <class StaticAnnotation, typename Label>
typename std::unique_ptr<StaticAnnotation>
convert(ColumnCompressed<Label>&& annotation);

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

} // namespace annotate

#endif // __ANNOTATION_CONVERTERS_HPP__
