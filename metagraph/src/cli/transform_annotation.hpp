#ifndef __TRANSFORM_ANNOTATION_HPP__
#define __TRANSFORM_ANNOTATION_HPP__


namespace mtg {
namespace cli {

class Config;

int transform_annotation(Config *config);

int merge_annotation(Config *config);

int relax_multi_brwt(Config *config);

} // namespace cli
} // namespace mtg

#endif // __TRANSFORM_ANNOTATION_HPP__
