#include <string>
#include <vector>

namespace models {

typedef struct 
{
    std::string name;
} Tag;

typedef struct
{
    std::vector<Tag> tags;
} Annotation;

typedef struct
{
    std::string name;
} Kmer;

}
