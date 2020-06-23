#ifndef __SEQ_IO_FORMATS__
#define __SEQ_IO_FORMATS__

#include <string>


namespace mtg {
namespace seq_io {

/**
 *  Returns the input file type, given a filename
 *  One of: [KMC|VCF|FASTQ|FASTA].
 *  If filetype is unknown, return empty string "".
 */
std::string file_format(const std::string &fname);

} // namespace seq_io
} // namespace mtg

#endif // __SEQ_IO_FORMATS__
