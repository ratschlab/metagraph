#ifndef __FILE_UTILS_HPP__
#define __FILE_UTILS_HPP__

#include <string>
#include <memory>
#include <fstream>


namespace utils {

bool check_if_writable(const std::string &filename);

/**
 *  Returns the input file type, given a filename
 *  One of: [KMC|VCF|FASTQ|FASTA].
 *  If filetype is unknown, return empty string "".
 */
std::string get_filetype(const std::string &fname);

class TempFile {
  public:
    // The state flow:
    //    init -> APPEND -> READ -> deinit
    enum State { APPEND, READ };

    TempFile(const std::string &tmp_dir = "");
    ~TempFile();

    std::ofstream& ofstream();
    std::ifstream& ifstream();

  private:
    std::string tmp_file_name_;
    State state_;
    std::unique_ptr<std::ofstream> tmp_ostream_;
    std::unique_ptr<std::ifstream> tmp_istream_;
};

} // namespace utils

#endif // __FILE_UTILS_HPP__
