#ifndef __METAGRAPH_VCFPARSER__
#define __METAGRAPH_VCFPARSER__


#include <htslib/vcf_sweep.h>
#include <htslib/faidx.h>
#include <string>
#include <vector>
#include <functional>
#include <iostream>


namespace mtg {
namespace seq_io {

class VCFParser {
  public:
    // reference_file -- file name of the FASTA file
    // vcf_file       -- file name of the VCF file
    VCFParser(const std::string &reference_file,
              const std::string &vcf_file,
              int k);

    ~VCFParser();

    void print_line(std::ostream &out = std::cout);

    void call_sequences(std::function<void(std::string&&)> callback,
                        std::function<bool()> terminate = []() { return false; });

    void call_annotated_sequences(std::function<void(std::string&&,
                                                     const std::vector<std::string>&)> callback,
                                  const std::vector<std::string> &annots,
                                  std::function<bool()> terminate = []() { return false; });

  private:
    faidx_t *reference_;
    bcf_sweep_t *sw_;
    uint32_t k_;

    std::string prefix_;
    std::string suffix_;

    bcf_hdr_t *hdr_;
    bcf1_t *rec_ = NULL;

    size_t current_line_ = 0;

    std::vector<const char *> seq_names_;

    std::vector<std::string> annotation_;

    bool read_next_line();

    void parse_annotation_from_vcf(const std::vector<std::string> &annots);
};

} // namespace seq_io
} // namespace mtg

#endif // __METAGRAPH_VCFPARSER__
