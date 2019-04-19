#ifndef __METAGRAPH_VCFPARSER__
#define __METAGRAPH_VCFPARSER__

#include <htslib/vcf_sweep.h>
#include <htslib/faidx.h>
#include <string>
#include <vector>


class vcf_parser {
  public:
    std::string seq;

    ~vcf_parser();

    // reference_file -- file name of the FASTA file
    // vcf_file       -- file name of the VCF file
    bool init(const std::string &reference_file,
              const std::string &vcf_file, int k);

    void print_line();

    //TODO: make a parallel version of this that outputs n lines at a time?
    bool get_seq(const std::vector<std::string> &annots,
                 std::vector<std::string> *annotation = NULL);

  private:
    faidx_t *reference_ = NULL;
    bcf_sweep_t *sw_ = NULL;
    uint32_t k_;

    std::string kmer1_;
    std::string kmer3_;

    bcf_hdr_t *hdr_;
    bcf1_t *rec_ = NULL;

    size_t curi = 0;

    std::vector<const char *> seq_names_;

    bool read_next_line();
};

#endif // __METAGRAPH_VCFPARSER__
