#ifndef METAGRAPH_VCFPARSE
#define METAGRAPH_VCFPARSE

#include <htslib/vcf_sweep.h>
#include <htslib/faidx.h>
#include <string>
#include <vector>


class vcf_parser {
  public:
    ~vcf_parser();

    bool init(const std::string &reference_file,
              const std::string &vcf_file, int k);

    void vcf_print_line();

    //TODO: make a parallel version of this that outputs n lines at a time?
    std::string vcf_get_seq(const char *annots[], size_t num_annots);

    bcf_sweep_t *sw;
    bcf_hdr_t *hdr;
    bcf1_t *rec;
    uint32_t curi;
    char *curalt;
    uint32_t curaltlen;
    kstring_t seq;
    char *kmer1;
    char *kmer3;
    uint32_t kmer1_l;
    uint32_t kmer3_l;
    uint32_t ref_callele_l;

  private:
    uint32_t k_;
    faidx_t *reference_;
    std::vector<const char *> seq_names_;

    void vcf_clean_kmers();
    int vcf_next_line();
};

#endif
