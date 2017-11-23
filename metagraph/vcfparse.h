#ifndef METAGRAPH_VCFPARSE
#define METAGRAPH_VCFPARSE

#include <htslib/vcf_sweep.h>
#include <htslib/faidx.h>
#include <string>


class vcf_parser {
  public:
    ~vcf_parser();

    bool init(const char *ref_, const char *vcf, int k_);

    void vcf_print_line();

    //TODO: make a parallel version of this that outputs n lines at a time?
    std::string vcf_get_seq(const char *annots[], size_t num_annots);

    faidx_t *ref;
    bcf_sweep_t *sw;
    bcf_hdr_t *hdr;
    bcf1_t *rec;
    uint32_t k;
    uint32_t curi;
    int32_t *curpos;
    const char *curkey;
    char *curalt;
    uint32_t curaltlen;
    kstring_t seq;
    char *kmer1;
    char *kmer3;
    uint32_t kmer1_l;
    uint32_t kmer3_l;
    uint32_t ref_callele_l;
    const char **seqnames = NULL;
    int nseq = 0;

  private:
    void vcf_clean_kmers();

    int vcf_next_line();
};

#endif
