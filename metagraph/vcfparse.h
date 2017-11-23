#ifndef METAGRAPH_VCFPARSE
#define METAGRAPH_VCFPARSE

#include <stdio.h>
#include <stdint.h>
#include <htslib/vcf_sweep.h>
#include <htslib/faidx.h>
#include <string>
#include <iostream>

#include "kseq.h"


static char passfilt[] = "PASS";


class vcf_parser {
  public:
    ~vcf_parser() {
        if (kmer1)
            vcf_clean_kmers();
        if (seq.s)
            free(seq.s);
        if (curalt)
            free(curalt);
        bcf_sweep_destroy(sw);
        fai_destroy(ref);
        free(seqnames);
    }

    bool init(const char *ref_, const char *vcf, int k_) {
        //load ref
        ref = fai_load(ref_);
        if (!ref) {
            fprintf(stderr, "Failed to read reference\n");
            exit(1);
        }
        //load vcf
        sw = bcf_sweep_init(vcf);
        if (!sw) {
            fprintf(stderr, "Failed to read VCF\n");
            exit(1);
        }
        k = k_;
        seq.s = (char*)calloc(1, sizeof(char));
        seq.l = 0;
        hdr = bcf_sweep_hdr(sw);
        kmer1 = NULL;
        rec = NULL;
        curalt = (char*)calloc(1, sizeof(char));
        ref_callele_l = 0;
        curaltlen = 0;
        if (!vcf_next_line()) {
            fprintf(stderr, "Empty VCF file\n");
            exit(1);
        }
        seqnames = bcf_hdr_seqnames(hdr, &(nseq));
        return true;
    }

    void vcf_print_line() {
        if (rec) {
            fprintf(stdout, "%s\t%d\t%d\t%d\t%d",
                    hdr->id[BCF_DT_CTG][rec->rid].key, //contig name
                    hdr->id[BCF_DT_CTG][rec->rid].val->info[0], //contig length
                    rec->rid, //contig id
                    rec->pos + 1, //1-based coordinate
                    rec->n_allele //number of alleles
                    );
        }
    }

    //TODO: make a parallel version of this that outputs n lines at a time?
    std::string vcf_get_seq(const char *annots[], size_t num_annots) {
        std::string annot = "";
        //std::string annot = "VCF:";
        //uint64_t annot;
        while (rec) {
            (curi)++;
            if ((curi >= rec->n_allele) || (!bcf_has_filter(hdr, rec, passfilt))) {
                vcf_next_line();
                continue;
            }
            /*
            if (seq.s) {
                free(seq.s);
                seq.s = NULL;
            }
            */
            char *temp;
            if (rec->d.allele[curi][0]=='<') {
                //if this is of the form <CN#>, then it's a copy number variation
                //otherwise, crash
                if (rec->d.allele[curi][1] != 'C' || rec->d.allele[curi][2] != 'N') {
                    //TODO: HANDLE RETROTRANSPOSONS PROPERLY
                    fprintf(stderr, "Can't handle this type of variant, skipping: %s\n",rec->d.allele[curi]);
                    vcf_next_line();
                    continue;
                    //exit(1);
                }
                //replace <CN#> with # copies of the ref allele
                //fprintf(stderr, "%s\n", rec->d.allele[curi]);
                rec->d.allele[curi][strlen(rec->d.allele[curi])-1]=0;
                //fprintf(stderr, "%s\n", rec->d.allele[curi]);
                size_t cn = atol(rec->d.allele[curi]+3);
                /*
                if (curalt) {
                    free(curalt);
                }
                */
                if (ref_callele_l * cn > curaltlen) {
                    temp = (char*)realloc(curalt, ref_callele_l * cn+1);
                    if (temp) {
                        curalt = temp;
                    } else {
                        free(curalt);
                        curalt = (char*)malloc(ref_callele_l * cn+1);
                    }
                }
                curalt[0] = 0;
                while (cn--) {
                    strcat(curalt, rec->d.allele[0]);
                }
                //fprintf(stderr,"\n%u\n%s\n%s\n",cn,rec->d.allele[0], curalt);
            } else {
                if (strlen(rec->d.allele[curi]) > curaltlen) {
                    temp = curalt = (char*)realloc(curalt, strlen(rec->d.allele[curi])+1);
                    if (temp) {
                        curalt = temp;
                    } else {
                        free(curalt);
                        curalt = (char*)malloc(strlen(rec->d.allele[curi])+1);
                    }
                }
                strcpy(curalt, rec->d.allele[curi]);
            }
            //curaltlen = strlen(rec->d.allele[curi]);
            curaltlen = strlen(curalt);
            if (kmer1_l + curaltlen + kmer3_l > seq.l) {
                temp = (char*)realloc(seq.s, kmer1_l + curaltlen + kmer3_l + 1);
                if (temp) {
                    seq.s = temp;
                } else {
                    free(seq.s);
                    seq.s = (char*)malloc(kmer1_l + curaltlen + kmer3_l + 1);
                }
            }
            //annotation is a bit vector indicating inclusion in the different ethnic groups
            //the first bit is always 1. if not, then the file is done and no sequence was output
            //TODO: check if these annots are part of the INFO, if not, check genotypes
            //ngt = bcf_get_genotypes(sr->readers[0].header, line0, &gt_arr, &ngt_arr); //get genotypes
            //annot=1;
            annot = std::string(seqnames[rec->rid]) + ":";
            for (size_t i=0; i < num_annots; ++i) {
                bcf_info_t *curinfo = bcf_get_info(hdr, rec, annots[i]);
                if (curinfo && curinfo->v1.i) {
                    annot += std::string(annots[i]) + ":";
                    //strcat(annot, annots[i]);
                    //strcat(annot, ":");
                    //fprintf(stdout, "%s\n", annot.c_str());
                }
                //annot[strlen(annot)-1]=0;
                //annot += (bool)(curinfo ? curinfo->v1.i : 0) * (1<<(i+1));
            }

            int32_t *gt_arr=NULL, ngt_arr = 0;
            int nsmpl = bcf_hdr_nsamples(hdr);
            bcf_fmt_t *curfmt;
            if ((curfmt = bcf_get_fmt(hdr, rec, "GT"))) {
                int ngt = bcf_get_genotypes(hdr, rec, &gt_arr, &ngt_arr);
                int ploidy = ngt/nsmpl;
                int cur;
                if (ngt > 0) {
                    for (int i = 0; i < ngt; i += ploidy) {
                        cur = 0;
                        for (int j = 0; j < ploidy; ++j) {
                            cur |= *(gt_arr + i + j);
                        }
                        if ((cur & 5) == 5) {
                            //at least one parent has this allele
                            annot += std::string(hdr->samples[i / ploidy]) + ":";
                        }
                    }
                }
            if (annot.length())
                free(gt_arr);
            }
            //if (annot[annot.length()-1] == ':')
            //if annot is of any length, assume that the last character is ':'
            if (annot.length())
                annot.pop_back();
            //std::cout << annot << "\n";
            strcpy(seq.s, kmer1);
            //strcat(seq.s, rec->d.allele[curi]);
            strcat(seq.s, curalt);
            strcat(seq.s, kmer3);
            seq.l = strlen(seq.s);
            //fprintf(stderr, "%s\n", seq.s);
            return annot.c_str();
        }
        return "";
    }

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
    void vcf_clean_kmers() {
        if (kmer1) {
            free(kmer1);
            kmer1 = NULL;
        }
        if (kmer3) {
            free(kmer3);
            kmer3 = NULL;
        }
        /*
        if (seq.s) {
            free(seq.s);
            seq.s = NULL;
        }
        */
        /*
        if (curalt) {
            free(curalt);
            curalt = NULL;
        }
        */
    }

    int vcf_next_line() {
        int len;
        if ( (rec = bcf_sweep_fwd(sw)) ) {
            bcf_unpack(rec, BCF_UN_FLT);
            if (kmer1)
                vcf_clean_kmers();
            curi = 0;
            ref_callele_l = strlen(rec->d.allele[0]);
            curkey = hdr->id[BCF_DT_CTG][rec->rid].key;
            curpos = &(rec->pos);
            kmer1 = faidx_fetch_seq(ref, hdr->id[BCF_DT_CTG][rec->rid].key,
                                         rec->pos - k,
                                         rec->pos - 1, &len);
            kmer3 = faidx_fetch_seq(ref, hdr->id[BCF_DT_CTG][rec->rid].key,
                                         rec->pos + ref_callele_l,
                                         rec->pos + ref_callele_l - 1 + k, &len);
            kmer1_l = strlen(kmer1);
            kmer3_l = strlen(kmer3);
            return 1;
        }
        return 0;
    }
};

#endif
