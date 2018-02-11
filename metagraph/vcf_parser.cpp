#include "vcf_parser.hpp"


static char passfilt[] = "PASS";


vcf_parser::~vcf_parser() {
    if (sw_)
        bcf_sweep_destroy(sw_);
    if (reference_)
        fai_destroy(reference_);
}

bool vcf_parser::init(const std::string &reference_file,
                      const std::string &vcf_file, int k) {
    assert(!reference_ && "Not initialized");

    if (!(reference_ = fai_load(reference_file.c_str()))) {
        fprintf(stderr, "Failed to read reference\n");
        return false;
    }

    if (!(sw_ = bcf_sweep_init(vcf_file.c_str()))) {
        fprintf(stderr, "Failed to read VCF\n");
        return false;
    }

    k_ = k;

    hdr_ = bcf_sweep_hdr(sw_);

    if (!read_next_line()) {
        fprintf(stderr, "Empty VCF file\n");
        return false;
    }

    int nseq = 0;
    auto **seq_names = bcf_hdr_seqnames(hdr_, &nseq);
    seq_names_.assign(seq_names, seq_names + nseq);
    free(seq_names);

    return true;
}

void vcf_parser::print_line() {
    if (!rec_)
        return;

    fprintf(stdout, "%s\t%d\t%d\t%d\t%d",
            hdr_->id[BCF_DT_CTG][rec_->rid].key, //contig name
            hdr_->id[BCF_DT_CTG][rec_->rid].val->info[0], //contig length
            rec_->rid, //contig id
            rec_->pos + 1, //1-based coordinate
            rec_->n_allele //number of alleles
    );
}

bool vcf_parser::get_seq(const std::vector<std::string> &annots,
                         std::string *sequence,
                         std::vector<std::string> &annotation) {
    while (rec_) {
        curi++;
        if (curi >= rec_->n_allele || !bcf_has_filter(hdr_, rec_, passfilt)) {
            read_next_line();
            continue;
        }
        std::string curalt;
        if (rec_->d.allele[curi][0] == '<') {
            //if this is of the form <CN#>, then it's a copy number variation
            //otherwise, crash
            if (rec_->d.allele[curi][1] != 'C' || rec_->d.allele[curi][2] != 'N') {
                //TODO: HANDLE RETROTRANSPOSONS PROPERLY
                fprintf(stderr, "Can't handle this type of variant, skipping: %s\n",
                                rec_->d.allele[curi]);
                read_next_line();
                continue;
                //exit(1);
            }
            //replace <CN#> with # copies of the ref allele
            //fprintf(stderr, "%s\n", rec_->d.allele[curi]);
            rec_->d.allele[curi][strlen(rec_->d.allele[curi]) - 1] = 0;
            //fprintf(stderr, "%s\n", rec_->d.allele[curi]);
            size_t cn = atol(rec_->d.allele[curi] + 3);

            while (cn--) {
                curalt += rec_->d.allele[0];
            }
            //fprintf(stderr, "\n%u\n%s\n%s\n", cn, rec_->d.allele[0], curalt);
        } else {
            curalt = rec_->d.allele[curi];
        }
        //annotation is a bit vector indicating inclusion in the different ethnic groups
        //the first bit is always 1. if not, then the file is done and no sequence was output
        //TODO: check if these annots are part of the INFO, if not, check genotypes
        //ngt = bcf_get_genotypes(sr->readers[0].header, line0, &gt_arr, &ngt_arr); //get genotypes
        //annot=1;
        //annotation->append(seq_names_[rec_->rid]);
        annotation.emplace_back(seq_names_[rec_->rid]);
        for (const auto &annot : annots) {
            bcf_info_t *curinfo = bcf_get_info(hdr_, rec_, annot.c_str());
            if (curinfo && curinfo->v1.i) {
                annotation.emplace_back(annot);
            }
        }

        if ((bcf_get_fmt(hdr_, rec_, "GT"))) {
            int32_t *gt_arr = NULL, ngt_arr = 0;
            int ngt = bcf_get_genotypes(hdr_, rec_, &gt_arr, &ngt_arr);
            int nsmpl = bcf_hdr_nsamples(hdr_);
            int ploidy = ngt / nsmpl;
            for (int i = 0; i < nsmpl; ++i) {
                int cur = 0;
                for (int j = 0; j < ploidy; ++j) {
                    cur |= gt_arr[i * ploidy + j];
                }
                if ((cur & 5) == 5) {
                    //at least one parent has this allele
                    annotation.emplace_back(hdr_->samples[i]);
                }
            }
            if (gt_arr)
                free(gt_arr);
        }

        *sequence = kmer1_ + curalt + kmer3_;
        return true;
    }
    return false;
}

bool vcf_parser::read_next_line() {
    if (!(rec_ = bcf_sweep_fwd(sw_)))
        return false;

    curi = 0;
    bcf_unpack(rec_, BCF_UN_FLT);

    uint32_t ref_callele_l = strlen(rec_->d.allele[0]);
    int len;
    auto first = faidx_fetch_seq(reference_,
                               hdr_->id[BCF_DT_CTG][rec_->rid].key,
                               rec_->pos - k_,
                               rec_->pos - 1, &len);
    kmer1_.assign(first, len);
    free(first);
    auto second = faidx_fetch_seq(reference_,
                               hdr_->id[BCF_DT_CTG][rec_->rid].key,
                               rec_->pos + ref_callele_l,
                               rec_->pos + ref_callele_l - 1 + k_, &len);
    kmer3_.assign(second, len);
    free(second);
    return true;
}
