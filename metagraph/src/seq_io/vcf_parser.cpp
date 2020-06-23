#include "vcf_parser.hpp"

#include <fstream>
#include <memory>
#include <filesystem>


namespace mtg {
namespace seq_io {

static char passfilt[] = "PASS";


VCFParser::~VCFParser() {
    bcf_sweep_destroy(sw_);
    fai_destroy(reference_);
}

VCFParser::VCFParser(const std::string &reference_file,
                     const std::string &vcf_file,
                     int k)
      : k_(k) {
    if (!std::filesystem::exists(reference_file))
        throw std::runtime_error("Failed to load reference file. Not found");

    if (!std::filesystem::exists(vcf_file))
        throw std::runtime_error("Failed to load VCF. Not found");

    if (!(reference_ = fai_load(reference_file.c_str())))
        throw std::ifstream::failure("Failed to read reference");

    if (!(sw_ = bcf_sweep_init(vcf_file.c_str())))
        throw std::ifstream::failure("Failed to initialize VCF");

    if (!(hdr_ = bcf_sweep_hdr(sw_)))
        throw std::runtime_error(
            "Failed to read VCF. If it is gzipped, recompress with bgzip"
        );

    read_next_line();

    int nseq = 0;

    // this is easier to read than the equivalent unique_ptr code
    auto **seq_names = bcf_hdr_seqnames(hdr_, &nseq);
    seq_names_.assign(seq_names, seq_names + nseq);
    free(seq_names);
}

void VCFParser::print_line(std::ostream &out) {
    if (rec_) {
        out << hdr_->id[BCF_DT_CTG][rec_->rid].key << "\t"          //contig name
            << hdr_->id[BCF_DT_CTG][rec_->rid].val->info[0] << "\t" //contig length
            << rec_->rid << "\t"                                    //contig id
            << rec_->pos + 1 << "\t"                                //1-based coordinate
            << rec_->n_allele;                                      //number of alleles
    }
}

void VCFParser::call_sequences(std::function<void(std::string&&)> callback,
                               std::function<bool()> terminate) {
    while (rec_ && !terminate()) {
        current_line_++;
        if (current_line_ >= rec_->n_allele || !bcf_has_filter(hdr_, rec_, passfilt)) {
            read_next_line();
            continue;
        }

        if (rec_->d.allele[current_line_][0] == '<') {
            std::string variant;
            //if this is of the form <CN#>, then it's a copy number variation
            //otherwise, crash
            if (rec_->d.allele[current_line_][1] != 'C'
                    || rec_->d.allele[current_line_][2] != 'N') {
                //TODO: HANDLE RETROTRANSPOSONS PROPERLY
                fprintf(stderr, "Can't handle this type of variant, skipping: %s\n",
                                rec_->d.allele[current_line_]);
                read_next_line();
                continue;
                //exit(1);
            }
            //replace <CN#> with # copies of the ref allele
            //fprintf(stderr, "%s\n", rec_->d.allele[current_line_]);
            rec_->d.allele[current_line_][strlen(rec_->d.allele[current_line_]) - 1] = 0;
            //fprintf(stderr, "%s\n", rec_->d.allele[curi]);
            size_t cn = atol(rec_->d.allele[current_line_] + 3);

            // TODO: better way to do this
            while (cn--) {
                variant += rec_->d.allele[0];
            }

            callback(prefix_ + variant + suffix_);
            //fprintf(stderr, "\n%u\n%s\n%s\n", cn, rec_->d.allele[0], variant);
        } else {
            callback(prefix_ + rec_->d.allele[current_line_] + suffix_);
        }
    }
}

void VCFParser
::call_annotated_sequences(std::function<void(std::string&&,
                                              const std::vector<std::string>&)> callback,
                           const std::vector<std::string> &annots,
                           std::function<bool()> terminate) {
    call_sequences(
        [&](auto&& sequence) {
            parse_annotation_from_vcf(annots);
            callback(std::move(sequence), annotation_);
        },
        terminate
    );
}

//annotation is a bit vector indicating inclusion in the different ethnic groups
//the first bit is always 1. if not, then the file is done and no sequence was output
//TODO: check if these annots are part of the INFO, if not, check genotypes
//ngt = bcf_get_genotypes(sr->readers[0].header, line0, &gt_arr, &ngt_arr); //get genotypes
//annot=1;
void VCFParser::parse_annotation_from_vcf(const std::vector<std::string> &annots) {
    annotation_.clear();
    annotation_.emplace_back(seq_names_.at(rec_->rid));
    for (const auto &annot : annots) {
        bcf_info_t *curinfo = bcf_get_info(hdr_, rec_, annot.c_str());
        if (curinfo && curinfo->v1.i) {
            annotation_.emplace_back(annot);
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
                annotation_.emplace_back(hdr_->samples[i]);
            }
        }
        if (gt_arr)
            free(gt_arr);
    }
}

bool VCFParser::read_next_line() {
    if (!(rec_ = bcf_sweep_fwd(sw_)))
        return false;

    current_line_ = 0;
    bcf_unpack(rec_, BCF_UN_FLT);

    uint32_t ref_callele_l = strlen(rec_->d.allele[0]);

    int len;
    auto first = faidx_fetch_seq(reference_,
                                 hdr_->id[BCF_DT_CTG][rec_->rid].key,
                                 rec_->pos - k_,
                                 rec_->pos - 1,
                                 &len);
    prefix_.assign(first, len);
    free(first);

    auto second = faidx_fetch_seq(reference_,
                                  hdr_->id[BCF_DT_CTG][rec_->rid].key,
                                  rec_->pos + ref_callele_l,
                                  rec_->pos + ref_callele_l - 1 + k_,
                                  &len);
    suffix_.assign(second, len);
    free(second);

    return true;
}

} // namespace seq_io
} // namespace mtg
