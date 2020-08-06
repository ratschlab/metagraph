#pragma once


class DiffAnnotation {
  private:
    sdsl::enc_vector<> diffs_;
    sdsl::bit_vector boundary_;
    sdsl::bit_vector terminal_;
    sdsl::select_support_mcl<> sboundary;

  public:
    const sdsl::bit_vector& terminal() const { return terminal_; }
    const sdsl::bit_vector& bundary() const { return boundary_; }
    const sdsl::enc_vector<>& diffs() const { return diffs_; }

    DiffAnnotation(const std::vector<uint64_t> &diffs,
                   const sdsl::bit_vector &boundary,
                   const sdsl::bit_vector &terminal)
            : diffs_(diffs), boundary_(boundary), terminal_(terminal) {
        sdsl::util::init_support(sboundary, &boundary);
    }

    std::vector<uint64_t> get_diff(uint64_t node_id) const {
        std::vector<uint64_t> result;
        uint64_t start_idx = sboundary.select(node_id) + 1;
        assert(boundary_[boundary_.size() - 1] == 1);
        while (boundary_[start_idx] == 0) {
            result.push_back(diffs_[start_idx]);
            start_idx++;
        }
        return result;
    }

    void serialize(const std::string &name) const {
        std::ofstream f(name);

        diffs_.serialize(f);
        boundary_.serialize(f);
        terminal_.serialize(f);

        f.close();
    }

    void load(const std::string &name) {
        std::ifstream f(name);

        diffs_.load(f);
        boundary_.load(f);
        terminal_.load(f);

        f.close();
    }
};
