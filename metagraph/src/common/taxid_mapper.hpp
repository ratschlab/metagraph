#ifndef __TAXID_MAPPER_HPP__
#define __TAXID_MAPPER_HPP__


#include <fstream>
#include <string>
#include <unordered_map>
#include <tsl/hopscotch_map.h>

class TaxIDMapper {
  public:
    typedef uint64_t taxid_t;

    bool parse_accession2taxid(const std::string &accession2taxid);
    bool parse_nodes(const std::string &nodes);
    bool parse_catalog(const std::string &catalog);
    static std::string parse_label(const std::string &gb);

    taxid_t gb_to_taxid(const std::string &gb) const;

    taxid_t get_parent(const taxid_t &id) const;

    std::string get_rank_label(const taxid_t &id) const;

    taxid_t get_ancestor_with_rank_label(taxid_t id, const std::string &rank_label) const;

    std::pair<taxid_t, std::string> get_ancestor_labeled(taxid_t id) const;

    bool load(std::ifstream &in);
    void serialize(std::ofstream &out) const;

  private:
    tsl::hopscotch_map<std::string, taxid_t> gb_to_taxid_;
    std::unordered_map<taxid_t, taxid_t> taxid_to_parent_;
    tsl::hopscotch_map<taxid_t, std::string> taxid_to_rank_label_;
};


#endif // __TAXID_MAPPER_HPP__
