#include "dbg_sshash.hpp"
#include <dictionary.hpp>


namespace mtg {
namespace graph {
DBGSSHash::~DBGSSHash() {}

DBGSSHash::DBGSSHash(size_t k):k_(k) {
    dict_ = std::make_unique<sshash::dictionary>();
}

DBGSSHash::DBGSSHash(std::string const& input_filename, size_t k):k_(k){
    sshash::build_configuration build_config;
    build_config.k = k;//
    // quick fix for value of m... k/2 but odd
    build_config.m = (k_+1)/2;
    if(build_config.m % 2 == 0) build_config.m++;
    dict_ = std::make_unique<sshash::dictionary>();
    dict_->build(input_filename, build_config);
}
    std::string DBGSSHash::file_extension() const  { return kExtension; }
    size_t DBGSSHash::get_k() const  { return k_; }
    DeBruijnGraph::Mode DBGSSHash::get_mode() const  { return BASIC; }

void DBGSSHash::add_sequence(std::string_view sequence,
                             const std::function<void(node_index)> &on_insertion) {
    // TODO: throw exception? :)
        throw std::runtime_error("adding sequences not implemented");

}

void DBGSSHash::map_to_nodes(std::string_view sequence,
                             const std::function<void(node_index)> &callback,
                             const std::function<bool()> &terminate) const {
    map_to_nodes_sequentially(sequence, callback, terminate);
}

void DBGSSHash::map_to_nodes_sequentially(std::string_view sequence,
                                           const std::function<void(node_index)> &callback,
                                           const std::function<bool()> &terminate) const {
    for (size_t i = 0; i + k_ <= sequence.size() && !terminate(); ++i) {
	callback(kmer_to_node(sequence.substr(i, k_)));
    }
}

DBGSSHash::node_index DBGSSHash::traverse(node_index node, char next_char) const {
    std::string kmer = DBGSSHash::get_node_sequence(node); 
    sshash::neighbourhood nb = dict_->kmer_forward_neighbours(&kmer[0]);
    uint64_t ssh_idx = -1; 
    switch (next_char) {
    case 'A':
        ssh_idx = nb.forward_A.kmer_id;
        break;
    case 'C':
        ssh_idx = nb.forward_C.kmer_id;
        break;
    case 'G':
        ssh_idx = nb.forward_G.kmer_id;
        break;
    case 'T':
        ssh_idx = nb.forward_T.kmer_id;
        break;
    default:
        break;
    }
    return ssh_idx + 1;
}

DBGSSHash::node_index DBGSSHash::traverse_back(node_index node, char prev_char) const {
    std::string kmer = DBGSSHash::get_node_sequence(node);
    sshash::neighbourhood nb = dict_->kmer_backward_neighbours(&kmer[0]);
    uint64_t ssh_idx = -1; 
    switch (prev_char) {
    case 'A':
        ssh_idx = nb.backward_A.kmer_id;
        break;
    case 'C':
        ssh_idx = nb.backward_C.kmer_id;
        break;
    case 'G':
        ssh_idx = nb.backward_G.kmer_id;
        break;
    case 'T':
        ssh_idx = nb.backward_T.kmer_id;
        break;
    default:
        break;
    }
    return ssh_idx + 1;
}

void DBGSSHash ::adjacent_outgoing_nodes(node_index node,
                                         const std::function<void(node_index)> &callback) const {
    assert(node > 0 && node <= num_nodes());
    call_outgoing_kmers(node, [&](auto child, char) { callback(child); });
    
}

void DBGSSHash ::adjacent_incoming_nodes(node_index node,
                                         const std::function<void(node_index)> &callback) const {
    assert(node > 0 && node <= num_nodes());
    call_incoming_kmers(node, [&](auto parent, char) { callback(parent); });
}

void DBGSSHash ::call_outgoing_kmers(node_index node,
                                     const OutgoingEdgeCallback &callback) const {
    assert(node > 0 && node <= num_nodes());

    auto prefix = get_node_sequence(node).substr(1);

    for (char c : alphabet_) {
        auto next = kmer_to_node(prefix + c);
        if (next != npos)
            callback(next, c);
    }
}


void DBGSSHash ::call_incoming_kmers(node_index node,
                                     const IncomingEdgeCallback &callback) const {
    assert(node > 0 && node <= num_nodes());

    std::string suffix = get_node_sequence(node);
    suffix.pop_back();

    for (char c : alphabet_) {
        auto prev = kmer_to_node(c + suffix);
        if (prev != npos)
            callback(prev, c);
    }
}

size_t DBGSSHash::outdegree(node_index node) const {
    std::string kmer = DBGSSHash::get_node_sequence(node);
    sshash::neighbourhood nb = dict_->kmer_forward_neighbours(&kmer[0]);
    size_t out_deg = bool(nb.forward_A.kmer_id + 1) // change to loop?
                    + bool(nb.forward_C.kmer_id + 1)
                    + bool(nb.forward_G.kmer_id + 1)
                    + bool(nb.forward_T.kmer_id + 1);
    return out_deg;
}

bool DBGSSHash::has_single_outgoing(node_index node) const {
    return DBGSSHash::outdegree(node) == 1;
}

bool DBGSSHash::has_multiple_outgoing(node_index node) const {
    return DBGSSHash::outdegree(node) > 1;
}

size_t DBGSSHash::indegree(node_index node) const {
    std::string kmer = DBGSSHash::get_node_sequence(node);
    sshash::neighbourhood nb = dict_->kmer_backward_neighbours(&kmer[0]);
    size_t in_deg = bool(nb.backward_A.kmer_id + 1) // change to loop?
                    + bool(nb.backward_C.kmer_id + 1)
                    + bool(nb.backward_G.kmer_id + 1)
                    + bool(nb.backward_T.kmer_id + 1);
    return in_deg;
}

bool DBGSSHash::has_no_incoming(node_index node) const {
    return DBGSSHash::indegree(node) == 0;
}

bool DBGSSHash::has_single_incoming(node_index node) const {
    return DBGSSHash::indegree(node) == 1;
}

void DBGSSHash::call_kmers(
        const std::function<void(node_index, const std::string &)> &callback) const {
    for (size_t node_idx = 1; node_idx <= num_nodes(); ++node_idx) {
        callback(node_idx, get_node_sequence(node_idx));
    }
}

DBGSSHash::node_index DBGSSHash::kmer_to_node_from_superkmer(std::string_view kmer) const {
    uint64_t ssh_idx = dict_->lookup(kmer.begin(), true);
    return ssh_idx + 1;
}

DBGSSHash::node_index DBGSSHash::kmer_to_node(std::string_view kmer) const {
    uint64_t ssh_idx = dict_->lookup(kmer.begin(), true);
    return ssh_idx + 1;
}

// superkmer experiment: use minimizer to get superkmer positions (offsets) -> get superkmer that contains kmer 
std::pair<DBGSSHash::node_index, uint64_t> DBGSSHash::kmer_to_superkmer_node(std::string_view kmer) const {
    auto [ssh_idx, superkmer_id]  = dict_->kmer_to_superkmer_idx(kmer.begin(), true);
    if(ssh_idx == sshash::constants::invalid_uint64){
        if(dict_->lookup(kmer.begin(), true)!= sshash::constants::invalid_uint64){
            std::cout<< kmer <<"\n***********  NOT FOUND WITH KMER_TO_SUPERKMER_IDX BUT FOUND WITH LOOKUP" <<std::endl;
        }
        return std::pair(npos, sshash::constants::invalid_uint64);
    }
    if(dict_->lookup(kmer.begin(), true) == sshash::constants::invalid_uint64){
        std::cout<< kmer <<"\n***********  found with kmer_to_index but not with lookup!! *********" <<std::endl;
    }

    return std::pair(ssh_idx + 1, superkmer_id);
}

std::string DBGSSHash::get_node_sequence(node_index node) const {
    std::string str_kmer = "";
    str_kmer.append(k_, ' ');
    uint64_t ssh_idx = node - 1; // switch back to sshash idx!!!
    dict_->access(ssh_idx, &str_kmer[0]);
    return str_kmer;
}

void DBGSSHash::serialize(std::ostream &out) const {
    //TODO
    throw std::runtime_error("serialize to stream not implemented");
}

void DBGSSHash::serialize(const std::string &filename) const {
    std::string suffixed_filename = utils::make_suffix(filename, kExtension);
    
    essentials::logger("saving data structure to disk...");
    essentials::save(*dict_, suffixed_filename.c_str());
    essentials::logger("DONE");
}

bool DBGSSHash::load(std::istream &in) {
    throw std::runtime_error("load from stream not implemented");
    return false;
}

bool DBGSSHash::load(const std::string &filename) {
    std::string suffixed_filename = utils::make_suffix(filename, kExtension);
    uint64_t num_bytes_read = essentials::load(*dict_, suffixed_filename.c_str());
    bool verbose = true; // temp
    if (verbose) {
        std::cout << "index size: " << essentials::convert(num_bytes_read, essentials::MB)
                  << " [MB] (" << (num_bytes_read * 8.0) / dict_->size() << " [bits/kmer])"
                  << std::endl;
        dict_->print_info();
    }
    k_ = dict_->k();

    std::string s_mask_name = utils::remove_suffix(filename, kExtension) + "_sk_mask";
    std::cout << "LOADING MASK! \n";
    load_superkmer_mask(s_mask_name);
    return true;
}

bool DBGSSHash::operator==(const DeBruijnGraph &other) const {
    throw std::runtime_error("operator== not implemented");
    return false;
}

const std::string DBGSSHash::alphabet_ = "ACGT";
const std::string &DBGSSHash::alphabet() const {
    return alphabet_;
}

uint64_t DBGSSHash::num_nodes() const { return dict_->size(); }
void write_stats_to_file(const std::vector<uint64_t>& km_skmer, const std::vector<uint64_t>& cc_skmer){
    std::ofstream output_file_kmers("./superkmer_stats_kmers.txt");
    output_file_kmers << "kmers_per_superkmer: "<<'\n';

    for(auto const& x : km_skmer){
        output_file_kmers << x << '\n';
    }
    output_file_kmers.close();
    
    std::ofstream output_file_colors("./superkmer_stats_color.txt");
    output_file_colors << "color_changes_per_superkmer: "<<'\n';
    for(auto const& y : cc_skmer){
        output_file_colors << y << '\n';
    }
    output_file_colors.close();
}
sdsl::bit_vector mask_into_bit_vec(const std::vector<bool>& mask){
    sdsl::bit_vector bv (mask.size());
    for(size_t idx = 0; idx < mask.size(); idx++){
        bv[idx] = mask[idx];
    }
    return bv;
}

void DBGSSHash::load_superkmer_mask(std::string file){
    loaded_mask = load_from_file(superkmer_mask, file);
    std::cout<< " successfully loaded " << file<<"?: " <<loaded_mask << std::endl;
}

void DBGSSHash::superkmer_statistics(const std::unique_ptr<AnnotatedDBG>& anno_graph, std::string file_sk_mask) const{
    std::cout<< "Computing superkmer statistics and building super kmer bit vector... \n";   
    
    //getting labels in batches
    size_t num_labels = anno_graph->get_annotator().num_labels();
    std::vector<uint64_t> superkmer_idxs = dict_->build_superkmer_bv([&anno_graph, num_labels](std::string_view sequence){
                    auto labels = anno_graph->get_top_label_signatures(sequence, num_labels);
                    // since get_top_label_signatures returns only labels that were found, 
                    // if any entry is zero not all the kmers share all the same labels -> return false
                    for(auto pair : labels){
                        sdsl::rank_support_v rs;
                        sdsl::util::init_support(rs,&pair.second);
                        size_t bit_vec_size = pair.second.size();
                        if(rs(bit_vec_size) != bit_vec_size) return false;
                    }
                    return true;
         });
    //std::vector<bool> superkmer_mask = dict_->build_superkmer_bv([&anno_graph](std::string_view str){return anno_graph->get_labels(str);});

    //sdsl::bit_vector non_mono_superkmer = mask_into_bit_vec(superkmer_mask);
    // print to check
    /*
    std::cout<< "printing sk_mask indeces: \n";
    for(size_t i = 0; i < superkmer_idxs.size(); i++){
        std::cout<< superkmer_idxs[i] << " ";
    }
    */
    // elias fano encoding and serialize
    sdsl::sd_vector<> ef_bv (superkmer_idxs.begin(), superkmer_idxs.end()); 
    
    std::cout << "serializing bit vector..." << std::endl;
    bool check = store_to_file(ef_bv, file_sk_mask);
    std::cout<< " successfully stored " << file_sk_mask<<"?: " <<check<<std::endl;


    /*
    uint64_t num_kmers = dict_->size();
    uint64_t one_pm_num_kmers = num_kmers/1000;
    //uint64_t num_super_kmers = dict_->num_superkmers();
    uint64_t dict_m = dict_->m();
    uint64_t dict_seed = dict_->seed();

    std::vector<uint64_t> color_changes_per_superkmer(0);
    std::vector<uint64_t> kmers_per_superkmer(0);

    sshash::dictionary::iterator it = dict_->begin();

    // first kmer
    uint64_t first_kmer_id = 0;
    std::string first_kmer_str = "";
    dict_->access(first_kmer_id, &(first_kmer_str[0]));
    sshash::kmer_t uint_kmer = sshash::util::string_to_uint_kmer(&(first_kmer_str[0]), k_);

    uint64_t count_labels = 0;
    uint64_t count_kmers = 1;

    uint64_t minim, new_minim;
    uint64_t contig_id, new_contig_id;

    // first contig
    contig_id = dict_->lookup_advanced_uint(uint_kmer, false).contig_id;
    new_contig_id = contig_id;
    //first minimizer
    minim = sshash::util::compute_minimizer(uint_kmer, k_, dict_m, dict_seed);
    new_minim = minim;

    //first labels
    std::vector<std::string> labels, new_labels;
    labels = anno_graph->get_labels(first_kmer_str);
    new_labels = labels;
    
    uint64_t kmer_id;
    std::string kmer_str;

    std::cout<< "iterating through graph... \n";
    while(it.has_next()){
        // step to next kmer        
        auto kmer_pair = it.next();
        kmer_str = kmer_pair.second;
        kmer_id = kmer_pair.first;
        uint_kmer = sshash::util::string_to_uint_kmer(&(kmer_str[0]), k_);
        new_minim = sshash::util::compute_minimizer(uint_kmer, k_, dict_m, dict_seed); // is this the correct minimizer?
        new_labels = anno_graph->get_labels(kmer_str);
        new_contig_id = dict_->lookup_advanced_uint(uint_kmer, false).contig_id;


        //next superkmer?
        if(new_minim != minim || new_contig_id != contig_id){//yes
            minim = new_minim;
            contig_id = new_contig_id;
            labels = new_labels;
            color_changes_per_superkmer.push_back(count_labels);
            kmers_per_superkmer.push_back(count_kmers);
            //reset counters
            count_labels = 0;
            count_kmers = 1;
        }else {//no
            if(!equal(new_labels, labels)){
                count_labels++;
                labels = new_labels;
            }
            
            count_kmers++;
        }
        if(kmer_id % one_pm_num_kmers == 0) std::cout<<'.'<<std::flush;
    }
    color_changes_per_superkmer.push_back(count_labels);
    kmers_per_superkmer.push_back(count_kmers);
    std::cout<< "done!\n";
    // sanity checks:
    // 1. per superkmer vectors have same size == num_super_kmers
    sanity_check_1(color_changes_per_superkmer, kmers_per_superkmer, 0);//num_super_kmers);
    // 2. sum of kmers_per_superkmer is num_kmers
    sanity_check_2(kmers_per_superkmer, num_kmers);


    // save stats
    write_stats_to_file(kmers_per_superkmer, color_changes_per_superkmer);
    */
}

bool DBGSSHash::equal(const std::vector<std::string>& input1, const std::vector<std::string>& input2)const {
    if(input1.size() != input2.size()){
        return false;
    }
    for(size_t i = 0; i < input1.size(); i++){
        if(input1.at(i) != input2.at(i)){
            return false;
        }
    }
    return true; 
}

void DBGSSHash::sanity_check_1(const std::vector<uint64_t>& cc_skmer, const std::vector<uint64_t>& km_skmer, uint64_t num_super_kmers)const {
    std::cout << "length of color changes vector: "<<cc_skmer.size()<< std::endl;
    std::cout << "length of kmers vector: "<<km_skmer.size()<< std::endl;
    std::cout << "number of superkmers: "<<num_super_kmers<< std::endl;
}
void DBGSSHash::sanity_check_2(const std::vector<uint64_t>& km_skmer, uint64_t num_kmers)const{
    uint64_t sum = std::accumulate(km_skmer.begin(), km_skmer.end(),0);
    std::cout << "sum of kmers in superkmers vector: "<<sum<< std::endl;
    std::cout << "total number of kmers: "<<num_kmers<< std::endl;
}

} // namespace graph
} // namespace mtg
