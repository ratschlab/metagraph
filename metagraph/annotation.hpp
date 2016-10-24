#ifndef __ANNOTATION_HPP__
#define __ANNOTATION_HPP__

#include <vector>
#include <algorithm>

#include "datatypes.hpp"

std::vector<uint32_t> get_curr_combination(std::vector<uint32_t> &combinations, uint32_t idx) {
    std::vector<uint32_t>::iterator start = combinations.begin() + idx;
    uint32_t length = *start;
    std::vector<uint32_t>::iterator end = combinations.begin() + idx + length + 1;
    start++;

    std::vector<uint32_t> curr_combo(length);
    std::copy(start, end, curr_combo.begin());

    return curr_combo;
};

/*bool in_combo(std::vector<uint32_t> &combo, uint32_t value) {
    
    std::cerr << "search for " << value << " at idx " << idx << std::endl;
    std::vector<uint32_t>::iterator start = anno.begin() + idx;
    std::cerr << "start val " << *start << std::endl;
    std::vector<uint32_t>::iterator end = anno.begin() + idx + *start;
    start++;

    return std::binary_search(start, end, value);
};*/

std::vector<uint32_t> add_to_combination(std::vector<uint32_t> &old_combination, uint32_t value) {

    std::vector<uint32_t> new_combination;

    if (old_combination.size() > 0) {
        std::vector<uint32_t>::iterator insert_idx = std::upper_bound(old_combination.begin(), old_combination.end(), value);
        std::vector<uint32_t>::iterator it;
        for (it = old_combination.begin(); it != insert_idx; ++it)
            new_combination.push_back(*it);
        new_combination.push_back(value);
        for (it = insert_idx; it != old_combination.end(); ++it)
            new_combination.push_back(*it);
    } else {
        new_combination.push_back(value);
    }

    return new_combination;
};

uint32_t insert_new_combination(std::vector<uint32_t> &combination_vector, std::vector<uint32_t> &new_combination) {

    uint32_t insert_pos = combination_vector.size();
    //std::cerr << "ins pos " << insert_pos << std::endl;

    // are we close to the next billion? -> pre-allocate another billion entries
    // this prevents the vector size from doubling with each re-allocation
    // we give up on the constant amortized cost, but this is ok for here

    if ((insert_pos + 1000) % 1000000000 == 0)
        combination_vector.reshape(insert_pos + 1000 + 1000000000);

    combination_vector.push_back(new_combination.size());
    for (std::vector<uint32_t>::iterator it = new_combination.begin(); it != new_combination.end(); ++it)
        combination_vector.push_back(*it);
    return insert_pos;
};
#endif
