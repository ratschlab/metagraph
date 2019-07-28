//
// Created by Jan Studený on 2019-07-17.
//

#ifndef METAGRAPH_CONFIG_HPP
#define METAGRAPH_CONFIG_HPP

#include <sdsl/rrr_vector.hpp>
#include <sdsl/wt_rlmn.hpp>

#define MASK_DUMMY_KMERS
using DefaultUnderlyingWaveletBV = sdsl::rrr_vector<>;
using DefaultWavelet = sdsl::wt_rlmn<DefaultUnderlyingWaveletBV,DefaultUnderlyingWaveletBV::rank_1_type,DefaultUnderlyingWaveletBV::select_1_type,sdsl::wt_huff<>>;
#endif //METAGRAPH_CONFIG_HPP
