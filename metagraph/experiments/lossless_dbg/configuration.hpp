//
// Created by Jan Studený on 2019-07-17.
//

#ifndef __CONFIGURATION_HPP__
#define __CONFIGURATION_HPP__

#include <sdsl/rrr_vector.hpp>
#include <sdsl/wt_rlmn.hpp>
// TODO: move to place where
//#define MASK_DUMMY_KMERS
using DefaultUnderlyingWaveletBV = sdsl::rrr_vector<>;
// the default
using WaveletTreeRLMN = sdsl::wt_rlmn<DefaultUnderlyingWaveletBV>; // TODO: inline it
#endif
