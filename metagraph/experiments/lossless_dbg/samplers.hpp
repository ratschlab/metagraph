//
//  samplers.hpp
//  lossless_dbg
//
//  Created by Jan Studený on 11/03/2019.
//  Copyright © 2019 Jan Studený. All rights reserved.
//

#ifndef __SAMPLERS_HPP__
#define __SAMPLERS_HPP__
#include <utility>
#include <iostream>
#include <map>
#include <vector>
#include <random>
#include <cassert>

using namespace std;

class SamplerConvenient {
public:
    virtual string sample(int64_t sample_length) = 0;
    virtual int64_t reference_size() = 0;
    virtual vector<string> sample_coverage(int64_t sample_length, double coverage) {
        int64_t count = ceil(reference_size()*coverage/sample_length);
        return sample(sample_length, count);
    }
    virtual vector<string> sample(int64_t sample_length, int64_t count) {
        auto res = vector<string>();
        for(int64_t i=0;i<count;i++) {
            res.push_back(sample(sample_length));
        }
        return res;
    }
};

class NoisySampler : public SamplerConvenient {
public:
    NoisySampler(string reference,
                 std::mt19937 generator,
                 double probability_of_error=0) : reference(std::move(reference)),
                                                  generator(generator),
                                                  probability_of_error(probability_of_error),
                                                  err({1-probability_of_error,
                                                       probability_of_error/3,
                                                       probability_of_error/3,
                                                       probability_of_error/3}){
    };
    string sample(int64_t sample_length) override {
        std::uniform_int_distribution<> dis(0, reference.length()-1-sample_length);
        string output = reference.substr(dis(generator),sample_length);
        for(auto& e : output) {
            map<char,int64_t> to_offset = {{'A',0,},{'C',1},{'G',2},{'T',3}};
            char to_char[] = {'A','C','G','T'};
            e = to_char[(to_offset[e] + err(generator)) % 4];
        }
        return output;
    }
    int64_t reference_size() override {
        return reference.size();
    }
private:
    string reference;
    std::mt19937 generator;
    double probability_of_error;
    std::discrete_distribution<> err;
};
class SubSampler : public NoisySampler {
public:
    SubSampler(const string& reference, int64_t subsample_size, std::mt19937 generator,double probability_of_error=0) :
        NoisySampler(reference.substr(
                std::uniform_int_distribution<>(0,reference.length()-min(subsample_size,(int64_t)reference.length()))(generator),
                min(subsample_size,(int64_t)reference.length())
                ),
                generator,
                probability_of_error)
                {}
};

class DeterministicSampler : public SamplerConvenient {
public:
    DeterministicSampler(vector<string> samples, int64_t reference_size) : _reference_size(reference_size), samples(std::move(samples)) {};
    string sample(int64_t sample_length) override {
        string sample = samples[current_sample];
        assert(sample_length==(int64_t )sample.length());
        current_sample = (current_sample + 1) % samples.size();
        return sample;
    }
    int64_t reference_size() override {
        return _reference_size;
    }
    int64_t _reference_size;
    vector<string> samples;
    int64_t current_sample = 0;
};

#endif /* __SAMPLERS_HPP__ */
