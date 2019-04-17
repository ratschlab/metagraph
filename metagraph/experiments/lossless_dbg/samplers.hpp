//
//  samplers.hpp
//  lossless_dbg
//
//  Created by Jan Studený on 11/03/2019.
//  Copyright © 2019 Jan Studený. All rights reserved.
//

#ifndef samplers_h
#define samplers_h
#include <utility>
#include <iostream>
#include <map>
#include <vector>
#include <random>
#include <cassert>

using namespace std;

class SamplerConvenient {
public:
    virtual string sample(int sample_length) = 0;
    virtual int reference_size() = 0;
    virtual vector<string> sample_coverage(int sample_length, double coverage) {
        int count = ceil(reference_size()*coverage/sample_length);
        return sample(sample_length, count);
    }
    virtual vector<string> sample(int sample_length, int count) {
        auto res = vector<string>();
        for(int i=0;i<count;i++) {
            res.push_back(sample(sample_length));
        }
        return res;
    }
};

class NoisySampler : public SamplerConvenient {
public:
    NoisySampler(string reference, std::mt19937 generator,double probability_of_error=0) : reference(std::move(reference)), generator(generator) {
    };
    string sample(int sample_length) override {
        std::uniform_int_distribution<> dis(0, reference.length()-1-sample_length);
        string output = reference.substr(dis(generator),sample_length);
        for(auto& e : output) {
            map<char,int> to_offset = {{'A',0,},{'C',1},{'G',2},{'T',3}};
            char to_char[] = {'A','C','G','T'};
            std::discrete_distribution<> err({1-probability_of_error,
                                              probability_of_error/3,
                                              probability_of_error/3,
                                              probability_of_error/3});
            e = to_char[(to_offset[e]+err(generator))%4];
        }
    }
    int reference_size() override {
        return reference.size();
    }
private:
    double probability_of_error;
    string reference;
    std::mt19937 generator;
};
class SubSampler : public NoisySampler {
public:
    SubSampler(const string& reference, int subsample_size, std::mt19937 generator,double probability_of_error) :
        NoisySampler(reference.substr(
                std::uniform_int_distribution<>(0,reference.length()-min(subsample_size,(int)reference.length()))(generator),
                min(subsample_size,(int)reference.length())
                ),
                generator,
                probability_of_error)
                {}
};

class DeterministicSampler : public SamplerConvenient {
public:
    DeterministicSampler(vector<string> samples, int reference_size) : _reference_size(reference_size), samples(std::move(samples)) {};
    string sample(int sample_length) override {
        string sample = samples[current_sample];
        assert(sample_length==sample.length());
        current_sample = (current_sample + 1) % samples.size();
        return sample;
    }
    int reference_size() override {
        return _reference_size;
    }
    vector<string> samples;
    int _reference_size;
    int current_sample = 0;
};

#endif /* samplers_h */
