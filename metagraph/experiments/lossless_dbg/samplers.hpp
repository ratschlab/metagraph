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

class SamplerConvenient {
public:
    virtual string sample(int length) = 0;
    virtual int reference_size() = 0;
    virtual vector<string> sample_coverage(int length, double coverage) {
        int count = ceil(reference_size()*coverage/length);
        return sample(length, count);
    }
    virtual vector<string> sample(int length, int count) {
        auto res = vector<string>();
        for(int i=0;i<count;i++) {
            res.push_back(sample(length));
        }
        return res;
    }
};

class Sampler : public SamplerConvenient {
public:
    Sampler(string reference, unsigned int seed) : reference(std::move(reference)) {
        generator = std::mt19937(seed); //Standard mersenne_twister_engine seeded with rd()
    };
    string sample(int length) override {
        std::uniform_int_distribution<> dis(0, reference.length()-1-length);
        return reference.substr(dis(generator),length);
    }
    int reference_size() override {
        return reference.size();
    }
private:
    string reference;
    std::mt19937 generator;
};

class DeterministicSampler : public SamplerConvenient {
public:
    DeterministicSampler(vector<string> samples, int reference_size) : _reference_size(reference_size), samples(std::move(samples)) {};
    string sample(int length) override {
        string sample = samples[current_sample];
        assert(length==sample.length());
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
