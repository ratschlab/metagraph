#ifndef __SERIALIZATION_HPP__
#define __SERIALIZATION_HPP__

#include <string>
#include <set>
#include <unordered_map>
#include <map>
#include <cmath>
#include <assert.h>

#include <libmaus2/util/NumberSerialisation.hpp>
#include <libmaus2/util/StringSerialisation.hpp>

void serialise_annotation(std::ostream & out, std::deque<uint32_t> & D) {
    libmaus2::util::NumberSerialisation::serialiseNumber(out, D.size()); 
    std::map<uint32_t, uint32_t> MD;
    std::map<uint32_t, uint32_t> DM;
    for (std::deque<std::uint32_t>::iterator it = D.begin(); it != D.end(); it++) {
        if (DM.find(*it) == DM.end()) {
            MD.insert(std::make_pair(DM.size(), *it));
            DM.insert(std::make_pair(*it, DM.size()));
        }
    }
    assert(DM.size() == MD.size());

    uint32_t bits = (uint32_t) std::ceil(std::log2((double) MD.size()));
    size_t chunks = 64 / bits;

    libmaus2::util::NumberSerialisation::serialiseNumber(out, bits); 

    uint64_t value = 0;
    size_t chunk = 0;
    for (std::deque<std::uint32_t>::iterator it = D.begin(); it != D.end(); it++) {
        if (chunk < chunks) {
            for (size_t i = 0; i < bits; i++)
                value |= ((1 << (chunk * bits + i)) & *it);
            chunk++;
        } else {
            libmaus2::util::NumberSerialisation::serialiseNumber(out, value);         
            value = 0;
            chunk = 0;
        }
    }
    if (chunk > 0)
        libmaus2::util::NumberSerialisation::serialiseNumber(out, value);         

    serilize_annotation_map(out, DM);
    serilize_annotation_map(out, MD);
}

void serialize_combination_vector(std::ostream & out, std::vector<std::uint32_t> const & V) {
    libmaus2::util::NumberSerialisation::serialiseNumber32Vector(out, V);
}

void serialize_annotation_map(std::ostream & out, std::map<uint32_t, uint32_t> const & M) {
    libmaus2::util::NumberSerialisation::serialiseNumber(out, M.size()); 
    for (std::map<uint32_t, uint32_t>::const_iterator it = M.begin(); it != M.end(); ++it) {
        libmaus2::util::NumberSerialisation::serialiseNumber32(out, it->first); 
        libmaus2::util::NumberSerialisation::serialiseNumber32(out, it->second); 
    }
}

/*void serialize_annotation_map(std::ostream & out, std::unordered_map<std::uint32_t, std::set<uint32_t> > const & M ) {
    libmaus2::util::NumberSerialisation::serialiseNumber(out, M.size()); 

    for (std::unordered_map<uint32_t, std::set<uint32_t> >::const_iterator it = M.begin(); it != M.end(); ++it) {
        libmaus2::util::NumberSerialisation::serialiseNumber32(out, it->first); 
        libmaus2::util::NumberSerialisation::serialiseNumber32Set(out, it->second);
        //for (std::set<std::string>::iterator itt = it->second.begin(); itt != it->second.end(); ++itt)
        //    std::cerr << *itt << ": " << sizeof(*itt) << " size: " << itt->size() << " capacity: " << itt->capacity() << std::endl;
        //std::cerr << "set size " << it->second.size() << std::endl;
    }
}*/

void serialize_label_to_id_map(std::ostream & out, std::unordered_map<std::string, uint32_t> const & M) {
    libmaus2::util::NumberSerialisation::serialiseNumber(out, M.size()); 
    for (std::unordered_map<std::string, uint32_t>::const_iterator it = M.begin(); it != M.end(); ++it) {
        libmaus2::util::StringSerialisation::serialiseString(out, it->first); 
        libmaus2::util::NumberSerialisation::serialiseNumber32(out, it->second);
    }
}

void serialize_annotation_id_vector(std::ostream & out, std::vector<std::string> const & V) {
    libmaus2::util::StringSerialisation::serialiseStringVector(out, V);
}

//void deserialize_annotation_map(std::istream & in, std::unordered_map<std::uint32_t, std::set<uint32_t> > & M) {
std::vector<uint32_t> deserialize_combination_vector(std::istream & in) {
    return libmaus2::util::NumberSerialisation::deserialiseNumber32Vector<uint32_t>(in);
}

void deserialize_annotation_map(std::istream & in, std::map<std::uint32_t, uint32_t> & M) {   
    size_t const num = libmaus2::util::NumberSerialisation::deserialiseNumber(in); 
    for (size_t i = 0; i < num; ++i) {
        uint32_t key = libmaus2::util::NumberSerialisation::deserialiseNumber32(in); 
        uint32_t value = libmaus2::util::NumberSerialisation::deserialiseNumber32(in); 
        M.insert(std::make_pair(key, value));
    }
};

/*    size_t const num = libmaus2::util::NumberSerialisation::deserialiseNumber(in); 
    uint64_t total_size = 0;
    std::cerr << "map has " << num << " elements" << std::endl;
    for (size_t i = 0; i < num; ++i) {
        if (i % 1000 == 0) {
            std::cerr << i << "/" << num << std::endl;
            std::cerr << "total size read: " << total_size << " byte / " << total_size / 1024 << " Kbyte / " << total_size / 1024 / 1024 << " MByte" << std::endl;
        }
        uint32_t key = libmaus2::util::NumberSerialisation::deserialiseNumber32(in); 
        std::set<uint32_t> value = libmaus2::util::NumberSerialisation::deserialiseNumber32Set<uint32_t>(in);
        //for (std::set<std::string>::iterator it = value.begin(); it != value.end(); ++ it)
        //    total_size += it->capacity();
        //total_size += sizeof(key);
        M.insert(std::make_pair(key, value));
        //M[key] = value;
    }
}*/

void deserialize_label_to_id_map(std::istream & in, std::unordered_map<std::string, uint32_t> & M) {
    size_t const num = libmaus2::util::NumberSerialisation::deserialiseNumber(in); 
    for (size_t i = 0; i < num; ++i) {
        std::string key = libmaus2::util::StringSerialisation::deserialiseString(in); 
        uint32_t value = libmaus2::util::NumberSerialisation::deserialiseNumber32(in); 
        M.insert(std::make_pair(key, value));
    }
}

std::vector<std::string> deserialize_annotation_id_vector(std::istream & in) {
    return libmaus2::util::StringSerialisation::deserialiseStringVector(in);
}


#endif
