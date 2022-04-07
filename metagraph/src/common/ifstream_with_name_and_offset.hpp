#pragma once
#include <fstream>

// mkokot should I enclose it in namespace?
class IfstreamWithNameAndOffset : public std::ifstream
{
    std::string fname;
    uint64_t offset{};
public:
    explicit IfstreamWithNameAndOffset(const std::string& fname, std::ios_base::openmode mode) : 
    std::ifstream(fname, mode), fname(fname) { }
    
    const std::string& GetFName() const {
        return fname;
    }
    
    void SetOffset(uint64_t offset) {
        this->offset = offset;
    }
    
    uint64_t GetOffset() const {
        return offset;
    }
};