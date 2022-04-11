#include "row_sparse_disk.hpp"

#include "common/logger.hpp"
#include "common/utils/file_utils.hpp"
#include "common/ifstream_with_name_and_offset.hpp"
#include "graph/representation/succinct/boss.hpp"


namespace mtg {
namespace annot {
namespace binmat {

using mtg::common::logger;

RowSparseDisk::RowSparseDisk(const std::function<void(const RowCallback &)> &/*call_rows*/,
                     uint64_t num_columns,
                     uint64_t num_rows,
                     uint64_t /*num_relations*/)
      : num_columns_(num_columns),
        num_rows_(num_columns > 0 ? num_rows : 0) {
    
    throw std::runtime_error("Currently not supported (mkokot)");
}

bool RowSparseDisk::get(Row row, Column column) const {
    SetBitPositions set_bits = get_row(row);
    SetBitPositions::iterator v = std::lower_bound(set_bits.begin(), set_bits.end(), column);
    return v != set_bits.end() && *v == column;
}

std::vector<BinaryMatrix::Row> RowSparseDisk::get_column(Column column) const {
    std::vector<Row> result;
    for (Row row = 0; row < num_rows(); ++row) {
        if (get(row, column))
            result.push_back(row);
    }
    return result;
}

BinaryMatrix::SetBitPositions RowSparseDisk::get_row(Row row) const {        
    assert(boundary_[boundary_.size() - 1] == 1);
    uint64_t start_idx = row == 0 ? 0 : boundary_.select1(row) + 1;
    uint64_t end_idx = boundary_.next1(start_idx);
    SetBitPositions result;
    result.reserve(end_idx - start_idx);
    // In each row, the first value in `set_bits_` stores the first set bit,
    // and all next ones store deltas pos[i] - pos[i-1].
    uint64_t last = 0;
    for (uint64_t i = start_idx; i != end_idx; ++i) {
        result.push_back(last += set_bits_[i - row]); 
    }
    return result;
}

bool RowSparseDisk::load(std::istream &f) {    
    auto _f = dynamic_cast<IfstreamWithNameAndOffset*>(&f);
    if (!_f) {        
        logger->error("RowSparseDisk::load requires IfstreamWithNameAndOffset object"); // mkokot fix this message
        return false;
    }   
    try {        
        f.read(reinterpret_cast<char *>(&num_columns_), sizeof(uint64_t));
        
        _f->SetOffset(static_cast<uint64_t>(f.tellg()));

        set_bits_ = sdsl::int_vector_buffer<>(_f->GetFName(), std::ios::in, 1024*1024, 0, false, _f->GetOffset());
        
        f.seekg(-sizeof(std::streampos), ios_base::end);
        std::streampos boundary_start;
        f.read(reinterpret_cast<char*>(&boundary_start), sizeof(boundary_start));
        //std::cerr << "boundary_start: " << boundary_start << "\n";
        f.seekg(boundary_start, ios_base::beg);
        
        boundary_.load(f);
        
        num_rows_ = boundary_.num_set_bits();

        // mkokot, TODO: remove        
        ////for(uint64_t i = 0 ; i < num_rows() ; ++i)
        //for(uint64_t i = 0 ; i < 1000 ; ++i)
        //{            
        //    uint64_t start_idx = i == 0 ? 0 : boundary_.select1(i) + 1;
        //    uint64_t end_idx = boundary_.next1(start_idx);
        //    std::cerr << i << " " << start_idx << " " << end_idx << ":";
        //    auto _row = get_row(i);
        //    for(auto x : _row)
        //        std::cerr << " " << x;
        //    std::cerr << "\n";
        //}
        
        
    } catch (...) {
        return false;
    }

    return true;
}

void RowSparseDisk::serialize(std::ostream &f) const  {
    throw std::runtime_error("Currently not supported (mkokot), the static variant should be used instead");
    f.write(reinterpret_cast<const char *>(&num_columns_), sizeof(uint64_t));
    //set_bits_.serialize(f); // mkokot consider this, I'm not sure if it is beeding called, throw runtime_error here
    boundary_.serialize(f);
}


void RowSparseDisk::serialize(const std::function<void(binmat::BinaryMatrix::RowCallback)> &call_rows,
                              const std::string& filename, uint64_t num_cols, uint64_t num_set_bits, uint64_t num_rows) {
    
    std::ofstream outstream(filename, std::ios::binary | std::ios::ate | std::ios::in);    //std::ios_base::in needed ...
    
    if (!outstream.good())
        throw std::ofstream::failure("Cannot write to file " + filename);

    outstream.write(reinterpret_cast<char *>(&num_cols), sizeof(uint64_t));

    const uint64_t iv_offs = outstream.tellp();
    outstream.close();

    auto outbuf = sdsl::int_vector_buffer<>(filename,
                                            std::ios::out | std::ios::binary,
                                            1024 * 1024,
                                            sdsl::bits::hi(num_cols) + 1,
                                            false,
                                            iv_offs);
    // mkokot, just recall: call_bits should be called for each argument that should be set to true (row boundaries)
    auto call_bits = [&](const std::function<void(uint64_t)> &call_bit) {
        uint64_t t = 0;
        call_rows([&](const auto &row) {
            uint64_t last = 0;
            for (auto val : row) {
                //outbuf.push_back(val); // mkokot, will the part below fix the issues?
                outbuf.push_back(val - last);
                last = val;
            }
            call_bit((t += row.size() + 1) - 1);                
        });
        outbuf.close();
    };

    bit_vector_small boundary(call_bits, num_rows + num_set_bits, num_rows);

    outstream.open(filename, std::ios::in | std::ios::out | std::ios::binary | std::ios::ate | std::ios::in); 
    std::streampos boundary_start = outstream.tellp();
    //std::cerr << "boundary_start: " << boundary_start << "\n";
    boundary.serialize(outstream);    
    outstream.write(reinterpret_cast<char*>(&boundary_start), sizeof(std::streampos));
}

} // namespace binmat
} // namespace annot
} // namespace mtg
