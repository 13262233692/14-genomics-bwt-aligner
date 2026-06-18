#pragma once

#include "sequence.h"
#include <string>
#include <vector>
#include <fstream>
#include <cstdint>

namespace bwt_aligner {

class FastaReader {
public:
    FastaReader();
    ~FastaReader();

    bool open(const std::string& filename);
    void close();
    bool read_all(std::vector<ReferenceSequence>& sequences);

    uint64_t total_length() const { return total_length_; }
    size_t num_sequences() const { return num_sequences_; }

private:
    std::ifstream file_;
    uint64_t total_length_;
    size_t num_sequences_;
    bool is_open_;
};

}
