#pragma once

#include <string>
#include <vector>
#include <cstdint>

namespace bwt_aligner {

struct ReferenceSequence {
    std::string name;
    std::string sequence;
    uint64_t length;
    uint64_t global_offset;
};

struct FastqRead {
    std::string id;
    std::string sequence;
    std::string quality;
    bool is_first;
    int32_t pair_id;
};

struct AlignmentResult {
    std::string read_id;
    std::string ref_name;
    int32_t ref_id;
    int64_t position;
    bool is_reverse;
    bool is_first;
    int32_t mapq;
    std::string cigar;
    std::string sequence;
    std::string quality;
    bool is_secondary;
    bool is_supplementary;
    bool is_paired;
    bool is_proper_pair;
    bool mate_is_reverse;
    int32_t mate_ref_id;
    int64_t mate_position;
    int64_t insert_size;
    bool is_unmapped;
    bool mate_is_unmapped;
    int32_t num_mismatches;
};

inline char complement_base(char c) {
    switch (c) {
        case 'A': case 'a': return 'T';
        case 'T': case 't': return 'A';
        case 'C': case 'c': return 'G';
        case 'G': case 'g': return 'C';
        case 'N': case 'n': return 'N';
        default: return 'N';
    }
}

inline std::string reverse_complement(const std::string& seq) {
    std::string result;
    result.reserve(seq.size());
    for (auto it = seq.rbegin(); it != seq.rend(); ++it) {
        result.push_back(complement_base(*it));
    }
    return result;
}

inline uint8_t base_to_index(char c) {
    switch (c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default: return 4;
    }
}

inline char index_to_base(uint8_t idx) {
    switch (idx) {
        case 0: return 'A';
        case 1: return 'C';
        case 2: return 'G';
        case 3: return 'T';
        default: return 'N';
    }
}

}
