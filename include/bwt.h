#pragma once

#include "suffix_array.h"
#include <string>
#include <vector>
#include <cstdint>

namespace bwt_aligner {

class BWT {
public:
    BWT();
    ~BWT();

    void build(const std::string& text, const SuffixArray& sa);

    const std::string& get_bwt() const { return bwt_str_; }
    char operator[](size_t i) const { return bwt_str_[i]; }
    size_t size() const { return bwt_str_.size(); }

    char get_primary() const { return primary_char_; }
    int64_t get_primary_pos() const { return primary_pos_; }

    void clear();

private:
    std::string bwt_str_;
    char primary_char_;
    int64_t primary_pos_;
};

}
