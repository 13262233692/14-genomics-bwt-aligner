#include "bwt.h"
#include <iostream>

namespace bwt_aligner {

BWT::BWT()
    : primary_char_('$')
    , primary_pos_(-1) {
}

BWT::~BWT() {
}

void BWT::build(const std::string& text, const SuffixArray& sa) {
    clear();

    int64_t n = static_cast<int64_t>(sa.size());
    bwt_str_.resize(static_cast<size_t>(n));

    for (int64_t i = 0; i < n; ++i) {
        int64_t sa_val = sa[static_cast<size_t>(i)];
        if (sa_val == 0) {
            bwt_str_[static_cast<size_t>(i)] = '$';
            primary_pos_ = i;
            primary_char_ = '$';
        } else {
            bwt_str_[static_cast<size_t>(i)] = text[static_cast<size_t>(sa_val - 1)];
        }
    }
}

void BWT::clear() {
    bwt_str_.clear();
    primary_char_ = '$';
    primary_pos_ = -1;
}

}
