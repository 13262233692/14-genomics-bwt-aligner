#include "suffix_array.h"
#include <algorithm>
#include <iostream>

namespace bwt_aligner {

SuffixArray::SuffixArray() : n_(0), k_(0) {
}

SuffixArray::~SuffixArray() {
}

bool SuffixArray::compare_sa(int64_t i, int64_t j) {
    if (rank_[static_cast<size_t>(i)] != rank_[static_cast<size_t>(j)])
        return rank_[static_cast<size_t>(i)] < rank_[static_cast<size_t>(j)];
    int64_t ri = i + k_ < n_ ? rank_[static_cast<size_t>(i + k_)] : -1;
    int64_t rj = j + k_ < n_ ? rank_[static_cast<size_t>(j + k_)] : -1;
    return ri < rj;
}

void SuffixArray::build_sa(const std::string& s) {
    n_ = static_cast<int64_t>(s.size());
    sa_.resize(static_cast<size_t>(n_));
    rank_.resize(static_cast<size_t>(n_));
    tmp_.resize(static_cast<size_t>(n_));

    for (int64_t i = 0; i < n_; ++i) {
        sa_[static_cast<size_t>(i)] = i;
        rank_[static_cast<size_t>(i)] = static_cast<int64_t>(static_cast<unsigned char>(s[static_cast<size_t>(i)]));
    }

    rank_[static_cast<size_t>(n_ - 1)] = -1;

    for (k_ = 1; k_ < n_; k_ *= 2) {
        std::sort(sa_.begin(), sa_.end(), [this](int64_t a, int64_t b) {
            return this->compare_sa(a, b);
        });

        tmp_[static_cast<size_t>(sa_[0])] = 0;
        for (int64_t i = 1; i < n_; ++i) {
            int64_t prev = static_cast<int64_t>(tmp_[static_cast<size_t>(sa_[static_cast<size_t>(i - 1)])]);
            tmp_[static_cast<size_t>(sa_[static_cast<size_t>(i)])] = prev + (compare_sa(sa_[static_cast<size_t>(i - 1)], sa_[static_cast<size_t>(i)]) ? 1 : 0);
        }
        for (int64_t i = 0; i < n_; ++i) {
            rank_[static_cast<size_t>(i)] = tmp_[static_cast<size_t>(i)];
        }
    }
}

void SuffixArray::build(const std::string& text) {
    clear();
    build_sa(text);
}

void SuffixArray::clear() {
    sa_.clear();
    rank_.clear();
    tmp_.clear();
    n_ = 0;
    k_ = 0;
}

}
