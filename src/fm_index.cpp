#include "fm_index.h"
#include <iostream>
#include <cstring>

namespace bwt_aligner {

FMIndex::FMIndex()
    : size_(0)
    , occ_bucket_count_(0) {
    c_table_.fill(0);
}

FMIndex::~FMIndex() {
}

int FMIndex::char_to_idx(char c) const {
    switch (c) {
        case '$': return 0;
        case '#': return 1;
        case 'A': case 'a': return 2;
        case 'C': case 'c': return 3;
        case 'G': case 'g': return 4;
        case 'T': case 't': return 5;
        default: return 1;
    }
}

void FMIndex::build_c_table() {
    std::array<int64_t, ALPHABET_SIZE> counts;
    counts.fill(0);

    const std::string& bwt = bwt_.get_bwt();
    for (char c : bwt) {
        int idx = char_to_idx(c);
        if (idx >= 0 && idx < ALPHABET_SIZE) {
            ++counts[idx];
        }
    }

    c_table_[0] = 0;
    int64_t cumulative = counts[0];
    for (int i = 1; i < ALPHABET_SIZE; ++i) {
        c_table_[i] = cumulative;
        cumulative += counts[i];
    }
}

void FMIndex::build_occ_table() {
    const std::string& bwt = bwt_.get_bwt();
    int64_t n = static_cast<int64_t>(bwt.size());

    occ_bucket_count_ = (n + OCC_INTERVAL - 1) / OCC_INTERVAL + 1;
    occ_packed_.clear();
    occ_packed_.resize(static_cast<size_t>(occ_bucket_count_ * ALPHABET_SIZE), 0);

    std::array<int64_t, ALPHABET_SIZE> running;
    running.fill(0);

    for (int64_t i = 0; i < n; ++i) {
        if (i % OCC_INTERVAL == 0) {
            int64_t bucket = i / OCC_INTERVAL;
            for (int c = 0; c < ALPHABET_SIZE; ++c) {
                uint32_t val = static_cast<uint32_t>(running[c]);
                occ_packed_[static_cast<size_t>(bucket * ALPHABET_SIZE + c)] = val;
            }
        }
        int idx = char_to_idx(bwt[static_cast<size_t>(i)]);
        if (idx >= 0 && idx < ALPHABET_SIZE) {
            ++running[idx];
        }
    }

    {
        int64_t bucket = occ_bucket_count_ - 1;
        for (int c = 0; c < ALPHABET_SIZE; ++c) {
            uint32_t val = static_cast<uint32_t>(running[c]);
            occ_packed_[static_cast<size_t>(bucket * ALPHABET_SIZE + c)] = val;
        }
    }
}

int64_t FMIndex::get_occ_from_packed(int64_t bucket_idx, int64_t char_idx) const {
    if (bucket_idx < 0) return 0;
    if (bucket_idx >= occ_bucket_count_) bucket_idx = occ_bucket_count_ - 1;
    size_t pos = static_cast<size_t>(bucket_idx * ALPHABET_SIZE + char_idx);
    if (pos >= occ_packed_.size()) return 0;
    return static_cast<int64_t>(occ_packed_[pos]);
}

void FMIndex::build_sa_samples() {
    int64_t n = static_cast<int64_t>(sa_.size());
    sa_samples_.clear();
    sa_samples_.reserve(static_cast<size_t>((n + SA_SAMPLE_RATE - 1) / SA_SAMPLE_RATE));

    for (int64_t i = 0; i < n; i += SA_SAMPLE_RATE) {
        sa_samples_.push_back(sa_[static_cast<size_t>(i)]);
    }
}

int64_t FMIndex::count_bwt_in_range(int64_t start, int64_t end, char c) const {
    int c_idx = char_to_idx(c);
    if (c_idx < 0 || c_idx >= ALPHABET_SIZE) return 0;

    int64_t count = 0;
    const std::string& bwt = bwt_.get_bwt();
    for (int64_t i = start; i < end && i < static_cast<int64_t>(bwt.size()); ++i) {
        if (char_to_idx(bwt[static_cast<size_t>(i)]) == c_idx) {
            ++count;
        }
    }
    return count;
}

void FMIndex::build(const std::string& text, const SuffixArray& sa, const BWT& bwt) {
    clear();

    text_ = text;
    size_ = static_cast<int64_t>(text_.size());
    sa_ = sa;
    bwt_ = bwt;

    build_c_table();
    build_occ_table();
    build_sa_samples();
}

void FMIndex::build(const std::string& text) {
    clear();

    SuffixArray sa;
    BWT bwt;

    std::string text_with_dollar = text + '$';
    sa.build(text_with_dollar);
    bwt.build(text_with_dollar, sa);

    build(text_with_dollar, sa, bwt);
}

int64_t FMIndex::get_c(char c) const {
    int idx = char_to_idx(c);
    if (idx < 0 || idx >= ALPHABET_SIZE) return 0;
    return c_table_[static_cast<size_t>(idx)];
}

int64_t FMIndex::get_occ(int64_t row, char c) const {
    if (row < 0) return 0;
    if (row >= size_) row = size_ - 1;

    int c_idx = char_to_idx(c);
    if (c_idx < 0 || c_idx >= ALPHABET_SIZE) return 0;

    int64_t target = row + 1;
    int64_t bucket = target / OCC_INTERVAL;
    int64_t rem = target % OCC_INTERVAL;

    int64_t count = get_occ_from_packed(bucket, c_idx);

    int64_t start = bucket * OCC_INTERVAL;
    count += count_bwt_in_range(start, start + rem, c);

    return count;
}

int64_t FMIndex::get_sa_value(int64_t sa_index) const {
    return sa_[static_cast<size_t>(sa_index)];
}

int64_t FMIndex::resolve_position(int64_t row) const {
    int64_t steps = 0;
    const std::string& bwt_str = bwt_.get_bwt();

    while (row >= 0 && row < size_ && !is_sa_sample(row) && steps < size_) {
        char c = bwt_str[static_cast<size_t>(row)];
        int64_t occ = (row > 0) ? get_occ(row - 1, c) : 0;
        row = get_c(c) + occ;
        ++steps;
    }

    if (row >= 0 && row < size_ && is_sa_sample(row)) {
        int64_t sample_idx = row / SA_SAMPLE_RATE;
        if (sample_idx >= 0 && sample_idx < static_cast<int64_t>(sa_samples_.size())) {
            int64_t sa_val = sa_samples_[static_cast<size_t>(sample_idx)];
            int64_t result = (sa_val + steps) % size_;
            if (result < 0) result += size_;
            return result;
        }
    }

    return -1;
}

size_t FMIndex::memory_usage() const {
    size_t total = 0;
    total += text_.capacity();
    total += sa_.size() * sizeof(int64_t);
    total += bwt_.size();
    total += occ_packed_.capacity() * sizeof(uint32_t);
    total += sa_samples_.capacity() * sizeof(int64_t);
    total += sizeof(c_table_);
    return total;
}

void FMIndex::clear() {
    text_.clear();
    text_.shrink_to_fit();
    sa_.clear();
    bwt_.clear();
    size_ = 0;
    c_table_.fill(0);
    occ_packed_.clear();
    occ_packed_.shrink_to_fit();
    occ_bucket_count_ = 0;
    sa_samples_.clear();
    sa_samples_.shrink_to_fit();
}

}
