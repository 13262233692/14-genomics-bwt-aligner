#pragma once

#include "bwt.h"
#include "suffix_array.h"
#include <string>
#include <vector>
#include <array>
#include <cstdint>

namespace bwt_aligner {

class FMIndex {
public:
    static constexpr int ALPHABET_SIZE = 6;
    static constexpr int OCC_INTERVAL = 512;
    static constexpr int SA_SAMPLE_RATE = 32;
    static constexpr int BWT_PACKED_DENSITY = 4;

    FMIndex();
    ~FMIndex();

    void build(const std::string& text);
    void build(const std::string& text, const SuffixArray& sa, const BWT& bwt);

    int64_t get_c(char c) const;
    int64_t get_occ(int64_t row, char c) const;

    const SuffixArray& get_sa() const { return sa_; }
    const BWT& get_bwt() const { return bwt_; }
    const std::string& get_text() const { return text_; }
    int64_t size() const { return size_; }

    int64_t get_sa_value(int64_t sa_index) const;
    bool is_sa_sample(int64_t index) const { return index % SA_SAMPLE_RATE == 0; }
    int64_t resolve_position(int64_t row) const;

    void clear();

    size_t memory_usage() const;

private:
    std::string text_;
    SuffixArray sa_;
    BWT bwt_;
    int64_t size_;

    std::array<int64_t, ALPHABET_SIZE> c_table_;
    std::vector<uint32_t> occ_packed_;
    int64_t occ_bucket_count_;
    std::vector<int64_t> sa_samples_;

    void build_c_table();
    void build_occ_table();
    void build_sa_samples();
    int64_t get_occ_from_packed(int64_t bucket_idx, int64_t char_idx) const;
    int64_t count_bwt_in_range(int64_t start, int64_t end, char c) const;
    int char_to_idx(char c) const;
};

}
