#pragma once

#include "fm_index.h"
#include "sequence.h"
#include <string>
#include <vector>
#include <cstdint>

namespace bwt_aligner {

struct SearchRange {
    int64_t low;
    int64_t high;
    bool valid() const { return low <= high; }
    int64_t count() const { return high - low + 1; }
};

struct SeedHit {
    int64_t ref_position;
    int32_t query_start;
    int32_t query_end;
    int32_t num_mismatches;
    bool is_reverse;
};

class BackwardSearch {
public:
    BackwardSearch();
    ~BackwardSearch();

    void set_fm_index(const FMIndex* fm_index);

    SearchRange search(const std::string& pattern) const;
    SearchRange search_reverse(const std::string& pattern) const;

    std::vector<SeedHit> find_seeds(const std::string& query,
                                    int32_t seed_length,
                                    int32_t max_mismatches = 0) const;

    std::vector<SeedHit> find_seeds_both_strands(const std::string& query,
                                                 int32_t seed_length,
                                                 int32_t max_mismatches = 0) const;

    bool verify_alignment(const std::string& query,
                          int64_t ref_position,
                          bool is_reverse,
                          int32_t max_mismatches,
                          int32_t& num_mismatches) const;

private:
    const FMIndex* fm_index_;

    SearchRange extend_left(SearchRange range, char c) const;
    int64_t lf_mapping(int64_t row, char c) const;
};

}
