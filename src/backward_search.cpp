#include "backward_search.h"
#include <iostream>
#include <algorithm>

namespace bwt_aligner {

BackwardSearch::BackwardSearch()
    : fm_index_(nullptr) {
}

BackwardSearch::~BackwardSearch() {
}

void BackwardSearch::set_fm_index(const FMIndex* fm_index) {
    fm_index_ = fm_index;
}

SearchRange BackwardSearch::extend_left(SearchRange range, char c) const {
    if (!fm_index_ || !range.valid()) {
        return {0, -1};
    }

    int64_t c_count = fm_index_->get_c(c);
    int64_t occ_low = (range.low > 0) ? fm_index_->get_occ(range.low - 1, c) : 0;
    int64_t occ_high = fm_index_->get_occ(range.high, c);

    return {c_count + occ_low, c_count + occ_high - 1};
}

int64_t BackwardSearch::lf_mapping(int64_t row, char c) const {
    if (!fm_index_) return -1;
    return fm_index_->get_c(c) + fm_index_->get_occ(row - 1, c);
}

SearchRange BackwardSearch::search(const std::string& pattern) const {
    if (!fm_index_ || pattern.empty()) {
        return {0, -1};
    }

    SearchRange range{0, fm_index_->size() - 1};

    for (auto it = pattern.rbegin(); it != pattern.rend() && range.valid(); ++it) {
        char c = *it;
        if (c == 'N' || c == 'n') {
            return {0, -1};
        }
        range = extend_left(range, c);
    }

    return range;
}

SearchRange BackwardSearch::search_reverse(const std::string& pattern) const {
    std::string rc = reverse_complement(pattern);
    return search(rc);
}

std::vector<SeedHit> BackwardSearch::find_seeds(const std::string& query,
                                                int32_t seed_length,
                                                int32_t max_mismatches) const {
    std::vector<SeedHit> hits;
    if (!fm_index_ || static_cast<int32_t>(query.size()) < seed_length) {
        return hits;
    }

    int32_t n = static_cast<int32_t>(query.size());
    for (int32_t start = 0; start + seed_length <= n; ++start) {
        std::string seed = query.substr(static_cast<size_t>(start), static_cast<size_t>(seed_length));
        SearchRange range = search(seed);

        if (range.valid()) {
            for (int64_t i = range.low; i <= range.high; ++i) {
                int64_t sa_val = fm_index_->get_sa_value(i);
                int64_t ref_pos;
                if (fm_index_->is_sa_sample(i)) {
                    ref_pos = sa_val;
                } else {
                    ref_pos = fm_index_->resolve_position(i);
                }
                if (ref_pos >= 0) {
                    hits.push_back({
                        ref_pos,
                        start,
                        start + seed_length,
                        0,
                        false
                    });
                }
            }
        }
    }

    return hits;
}

std::vector<SeedHit> BackwardSearch::find_seeds_both_strands(const std::string& query,
                                                             int32_t seed_length,
                                                             int32_t max_mismatches) const {
    std::vector<SeedHit> fwd_hits = find_seeds(query, seed_length, max_mismatches);
    std::vector<SeedHit> rev_hits = find_seeds(reverse_complement(query), seed_length, max_mismatches);

    for (auto& hit : rev_hits) {
        hit.is_reverse = true;
        int32_t qs = hit.query_start;
        int32_t qe = hit.query_end;
        int32_t n = static_cast<int32_t>(query.size());
        hit.query_start = n - qe;
        hit.query_end = n - qs;
    }

    fwd_hits.insert(fwd_hits.end(), rev_hits.begin(), rev_hits.end());
    return fwd_hits;
}

bool BackwardSearch::verify_alignment(const std::string& query,
                                      int64_t ref_position,
                                      bool is_reverse,
                                      int32_t max_mismatches,
                                      int32_t& num_mismatches) const {
    if (!fm_index_) return false;

    const std::string& text = fm_index_->get_text();
    int32_t query_len = static_cast<int32_t>(query.size());

    if (ref_position < 0 || ref_position + query_len > static_cast<int64_t>(text.size())) {
        return false;
    }

    std::string q = is_reverse ? reverse_complement(query) : query;
    num_mismatches = 0;

    for (int32_t i = 0; i < query_len; ++i) {
        char qc = q[static_cast<size_t>(i)];
        char rc = text[static_cast<size_t>(ref_position + i)];

        if (qc == 'N' || qc == 'n' || rc == 'N' || rc == 'n') {
            continue;
        }
        if (qc != rc && qc != std::toupper(rc) && std::toupper(qc) != rc) {
            ++num_mismatches;
            if (num_mismatches > max_mismatches) {
                return false;
            }
        }
    }

    return true;
}

}
