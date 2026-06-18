#pragma once

#include <string>
#include <vector>
#include <cstdint>

namespace bwt_aligner {

struct SWConfig {
    int32_t match_score;
    int32_t mismatch_score;
    int32_t gap_open;
    int32_t gap_extend;
    int32_t score_threshold;

    SWConfig()
        : match_score(2)
        , mismatch_score(-2)
        , gap_open(3)
        , gap_extend(1)
        , score_threshold(8)
    {}
};

struct SWResult {
    int32_t score;
    int32_t ref_start;
    int32_t ref_end;
    int32_t query_start;
    int32_t query_end;
    int32_t num_mismatches;
    int32_t num_ins;
    int32_t num_del;
    std::string cigar;
    bool valid;

    SWResult()
        : score(0)
        , ref_start(0)
        , ref_end(0)
        , query_start(0)
        , query_end(0)
        , num_mismatches(0)
        , num_ins(0)
        , num_del(0)
        , valid(false)
    {}
};

class SmithWaterman {
public:
    SmithWaterman();
    ~SmithWaterman();

    void set_config(const SWConfig& config) { config_ = config; }
    const SWConfig& get_config() const { return config_; }

    SWResult align(const std::string& query,
                   const std::string& ref,
                   int32_t ref_start = 0,
                   int32_t ref_end = -1) const;

    SWResult align_avx2(const std::string& query,
                        const std::string& ref,
                        int32_t ref_start = 0,
                        int32_t ref_end = -1) const;

    static std::string build_cigar(const std::vector<int8_t>& traceback,
                                   int32_t query_len,
                                   int32_t ref_len,
                                   int32_t& out_mismatches,
                                   int32_t& out_ins,
                                   int32_t& out_del);

private:
    SWConfig config_;

    SWResult align_basic(const std::string& query,
                         const std::string& ref) const;

    static bool has_avx2();
};

}
