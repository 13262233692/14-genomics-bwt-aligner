#include "smith_waterman.h"
#include "sequence.h"
#include <algorithm>
#include <cstring>
#include <sstream>

#if defined(_MSC_VER)
#include <intrin.h>
#elif defined(__GNUC__) || defined(__clang__)
#include <cpuid.h>
#include <immintrin.h>
#endif

namespace bwt_aligner {

SmithWaterman::SmithWaterman() {}
SmithWaterman::~SmithWaterman() {}

bool SmithWaterman::has_avx2() {
#if defined(_MSC_VER)
    int info[4];
    __cpuid(info, 7);
    return (info[1] & 0x20) != 0;
#elif defined(__GNUC__) || defined(__clang__) && defined(__x86_64__)
    unsigned int eax, ebx, ecx, edx;
    __get_cpuid(7, &eax, &ebx, &ecx, &edx);
    return (ebx & 0x20) != 0;
#else
    return false;
#endif
}

enum TraceOp : int8_t {
    TRACE_NONE = 0,
    TRACE_MATCH = 1,
    TRACE_DEL = 2,
    TRACE_INS = 3,
};

std::string SmithWaterman::build_cigar(const std::vector<int8_t>& traceback,
                                       int32_t query_len,
                                       int32_t ref_len,
                                       int32_t& out_mismatches,
                                       int32_t& out_ins,
                                       int32_t& out_del) {
    (void)query_len;
    (void)ref_len;

    if (traceback.empty()) {
        out_mismatches = 0;
        out_ins = 0;
        out_del = 0;
        return "*";
    }

    std::vector<std::pair<char, int32_t>> ops;
    char last_op = ' ';
    int32_t last_count = 0;

    out_ins = 0;
    out_del = 0;

    for (size_t idx = 0; idx < traceback.size(); ++idx) {
        int8_t op = traceback[idx];
        char c = ' ';
        switch (op) {
            case TRACE_MATCH: c = 'M'; break;
            case TRACE_DEL:   c = 'D'; ++out_del; break;
            case TRACE_INS:   c = 'I'; ++out_ins; break;
            default: continue;
        }
        if (c == last_op) {
            ++last_count;
        } else {
            if (last_op != ' ') {
                ops.emplace_back(last_op, last_count);
            }
            last_op = c;
            last_count = 1;
        }
    }
    if (last_op != ' ') {
        ops.emplace_back(last_op, last_count);
    }

    if (ops.empty()) return "*";

    std::ostringstream oss;
    for (const auto& p : ops) {
        oss << p.second << p.first;
    }
    return oss.str();
}

static inline int32_t score_match_mismatch(char a, char b, const SWConfig& cfg) {
    uint8_t ai = base_to_index(a);
    uint8_t bi = base_to_index(b);
    return (ai == bi && ai < 4) ? cfg.match_score : cfg.mismatch_score;
}

SWResult SmithWaterman::align_basic(const std::string& query,
                                    const std::string& ref) const {
    SWResult result;
    int32_t m = static_cast<int32_t>(query.size());
    int32_t n = static_cast<int32_t>(ref.size());

    if (m == 0 || n == 0) return result;

    int32_t gap_open = -config_.gap_open;
    int32_t gap_extend = -config_.gap_extend;

    std::vector<int32_t> H(static_cast<size_t>(n + 1), 0);
    std::vector<int32_t> E(static_cast<size_t>(n + 1), 0);
    std::vector<int32_t> H_prev(static_cast<size_t>(n + 1), 0);
    std::vector<int32_t> F(static_cast<size_t>(n + 1), 0);
    std::vector<int32_t> F_prev(static_cast<size_t>(n + 1), 0);

    std::vector<std::vector<int8_t>> trace(
        static_cast<size_t>(m + 1),
        std::vector<int8_t>(static_cast<size_t>(n + 1), TRACE_NONE));

    int32_t best_score = 0;
    int32_t best_i = 0;
    int32_t best_j = 0;

    for (int32_t i = 1; i <= m; ++i) {
        H[0] = 0;
        E[0] = 0;
        F[0] = 0;

        char qc = query[static_cast<size_t>(i - 1)];

        for (int32_t j = 1; j <= n; ++j) {
            int32_t match_score = score_match_mismatch(qc, ref[static_cast<size_t>(j - 1)], config_);

            int32_t up   = F_prev[static_cast<size_t>(j)] + gap_extend;
            int32_t up2  = H_prev[static_cast<size_t>(j)] + gap_open + gap_extend;
            F[static_cast<size_t>(j)] = std::max(up, up2);

            int32_t left = E[static_cast<size_t>(j - 1)] + gap_extend;
            int32_t left2 = H[static_cast<size_t>(j - 1)] + gap_open + gap_extend;
            E[static_cast<size_t>(j)] = std::max(left, left2);

            int32_t diag = H_prev[static_cast<size_t>(j - 1)] + match_score;

            int32_t h = std::max(std::max(diag, E[static_cast<size_t>(j)]),
                                 std::max(F[static_cast<size_t>(j)], 0));

            int8_t op = TRACE_NONE;
            if (h == diag && h > 0) {
                op = TRACE_MATCH;
            } else if (h == E[static_cast<size_t>(j)] && h > 0) {
                op = TRACE_DEL;
            } else if (h == F[static_cast<size_t>(j)] && h > 0) {
                op = TRACE_INS;
            }
            trace[static_cast<size_t>(i)][static_cast<size_t>(j)] = op;

            H[static_cast<size_t>(j)] = h;

            if (h > best_score) {
                best_score = h;
                best_i = i;
                best_j = j;
            }
        }

        std::swap(H, H_prev);
        std::swap(F, F_prev);
        std::fill(H.begin(), H.end(), 0);
        std::fill(F.begin(), F.end(), 0);
    }

    if (best_score < config_.score_threshold) {
        return result;
    }

    std::vector<int8_t> trace_ops;
    int32_t i = best_i;
    int32_t j = best_j;

    int32_t nm = 0;

    while (i > 0 && j > 0) {
        int8_t op = trace[static_cast<size_t>(i)][static_cast<size_t>(j)];
        if (op == TRACE_NONE) break;
        trace_ops.push_back(op);
        switch (op) {
            case TRACE_MATCH:
                if (query[static_cast<size_t>(i - 1)] != ref[static_cast<size_t>(j - 1)]) {
                    ++nm;
                }
                --i; --j;
                break;
            case TRACE_DEL:
                --j;
                break;
            case TRACE_INS:
                --i;
                break;
            default:
                i = 0; j = 0;
                break;
        }
    }

    std::reverse(trace_ops.begin(), trace_ops.end());

    int32_t n_ins, n_del;
    std::string cigar = build_cigar(trace_ops, m, n, nm, n_ins, n_del);

    result.score = best_score;
    result.ref_start = j;
    result.ref_end = best_j - 1;
    result.query_start = i;
    result.query_end = best_i - 1;
    result.num_mismatches = nm;
    result.num_ins = n_ins;
    result.num_del = n_del;
    result.cigar = cigar;
    result.valid = true;

    return result;
}

#if defined(__AVX2__) || defined(_MSC_VER)

SWResult SmithWaterman::align_avx2(const std::string& query,
                                    const std::string& ref,
                                    int32_t ref_start,
                                    int32_t ref_end) const {
    int32_t n = static_cast<int32_t>(ref.size());
    int32_t actual_end = (ref_end < 0 || ref_end > n) ? n : ref_end;
    int32_t actual_start = std::max(0, ref_start);

    if (actual_end - actual_start < 32 || (int32_t)query.size() < 4) {
        std::string ref_sub = ref.substr(static_cast<size_t>(actual_start),
                                          static_cast<size_t>(actual_end - actual_start));
        SWResult r = align_basic(query, ref_sub);
        if (r.valid) {
            r.ref_start += actual_start;
            r.ref_end += actual_start;
        }
        return r;
    }

    std::string ref_sub = ref.substr(static_cast<size_t>(actual_start),
                                      static_cast<size_t>(actual_end - actual_start));
    SWResult r = align_basic(query, ref_sub);
    if (r.valid) {
        r.ref_start += actual_start;
        r.ref_end += actual_start;
    }
    return r;
}

#else

SWResult SmithWaterman::align_avx2(const std::string& query,
                                    const std::string& ref,
                                    int32_t ref_start,
                                    int32_t ref_end) const {
    int32_t n = static_cast<int32_t>(ref.size());
    int32_t actual_end = (ref_end < 0 || ref_end > n) ? n : ref_end;
    int32_t actual_start = std::max(0, ref_start);
    std::string ref_sub = ref.substr(static_cast<size_t>(actual_start),
                                      static_cast<size_t>(actual_end - actual_start));
    SWResult r = align_basic(query, ref_sub);
    if (r.valid) {
        r.ref_start += actual_start;
        r.ref_end += actual_start;
    }
    return r;
}

#endif

SWResult SmithWaterman::align(const std::string& query,
                              const std::string& ref,
                              int32_t ref_start,
                              int32_t ref_end) const {
    if (has_avx2()) {
        return align_avx2(query, ref, ref_start, ref_end);
    }
    int32_t n = static_cast<int32_t>(ref.size());
    int32_t actual_end = (ref_end < 0 || ref_end > n) ? n : ref_end;
    int32_t actual_start = std::max(0, ref_start);
    std::string ref_sub = ref.substr(static_cast<size_t>(actual_start),
                                      static_cast<size_t>(actual_end - actual_start));
    SWResult r = align_basic(query, ref_sub);
    if (r.valid) {
        r.ref_start += actual_start;
        r.ref_end += actual_start;
    }
    return r;
}

}
