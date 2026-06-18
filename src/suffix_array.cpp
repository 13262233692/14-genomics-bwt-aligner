#include "suffix_array.h"
#include <algorithm>
#include <cstring>

namespace bwt_aligner {

SuffixArray::SuffixArray() {}
SuffixArray::~SuffixArray() {}

static void sais_get_buckets(const int64_t* s, int64_t* bkt, int64_t n, int64_t k, bool end) {
    std::memset(bkt, 0, static_cast<size_t>(k) * sizeof(int64_t));
    for (int64_t i = 0; i < n; ++i) ++bkt[static_cast<size_t>(s[i])];
    int64_t sum = 0;
    for (int64_t i = 0; i < k; ++i) {
        sum += bkt[i];
        bkt[i] = end ? sum : sum - bkt[i];
    }
}

static void sais_induce_l(int64_t* sa, const int64_t* s, const uint8_t* t,
                           int64_t* bkt, int64_t n, int64_t k) {
    sais_get_buckets(s, bkt, n, k, false);
    for (int64_t i = 0; i < n; ++i) {
        int64_t j = sa[i] - 1;
        if (j >= 0 && !t[j]) {
            sa[bkt[s[j]]++] = j;
        }
    }
}

static void sais_induce_s(int64_t* sa, const int64_t* s, const uint8_t* t,
                           int64_t* bkt, int64_t n, int64_t k) {
    sais_get_buckets(s, bkt, n, k, true);
    for (int64_t i = n - 1; i >= 0; --i) {
        int64_t j = sa[i] - 1;
        if (j >= 0 && t[j]) {
            sa[--bkt[s[j]]] = j;
        }
    }
}

void SuffixArray::sais(const int64_t* s, int64_t* sa, int64_t n, int64_t k) {
    if (n <= 1) {
        if (n == 1) sa[0] = 0;
        return;
    }

    std::vector<uint8_t> t(static_cast<size_t>(n), 0);
    t[static_cast<size_t>(n - 1)] = 1;
    for (int64_t i = n - 2; i >= 0; --i) {
        if (s[i] < s[i + 1]) t[i] = 1;
        else if (s[i] == s[i + 1]) t[i] = t[i + 1];
        else t[i] = 0;
    }

    std::vector<int64_t> bkt(static_cast<size_t>(k));

    std::fill(sa, sa + n, static_cast<int64_t>(-1));
    sais_get_buckets(s, bkt.data(), n, k, true);
    for (int64_t i = 1; i < n; ++i) {
        if (t[i] && !t[i - 1]) {
            sa[--bkt[s[i]]] = i;
        }
    }

    sais_induce_l(sa, s, t.data(), bkt.data(), n, k);
    sais_induce_s(sa, s, t.data(), bkt.data(), n, k);

    int64_t n1 = 0;
    for (int64_t i = 0; i < n; ++i) {
        if (sa[i] > 0 && t[sa[i]] && !t[sa[i] - 1]) {
            sa[n1++] = sa[i];
        }
    }

    if (n1 == 0) return;

    std::fill(sa + n1, sa + n, static_cast<int64_t>(-1));

    int64_t name = 0;
    int64_t prev = -1;
    for (int64_t i = 0; i < n1; ++i) {
        int64_t pos = sa[i];
        bool diff = false;
        if (prev < 0) {
            diff = true;
        } else {
            for (int64_t d = 0; ; ++d) {
                if (s[pos + d] != s[prev + d]) { diff = true; break; }
                bool pos_lms = (d > 0 && t[pos + d] && !t[pos + d - 1]);
                bool prev_lms = (d > 0 && t[prev + d] && !t[prev + d - 1]);
                if (pos_lms || prev_lms) {
                    diff = (pos_lms != prev_lms);
                    break;
                }
            }
        }
        if (diff) { ++name; prev = pos; }
        int64_t p = (pos % 2 == 0) ? pos / 2 : (pos - 1) / 2;
        sa[n1 + p] = name - 1;
    }

    int64_t j_idx = n - 1;
    for (int64_t i = n - 1; i >= n1; --i) {
        if (sa[i] >= 0) {
            sa[j_idx--] = sa[i];
        }
    }

    int64_t* sa1 = sa;
    int64_t* s1 = sa + n - n1;

    if (name < n1) {
        this->sais(s1, sa1, n1, name);
    } else {
        for (int64_t i = 0; i < n1; ++i) {
            sa1[s1[i]] = i;
        }
    }

    std::vector<int64_t> p1(static_cast<size_t>(n1));
    int64_t cnt = 0;
    for (int64_t i = 1; i < n; ++i) {
        if (t[i] && !t[i - 1]) {
            p1[static_cast<size_t>(cnt++)] = i;
        }
    }

    for (int64_t i = 0; i < n1; ++i) {
        sa1[i] = p1[static_cast<size_t>(sa1[i])];
    }

    std::fill(sa + n1, sa + n, static_cast<int64_t>(-1));
    sais_get_buckets(s, bkt.data(), n, k, true);
    for (int64_t i = n1 - 1; i >= 0; --i) {
        int64_t jj = sa[i];
        sa[i] = -1;
        sa[--bkt[s[jj]]] = jj;
    }

    sais_induce_l(sa, s, t.data(), bkt.data(), n, k);
    sais_induce_s(sa, s, t.data(), bkt.data(), n, k);
}

void SuffixArray::build(const std::string& text) {
    clear();

    int64_t n = static_cast<int64_t>(text.size());
    if (n == 0) return;

    int64_t k = 256;
    std::vector<int64_t> s_int(static_cast<size_t>(n));
    for (int64_t i = 0; i < n; ++i) {
        s_int[static_cast<size_t>(i)] = static_cast<int64_t>(
            static_cast<unsigned char>(text[static_cast<size_t>(i)]));
    }

    sa_.resize(static_cast<size_t>(n), -1);
    sais(s_int.data(), sa_.data(), n, k);
}

void SuffixArray::clear() {
    sa_.clear();
    sa_.shrink_to_fit();
}

}
