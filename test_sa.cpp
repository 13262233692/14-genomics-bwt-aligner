#include "suffix_array.h"
#include <iostream>
#include <cstring>
#include <string>
#include <vector>
#include <algorithm>

using namespace bwt_aligner;

static void verify_sa(const std::string& text, const std::vector<int64_t>& sa) {
    int64_t n = static_cast<int64_t>(text.size());
    std::vector<int64_t> expected(static_cast<size_t>(n));
    for (int64_t i = 0; i < n; ++i) expected[i] = i;
    std::sort(expected.begin(), expected.end(), [&text](int64_t a, int64_t b) {
        return text.substr(a) < text.substr(b);
    });

    bool ok = true;
    for (int64_t i = 0; i < n; ++i) {
        if (sa[i] != expected[i]) {
            std::cout << "MISMATCH at " << i << ": got " << sa[i] << " expected " << expected[i] << std::endl;
            ok = false;
        }
    }
    if (ok) std::cout << "SA CORRECT!" << std::endl;
}

int main() {
    std::vector<std::string> tests = {
        "a$",
        "ab$",
        "aba$",
        "ACGT$",
        "aaaaaa$",
        "ACGTACGT$",
        "ACGTACGTGATTACA$",
        "banana$",
        "mississippi$",
        "abababab$",
    };

    int pass = 0, fail = 0;
    for (const auto& text : tests) {
        std::cout << "Testing: " << text << " (len=" << text.size() << ")" << std::endl;

        std::vector<int64_t> s_int(text.size());
        for (size_t i = 0; i < text.size(); ++i) {
            s_int[i] = (int64_t)(unsigned char)text[i];
        }

        std::vector<uint8_t> t(text.size(), 0);
        t[text.size()-1] = 1;
        for (int i = (int)text.size()-2; i >= 0; --i) {
            if (s_int[i] < s_int[i+1]) t[i] = 1;
            else if (s_int[i] == s_int[i+1]) t[i] = t[i+1];
            else t[i] = 0;
        }
        std::cout << "  t = [";
        for (size_t i = 0; i < t.size(); ++i) std::cout << (int)t[i] << (i+1<t.size()?",":"");
        std::cout << "]" << std::endl;

        std::cout << "  LMS positions: ";
        for (size_t i = 1; i < text.size(); ++i) {
            if (t[i] && !t[i-1])
                std::cout << i << " ";
        }
        std::cout << std::endl;

        try {
            SuffixArray sa;
            sa.build(text);
            std::cout << "  SA size: " << sa.size() << std::endl;
            for (size_t i = 0; i < sa.size(); ++i) {
                std::cout << "  SA[" << i << "]=" << sa[i];
                if (sa[i] >= 0 && sa[i] < static_cast<int64_t>(text.size())) {
                    std::cout << "  '" << text.substr(sa[i], std::min(size_t(8), text.size() - sa[i])) << "'";
                }
                std::cout << std::endl;
            }
            std::vector<int64_t> sa_vec(sa.size());
            for (size_t i = 0; i < sa.size(); ++i) sa_vec[i] = sa[i];
            verify_sa(text, sa_vec);
            bool ok = true;
            for (size_t i = 0; i < sa_vec.size(); ++i) {
                std::vector<int64_t> expected(sa_vec.size());
                for (size_t j = 0; j < expected.size(); ++j) expected[j] = j;
                std::sort(expected.begin(), expected.end(), [&text](int64_t a, int64_t b) {
                    return text.substr(a) < text.substr(b);
                });
                if (sa_vec[i] != expected[i]) { ok = false; break; }
            }
            if (ok) ++pass; else ++fail;
        } catch (...) {
            std::cout << "  EXCEPTION!" << std::endl;
            ++fail;
        }
        std::cout << std::endl;
    }
    std::cout << "Results: " << pass << " passed, " << fail << " failed" << std::endl;
    return fail > 0 ? 1 : 0;
}
