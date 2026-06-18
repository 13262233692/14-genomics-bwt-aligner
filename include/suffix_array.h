#pragma once

#include <string>
#include <vector>
#include <cstdint>

namespace bwt_aligner {

class SuffixArray {
public:
    SuffixArray();
    ~SuffixArray();

    void build(const std::string& text);

    const std::vector<int64_t>& get_sa() const { return sa_; }
    int64_t operator[](size_t i) const { return sa_[i]; }
    size_t size() const { return sa_.size(); }

    void clear();

private:
    std::vector<int64_t> sa_;

    void sais(const int64_t* s, int64_t* sa, int64_t n, int64_t k);
};

}
