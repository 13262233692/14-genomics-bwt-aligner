#pragma once

#include "sequence.h"
#include <string>
#include <vector>
#include <fstream>
#include <mutex>
#include <cstdint>
#include <atomic>

namespace bwt_aligner {

class FastqReader {
public:
    FastqReader();
    ~FastqReader();

    bool open(const std::string& filename, bool is_first_end = true);
    void close();

    bool read_next(FastqRead& read);
    size_t read_batch(std::vector<FastqRead>& reads, size_t batch_size);

    uint64_t total_reads() const { return total_reads_.load(); }
    bool is_open() const { return is_open_; }

private:
    std::ifstream file_;
    std::mutex mutex_;
    std::atomic<uint64_t> total_reads_;
    int32_t current_pair_id_;
    bool is_first_end_;
    bool is_open_;
};

class PairedFastqReader {
public:
    PairedFastqReader();
    ~PairedFastqReader();

    bool open(const std::string& filename_r1, const std::string& filename_r2);
    void close();

    size_t read_batch_paired(std::vector<FastqRead>& reads, size_t batch_size);
    uint64_t total_reads() const { return total_reads_.load(); }

private:
    FastqReader reader1_;
    FastqReader reader2_;
    std::mutex mutex_;
    std::atomic<uint64_t> total_reads_;
    bool is_open_;
};

}
