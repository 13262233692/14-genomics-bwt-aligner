#pragma once

#include "sequence.h"
#include <string>
#include <vector>
#include <fstream>
#include <mutex>
#include <cstdint>

namespace bwt_aligner {

class SamWriter {
public:
    SamWriter();
    ~SamWriter();

    bool open(const std::string& filename);
    void close();

    void write_header(const std::vector<ReferenceSequence>& references,
                      const std::string& program_id = "bwt_aligner",
                      const std::string& program_version = "1.0.0");

    void write_alignment(const AlignmentResult& aln);
    void write_alignments(const std::vector<AlignmentResult>& alns);

    uint64_t total_written() const { return total_written_; }

private:
    std::ofstream file_;
    std::mutex mutex_;
    uint64_t total_written_;
    bool is_open_;
    bool header_written_;

    uint32_t compute_flag(const AlignmentResult& aln) const;
    std::string format_cigar(int64_t read_length) const;
};

}
