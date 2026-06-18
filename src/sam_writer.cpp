#include "sam_writer.h"
#include <sstream>
#include <iomanip>

namespace bwt_aligner {

SamWriter::SamWriter()
    : total_written_(0)
    , is_open_(false)
    , header_written_(false) {
}

SamWriter::~SamWriter() {
    close();
}

bool SamWriter::open(const std::string& filename) {
    close();
    file_.open(filename, std::ios::out);
    if (file_.is_open()) {
        is_open_ = true;
        total_written_ = 0;
        header_written_ = false;
        return true;
    }
    return false;
}

void SamWriter::close() {
    if (is_open_) {
        file_.close();
        is_open_ = false;
    }
}

void SamWriter::write_header(const std::vector<ReferenceSequence>& references,
                             const std::string& program_id,
                             const std::string& program_version) {
    if (!is_open_) return;
    std::lock_guard<std::mutex> lock(mutex_);

    file_ << "@HD\tVN:1.6\tSO:coordinate\n";

    for (const auto& ref : references) {
        file_ << "@SQ\tSN:" << ref.name << "\tLN:" << ref.length << "\n";
    }

    file_ << "@PG\tID:" << program_id
          << "\tPN:" << program_id
          << "\tVN:" << program_version << "\n";

    header_written_ = true;
}

uint32_t SamWriter::compute_flag(const AlignmentResult& aln) const {
    uint32_t flag = 0;

    if (aln.is_paired) flag |= 0x1;
    if (aln.is_proper_pair) flag |= 0x2;
    if (aln.is_unmapped) flag |= 0x4;
    if (aln.mate_is_unmapped) flag |= 0x8;
    if (aln.is_reverse) flag |= 0x10;
    if (aln.mate_is_reverse) flag |= 0x20;
    if (aln.is_paired && aln.is_first) flag |= 0x40;
    if (aln.is_paired && !aln.is_first) flag |= 0x80;
    if (aln.is_secondary) flag |= 0x100;
    if (aln.is_supplementary) flag |= 0x800;

    return flag;
}

std::string SamWriter::format_cigar(int64_t read_length) const {
    if (read_length <= 0) return "*";
    return std::to_string(read_length) + "M";
}

void SamWriter::write_alignment(const AlignmentResult& aln) {
    if (!is_open_) return;
    std::lock_guard<std::mutex> lock(mutex_);

    std::stringstream ss;

    ss << aln.read_id << "\t";
    ss << compute_flag(aln) << "\t";

    if (aln.is_unmapped) {
        ss << "*\t0\t0\t*\t*\t0\t0\t";
    } else {
        ss << aln.ref_name << "\t";
        ss << (aln.position + 1) << "\t";
        ss << aln.mapq << "\t";
        ss << (aln.cigar.empty() ? format_cigar(static_cast<int64_t>(aln.sequence.size())) : aln.cigar) << "\t";

        if (aln.is_paired) {
            if (aln.mate_is_unmapped) {
                ss << "*\t0\t0\t";
            } else {
                ss << (aln.ref_name == aln.ref_name ? "=" : (aln.mate_ref_id >= 0 ? aln.ref_name : "*")) << "\t";
                ss << (aln.mate_position + 1) << "\t";
                ss << aln.insert_size << "\t";
            }
        } else {
            ss << "*\t0\t0\t";
        }
    }

    if (aln.is_reverse) {
        ss << reverse_complement(aln.sequence) << "\t";
        std::string rev_qual(aln.quality.rbegin(), aln.quality.rend());
        ss << rev_qual;
    } else {
        ss << aln.sequence << "\t";
        ss << aln.quality;
    }

    if (!aln.is_unmapped) {
        ss << "\tNM:i:" << aln.num_mismatches;
    }

    ss << "\n";

    file_ << ss.str();
    ++total_written_;
}

void SamWriter::write_alignments(const std::vector<AlignmentResult>& alns) {
    for (const auto& aln : alns) {
        write_alignment(aln);
    }
}

}
