#include "fastq_reader.h"
#include <algorithm>

namespace bwt_aligner {

FastqReader::FastqReader()
    : total_reads_(0)
    , current_pair_id_(0)
    , is_first_end_(true)
    , is_open_(false) {
}

FastqReader::~FastqReader() {
    close();
}

bool FastqReader::open(const std::string& filename, bool is_first_end) {
    close();
    file_.open(filename, std::ios::in);
    if (file_.is_open()) {
        is_open_ = true;
        is_first_end_ = is_first_end;
        total_reads_ = 0;
        current_pair_id_ = 0;
        return true;
    }
    return false;
}

void FastqReader::close() {
    if (is_open_) {
        file_.close();
        is_open_ = false;
    }
}

bool FastqReader::read_next(FastqRead& read) {
    if (!is_open_) return false;

    std::lock_guard<std::mutex> lock(mutex_);

    std::string line1, line2, line3, line4;
    if (!std::getline(file_, line1)) return false;
    if (!std::getline(file_, line2)) return false;
    if (!std::getline(file_, line3)) return false;
    if (!std::getline(file_, line4)) return false;

    if (!line1.empty() && line1[0] == '@') {
        read.id = line1.substr(1);
    } else {
        read.id = line1;
    }

    size_t space_pos = read.id.find_first_of(" \t");
    if (space_pos != std::string::npos) {
        read.id = read.id.substr(0, space_pos);
    }

    if (read.id.size() >= 2 && (read.id.substr(read.id.size() - 2) == "/1" || read.id.substr(read.id.size() - 2) == "/2")) {
        read.id.resize(read.id.size() - 2);
    }

    read.sequence = line2;
    read.quality = line4;
    read.is_first = is_first_end_;
    read.pair_id = current_pair_id_++;

    ++total_reads_;

    return true;
}

size_t FastqReader::read_batch(std::vector<FastqRead>& reads, size_t batch_size) {
    reads.clear();
    reads.reserve(batch_size);

    FastqRead read;
    size_t count = 0;

    while (count < batch_size && read_next(read)) {
        reads.push_back(std::move(read));
        ++count;
    }

    return count;
}

PairedFastqReader::PairedFastqReader()
    : total_reads_(0)
    , is_open_(false) {
}

PairedFastqReader::~PairedFastqReader() {
    close();
}

bool PairedFastqReader::open(const std::string& filename_r1, const std::string& filename_r2) {
    close();
    total_reads_.store(0);
    bool r1_ok = reader1_.open(filename_r1, true);
    bool r2_ok = reader2_.open(filename_r2, false);
    is_open_ = r1_ok && r2_ok;
    return is_open_;
}

void PairedFastqReader::close() {
    reader1_.close();
    reader2_.close();
    is_open_ = false;
}

size_t PairedFastqReader::read_batch_paired(std::vector<FastqRead>& reads, size_t batch_size) {
    reads.clear();
    reads.reserve(batch_size);

    std::lock_guard<std::mutex> lock(mutex_);

    size_t count = 0;
    FastqRead r1, r2;
    size_t pair_count = 0;
    size_t target_pairs = batch_size / 2;

    while (pair_count < target_pairs) {
        bool has_r1 = reader1_.read_next(r1);
        bool has_r2 = reader2_.read_next(r2);

        if (has_r1 && has_r2) {
            r2.pair_id = r1.pair_id;
            reads.push_back(std::move(r1));
            reads.push_back(std::move(r2));
            ++pair_count;
            count += 2;
            total_reads_ += 2;
        } else {
            break;
        }
    }

    return count;
}

}
