#include "fasta_reader.h"
#include <sstream>
#include <algorithm>
#include <cctype>

namespace bwt_aligner {

FastaReader::FastaReader()
    : total_length_(0)
    , num_sequences_(0)
    , is_open_(false) {
}

FastaReader::~FastaReader() {
    close();
}

bool FastaReader::open(const std::string& filename) {
    close();
    file_.open(filename, std::ios::in);
    if (file_.is_open()) {
        is_open_ = true;
        total_length_ = 0;
        num_sequences_ = 0;
        return true;
    }
    return false;
}

void FastaReader::close() {
    if (is_open_) {
        file_.close();
        is_open_ = false;
    }
}

static void trim_right(std::string& s) {
    while (!s.empty()) {
        char c = s.back();
        if (c == ' ' || c == '\t' || c == '\n' || c == '\r') {
            s.pop_back();
        } else {
            break;
        }
    }
}

bool FastaReader::read_all(std::vector<ReferenceSequence>& sequences) {
    if (!is_open_) return false;

    sequences.clear();
    std::string line;
    ReferenceSequence current;
    bool in_sequence = false;
    uint64_t global_offset = 0;

    while (std::getline(file_, line)) {
        if (line.empty()) continue;

        trim_right(line);
        if (line.empty()) continue;

        if (line[0] == '>') {
            if (in_sequence) {
                current.length = current.sequence.size();
                current.global_offset = global_offset;
                global_offset += current.length + 1;
                total_length_ += current.length;
                sequences.push_back(std::move(current));
                ++num_sequences_;
                current = ReferenceSequence();
            }

            std::string name_part = line.substr(1);
            size_t space_pos = name_part.find_first_of(" \t");
            if (space_pos != std::string::npos) {
                current.name = name_part.substr(0, space_pos);
            } else {
                current.name = name_part;
            }
            in_sequence = true;
        } else {
            for (char& c : line) {
                c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
            }
            current.sequence += line;
        }
    }

    if (in_sequence && !current.name.empty()) {
        current.length = current.sequence.size();
        current.global_offset = global_offset;
        total_length_ += current.length;
        sequences.push_back(std::move(current));
        ++num_sequences_;
    }

    return !sequences.empty();
}

}
