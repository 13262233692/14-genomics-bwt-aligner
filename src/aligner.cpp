#include "aligner.h"
#include <iostream>
#include <thread>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <unordered_map>

namespace bwt_aligner {

Aligner::Aligner()
    : index_built_(false) {
    SWConfig sw_cfg;
    smith_waterman_.set_config(sw_cfg);
}

Aligner::~Aligner() {
}

void Aligner::set_config(const AlignerConfig& config) {
    config_ = config;
    SWConfig sw_cfg;
    sw_cfg.score_threshold = config_.sw_score_threshold;
    smith_waterman_.set_config(sw_cfg);
}

int64_t Aligner::get_ref_id_from_position(int64_t global_pos, int64_t& local_pos) const {
    if (ref_offsets_.empty()) { local_pos = -1; return -1; }

    auto it = std::upper_bound(ref_offsets_.begin(), ref_offsets_.end(), global_pos);
    if (it == ref_offsets_.begin()) { local_pos = -1; return -1; }

    int64_t idx = static_cast<int64_t>(it - ref_offsets_.begin()) - 1;
    if (idx < 0 || idx >= static_cast<int64_t>(references_.size())) {
        local_pos = -1;
        return -1;
    }

    const auto& ref = references_[static_cast<size_t>(idx)];
    local_pos = global_pos - ref.global_offset;

    if (local_pos < 0 || local_pos >= static_cast<int64_t>(ref.length)) {
        local_pos = -1;
        return -1;
    }

    return idx;
}

bool Aligner::load_reference(const std::string& fasta_filename) {
    FastaReader reader;
    if (!reader.open(fasta_filename)) {
        std::cerr << "Error: Cannot open FASTA file: " << fasta_filename << std::endl;
        return false;
    }

    if (!reader.read_all(references_)) {
        std::cerr << "Error: Cannot read sequences from FASTA file" << std::endl;
        return false;
    }

    std::string concatenated_text;
    concatenated_text.reserve(
        std::accumulate(references_.begin(), references_.end(), size_t(0),
            [](size_t s, const ReferenceSequence& r) { return s + r.length; })
        + references_.size()
    );

    ref_offsets_.clear();
    ref_offsets_.reserve(references_.size());

    for (size_t i = 0; i < references_.size(); ++i) {
        references_[i].global_offset = concatenated_text.size();
        ref_offsets_.push_back(static_cast<int64_t>(concatenated_text.size()));
        concatenated_text += references_[i].sequence;
        if (i < references_.size() - 1) {
            concatenated_text += '#';
        }
    }

    std::cout << "Loaded " << references_.size() << " reference sequences, "
              << "total length: " << concatenated_text.size() << " bp" << std::endl;

    fm_index_.build(concatenated_text);
    backward_search_.set_fm_index(&fm_index_);

    concatenated_text.clear();
    concatenated_text.shrink_to_fit();

    index_built_ = true;
    return true;
}

bool Aligner::build_index(const std::string& fasta_filename) {
    std::cout << "Loading reference genome and building index..." << std::endl;
    if (!load_reference(fasta_filename)) {
        return false;
    }

    std::cout << "FM-Index built successfully." << std::endl;
    std::cout << "Memory usage: " << (fm_index_.memory_usage() / (1024.0 * 1024.0))
              << " MB" << std::endl;
    return true;
}

AlignmentResult Aligner::align_read(const FastqRead& read) const {
    AlignmentResult aln;
    aln.read_id = read.id;
    aln.sequence = read.sequence;
    aln.quality = read.quality;
    aln.is_first = read.is_first;
    aln.is_paired = false;
    aln.is_unmapped = true;
    aln.mate_is_unmapped = true;
    aln.is_reverse = false;
    aln.mapq = 0;
    aln.position = -1;
    aln.ref_id = -1;
    aln.ref_name = "*";
    aln.num_mismatches = 0;
    aln.cigar = "*";
    aln.is_secondary = false;
    aln.is_supplementary = false;
    aln.is_proper_pair = false;
    aln.mate_is_reverse = false;
    aln.mate_ref_id = -1;
    aln.mate_position = -1;
    aln.insert_size = 0;

    if (!index_built_ || read.sequence.empty()) {
        return aln;
    }

    int32_t read_len = static_cast<int32_t>(read.sequence.size());
    int32_t seed_len = std::min(config_.seed_length, read_len);

    int32_t best_mismatches = config_.max_mismatches + 1;

    {
        SearchRange range = backward_search_.search(read.sequence);
        if (range.valid() && range.count() > 0 && range.count() < 10000) {
            for (int64_t i = range.low; i <= range.high; ++i) {
                int64_t global_pos = fm_index_.resolve_position(i);
                int64_t local_pos;
                int64_t ref_id = get_ref_id_from_position(global_pos, local_pos);

                if (ref_id >= 0) {
                    int32_t mm;
                    if (backward_search_.verify_alignment(read.sequence, global_pos, false,
                                                          config_.max_mismatches, mm)) {
                        aln.is_unmapped = false;
                        aln.is_reverse = false;
                        aln.position = local_pos;
                        aln.ref_id = static_cast<int32_t>(ref_id);
                        aln.ref_name = references_[static_cast<size_t>(ref_id)].name;
                        aln.num_mismatches = mm;
                        aln.cigar = std::to_string(read_len) + "M";
                        aln.mapq = static_cast<int32_t>(std::min(60.0,
                            40.0 - 10.0 * std::log10(std::max(1.0, (double)range.count()))));
                        if (aln.mapq < 0) aln.mapq = 0;
                        return aln;
                    }
                }
            }
        }
    }

    {
        std::string rc = reverse_complement(read.sequence);
        SearchRange range = backward_search_.search(rc);
        if (range.valid() && range.count() > 0 && range.count() < 10000) {
            for (int64_t i = range.low; i <= range.high; ++i) {
                int64_t global_pos = fm_index_.resolve_position(i);
                int64_t local_pos;
                int64_t ref_id = get_ref_id_from_position(global_pos, local_pos);

                if (ref_id >= 0) {
                    int32_t mm;
                    if (backward_search_.verify_alignment(read.sequence, global_pos, true,
                                                          config_.max_mismatches, mm)) {
                        aln.is_unmapped = false;
                        aln.is_reverse = true;
                        aln.position = local_pos;
                        aln.ref_id = static_cast<int32_t>(ref_id);
                        aln.ref_name = references_[static_cast<size_t>(ref_id)].name;
                        aln.num_mismatches = mm;
                        aln.cigar = std::to_string(read_len) + "M";
                        aln.mapq = static_cast<int32_t>(std::min(60.0,
                            40.0 - 10.0 * std::log10(std::max(1.0, (double)range.count()))));
                        if (aln.mapq < 0) aln.mapq = 0;
                        return aln;
                    }
                }
            }
        }
    }

    {
        std::vector<SeedHit> seeds = backward_search_.find_seeds_both_strands(
            read.sequence, seed_len, 0);

        std::unordered_map<int64_t, std::vector<SeedHit>> pos_groups;
        for (const auto& seed : seeds) {
            int64_t aligned_pos = seed.ref_position - seed.query_start;
            pos_groups[aligned_pos].push_back(seed);
        }

        for (const auto& kv : pos_groups) {
            if (static_cast<int32_t>(kv.second.size()) >= config_.min_seed_hits) {
                int64_t global_pos = kv.first;
                bool is_rev = kv.second.front().is_reverse;

                int64_t local_pos;
                int64_t ref_id = get_ref_id_from_position(global_pos, local_pos);

                if (ref_id >= 0) {
                    int32_t mm;
                    if (backward_search_.verify_alignment(read.sequence, global_pos, is_rev,
                                                          config_.max_mismatches, mm)) {
                        if (mm < best_mismatches) {
                            best_mismatches = mm;
                            aln.is_unmapped = false;
                            aln.is_reverse = is_rev;
                            aln.position = local_pos;
                            aln.ref_id = static_cast<int32_t>(ref_id);
                            aln.ref_name = references_[static_cast<size_t>(ref_id)].name;
                            aln.num_mismatches = mm;
                            aln.cigar = std::to_string(read_len) + "M";
                            aln.mapq = static_cast<int32_t>(std::min(60.0,
                                30.0 - 10.0 * std::log10(std::max(1.0, (double)mm + 1))));
                            if (aln.mapq < 0) aln.mapq = 0;
                        }
                    }
                }
            }
        }
    }

    if (config_.enable_sw_fallback && aln.is_unmapped && !read.sequence.empty()) {
        std::vector<SeedHit> seeds = backward_search_.find_seeds_both_strands(
            read.sequence, seed_len, 0);

        std::unordered_map<int64_t, std::vector<SeedHit>> pos_groups;
        for (const auto& seed : seeds) {
            int64_t aligned_pos = seed.ref_position - seed.query_start;
            pos_groups[aligned_pos].push_back(seed);
        }

        int32_t best_score = config_.sw_score_threshold - 1;
        SWResult best_sw;
        int64_t best_global_pos = -1;
        int64_t best_ref_id = -1;
        bool best_is_rev = false;
        int64_t best_local_pos = -1;

        const ReferenceSequence* best_ref = nullptr;

        for (const auto& kv : pos_groups) {
            if (static_cast<int32_t>(kv.second.size()) >= config_.min_seed_hits) {
                int64_t global_pos = kv.first;
                bool is_rev = kv.second.front().is_reverse;

                int64_t local_pos;
                int64_t ref_id = get_ref_id_from_position(global_pos, local_pos);

                if (ref_id < 0) continue;

                const auto& ref = references_[static_cast<size_t>(ref_id)];

                int32_t flank = config_.sw_flank;
                int32_t ref_start = std::max(static_cast<int32_t>(0),
                    static_cast<int32_t>(local_pos) - flank);
                int32_t ref_end = std::min(static_cast<int32_t>(ref.length),
                    static_cast<int32_t>(local_pos) + read_len + flank);

                if (ref_end <= ref_start) continue;

                std::string query_seq = is_rev ?
                    reverse_complement(read.sequence) : read.sequence;

                SWResult sw = smith_waterman_.align(
                    query_seq, ref.sequence, ref_start, ref_end);

                if (sw.valid && sw.score > best_score) {
                    best_score = sw.score;
                    best_sw = sw;
                    best_global_pos = global_pos;
                    best_ref_id = ref_id;
                    best_is_rev = is_rev;
                    best_local_pos = sw.ref_start;
                    best_ref = &ref;
                }
            }
        }

        if (best_sw.valid && best_ref) {
            int32_t total_err = best_sw.num_mismatches + best_sw.num_ins + best_sw.num_del;
            int32_t aligned_len = best_sw.query_end - best_sw.query_start + 1;

            aln.is_unmapped = false;
            aln.is_reverse = best_is_rev;
            aln.position = best_local_pos;
            aln.ref_id = static_cast<int32_t>(best_ref_id);
            aln.ref_name = best_ref->name;
            aln.num_mismatches = total_err;
            aln.cigar = best_sw.cigar;
            aln.mapq = static_cast<int32_t>(std::min(60.0,
                25.0 - 10.0 * std::log10(std::max(1.0, (double)total_err + 1))));
            if (aln.mapq < 0) aln.mapq = 0;

            if (best_sw.query_start > 0 || best_sw.query_end + 1 < read_len) {
                std::string new_cigar;
                if (best_sw.query_start > 0) {
                    new_cigar += std::to_string(best_sw.query_start) + "S";
                }
                new_cigar += best_sw.cigar;
                if (best_sw.query_end + 1 < read_len) {
                    new_cigar += std::to_string(read_len - best_sw.query_end - 1) + "S";
                }
                aln.cigar = new_cigar;
            }
        }
    }

    return aln;
}

std::vector<AlignmentResult> Aligner::align_read_batch(const std::vector<FastqRead>& reads) const {
    std::vector<AlignmentResult> results;
    results.reserve(reads.size());

    for (const auto& read : reads) {
        results.push_back(align_read(read));
    }

    return results;
}

bool Aligner::resolve_pair(AlignmentResult& aln1, AlignmentResult& aln2,
                           int64_t max_insert_size) const {
    if (aln1.is_unmapped || aln2.is_unmapped) {
        return false;
    }

    if (aln1.ref_id != aln2.ref_id) {
        return false;
    }

    aln1.is_paired = true;
    aln2.is_paired = true;
    aln1.mate_ref_id = aln2.ref_id;
    aln2.mate_ref_id = aln1.ref_id;
    aln1.mate_position = aln2.position;
    aln2.mate_position = aln1.position;
    aln1.mate_is_reverse = aln2.is_reverse;
    aln2.mate_is_reverse = aln1.is_reverse;
    aln1.mate_is_unmapped = false;
    aln2.mate_is_unmapped = false;

    int64_t pos1 = aln1.position;
    int64_t pos2 = aln2.position;
    int64_t end1 = pos1 + aln1.sequence.size();
    int64_t end2 = pos2 + aln2.sequence.size();

    int64_t insert_start = std::min(pos1, pos2);
    int64_t insert_end = std::max(end1, end2);
    int64_t insert_size = insert_end - insert_start;

    aln1.insert_size = static_cast<int32_t>(aln1.position < aln2.position ? insert_size : -insert_size);
    aln2.insert_size = static_cast<int32_t>(aln2.position < aln1.position ? insert_size : -insert_size);

    bool proper = (aln1.is_reverse != aln2.is_reverse);
    if (proper && insert_size > 0 && insert_size <= max_insert_size) {
        aln1.is_proper_pair = true;
        aln2.is_proper_pair = true;
        return true;
    }

    return false;
}

void Aligner::worker_thread_single(FastqReader* reader,
                                   SamWriter* writer,
                                   std::atomic<bool>* stop_flag) {
    std::vector<FastqRead> batch;

    while (!stop_flag->load()) {
        size_t count = reader->read_batch(batch, static_cast<size_t>(config_.batch_size));
        if (count == 0) break;

        auto results = align_read_batch(batch);
        writer->write_alignments(results);
    }
}

void Aligner::worker_thread_paired(PairedFastqReader* reader,
                                   SamWriter* writer,
                                   std::atomic<bool>* stop_flag) {
    std::vector<FastqRead> batch;

    while (!stop_flag->load()) {
        size_t count = reader->read_batch_paired(batch, static_cast<size_t>(config_.batch_size));
        if (count == 0) break;

        std::vector<AlignmentResult> results;
        results.reserve(batch.size());

        for (size_t i = 0; i < batch.size(); i += 2) {
            if (i + 1 >= batch.size()) break;

            const FastqRead& r1 = batch[i];
            const FastqRead& r2 = batch[i + 1];

            AlignmentResult aln1 = align_read(r1);
            AlignmentResult aln2 = align_read(r2);

            aln1.is_paired = true;
            aln2.is_paired = true;
            aln1.is_first = true;
            aln2.is_first = false;

            resolve_pair(aln1, aln2);

            results.push_back(std::move(aln1));
            results.push_back(std::move(aln2));
        }

        writer->write_alignments(results);
    }
}

bool Aligner::align_single(const std::string& fastq_filename,
                           const std::string& sam_filename) {
    if (!index_built_) {
        std::cerr << "Error: Index not built. Call build_index() first." << std::endl;
        return false;
    }

    FastqReader reader;
    if (!reader.open(fastq_filename, true)) {
        std::cerr << "Error: Cannot open FASTQ file: " << fastq_filename << std::endl;
        return false;
    }

    SamWriter writer;
    if (!writer.open(sam_filename)) {
        std::cerr << "Error: Cannot create SAM file: " << sam_filename << std::endl;
        return false;
    }

    writer.write_header(references_);

    std::cout << "Starting single-end alignment with " << config_.num_threads << " threads..." << std::endl;

    std::atomic<bool> stop_flag(false);
    std::vector<std::thread> threads;

    for (int32_t t = 0; t < config_.num_threads; ++t) {
        threads.emplace_back(&Aligner::worker_thread_single, this,
                            &reader, &writer, &stop_flag);
    }

    for (auto& th : threads) {
        th.join();
    }

    reader.close();
    writer.close();

    std::cout << "Alignment complete. Total reads aligned: " << writer.total_written() << std::endl;

    return true;
}

bool Aligner::align_paired(const std::string& fastq_r1,
                           const std::string& fastq_r2,
                           const std::string& sam_filename) {
    if (!index_built_) {
        std::cerr << "Error: Index not built. Call build_index() first." << std::endl;
        return false;
    }

    PairedFastqReader reader;
    if (!reader.open(fastq_r1, fastq_r2)) {
        std::cerr << "Error: Cannot open paired FASTQ files" << std::endl;
        return false;
    }

    SamWriter writer;
    if (!writer.open(sam_filename)) {
        std::cerr << "Error: Cannot create SAM file: " << sam_filename << std::endl;
        return false;
    }

    writer.write_header(references_);

    std::cout << "Starting paired-end alignment with " << config_.num_threads << " threads..." << std::endl;

    std::atomic<bool> stop_flag(false);
    std::vector<std::thread> threads;

    for (int32_t t = 0; t < config_.num_threads; ++t) {
        threads.emplace_back(&Aligner::worker_thread_paired, this,
                            &reader, &writer, &stop_flag);
    }

    for (auto& th : threads) {
        th.join();
    }

    reader.close();
    writer.close();

    std::cout << "Alignment complete. Total reads aligned: " << writer.total_written() << std::endl;

    return true;
}

size_t Aligner::memory_usage() const {
    size_t total = fm_index_.memory_usage();
    total += ref_offsets_.capacity() * sizeof(int64_t);
    for (const auto& ref : references_) {
        total += ref.name.capacity() + ref.sequence.capacity();
    }
    return total;
}

}
