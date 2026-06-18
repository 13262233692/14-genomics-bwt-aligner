#pragma once

#include "fm_index.h"
#include "backward_search.h"
#include "fasta_reader.h"
#include "fastq_reader.h"
#include "sam_writer.h"
#include "sequence.h"
#include "smith_waterman.h"
#include <string>
#include <vector>
#include <thread>
#include <atomic>
#include <memory>

namespace bwt_aligner {

struct AlignerConfig {
    int32_t num_threads;
    int32_t seed_length;
    int32_t max_mismatches;
    int32_t min_seed_hits;
    int32_t batch_size;
    int32_t min_mapq;
    int32_t sw_flank;
    int32_t sw_score_threshold;
    bool enable_sw_fallback;

    AlignerConfig()
        : num_threads(std::thread::hardware_concurrency())
        , seed_length(16)
        , max_mismatches(2)
        , min_seed_hits(1)
        , batch_size(10000)
        , min_mapq(0)
        , sw_flank(50)
        , sw_score_threshold(8)
        , enable_sw_fallback(true)
    {}
};

class Aligner {
public:
    Aligner();
    ~Aligner();

    void set_config(const AlignerConfig& config);
    const AlignerConfig& get_config() const { return config_; }

    bool build_index(const std::string& fasta_filename);
    bool load_reference(const std::string& fasta_filename);

    bool align_single(const std::string& fastq_filename,
                      const std::string& sam_filename);

    bool align_paired(const std::string& fastq_r1,
                      const std::string& fastq_r2,
                      const std::string& sam_filename);

    const std::vector<ReferenceSequence>& get_references() const { return references_; }
    const FMIndex& get_fm_index() const { return fm_index_; }

    int64_t get_ref_id_from_position(int64_t global_pos, int64_t& local_pos) const;

    AlignmentResult align_read(const FastqRead& read) const;
    std::vector<AlignmentResult> align_read_batch(const std::vector<FastqRead>& reads) const;

    size_t memory_usage() const;

private:
    AlignerConfig config_;
    FMIndex fm_index_;
    BackwardSearch backward_search_;
    std::vector<ReferenceSequence> references_;
    std::vector<int64_t> ref_offsets_;
    mutable SmithWaterman smith_waterman_;
    bool index_built_;

    void worker_thread_single(FastqReader* reader,
                              SamWriter* writer,
                              std::atomic<bool>* stop_flag);

    void worker_thread_paired(PairedFastqReader* reader,
                              SamWriter* writer,
                              std::atomic<bool>* stop_flag);

    bool resolve_pair(AlignmentResult& aln1, AlignmentResult& aln2,
                      int64_t max_insert_size = 1000) const;
};

}
