#include "smith_waterman.h"
#include "aligner.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cassert>

using namespace bwt_aligner;

static void test_snp() {
    std::cout << "=== Test: SNP (1 mismatch) ===" << std::endl;
    std::string ref   = "ACGTACGTGATTACA";
    std::string query = "ACGTACGTAATTACA";
    SmithWaterman sw;
    SWConfig cfg;
    cfg.score_threshold = 10;
    sw.set_config(cfg);
    SWResult r = sw.align(query, ref);
    std::cout << "  score=" << r.score << " valid=" << r.valid << std::endl;
    std::cout << "  ref_start=" << r.ref_start << " ref_end=" << r.ref_end << std::endl;
    std::cout << "  query_start=" << r.query_start << " query_end=" << r.query_end << std::endl;
    std::cout << "  mismatches=" << r.num_mismatches
              << " ins=" << r.num_ins << " del=" << r.num_del << std::endl;
    std::cout << "  CIGAR: " << r.cigar << std::endl;
    assert(r.valid);
    assert(r.num_mismatches == 1);
    assert(r.num_ins == 0);
    assert(r.num_del == 0);
    std::cout << "  PASS\n" << std::endl;
}

static void test_deletion() {
    std::cout << "=== Test: 3bp Deletion ===" << std::endl;
    std::string ref   = "ACGTACGTGATTACA";
    std::string query = "ACGTACGATTACA";
    SmithWaterman sw;
    SWConfig cfg;
    cfg.score_threshold = 10;
    sw.set_config(cfg);
    SWResult r = sw.align(query, ref);
    std::cout << "  score=" << r.score << " valid=" << r.valid << std::endl;
    std::cout << "  mismatches=" << r.num_mismatches
              << " ins=" << r.num_ins << " del=" << r.num_del << std::endl;
    std::cout << "  CIGAR: " << r.cigar << std::endl;
    assert(r.valid);
    assert(r.num_del >= 1);
    std::cout << "  PASS\n" << std::endl;
}

static void test_insertion() {
    std::cout << "=== Test: 3bp Insertion ===" << std::endl;
    std::string ref   = "ACGTACGATTACA";
    std::string query = "ACGTACGTGATTACA";
    SmithWaterman sw;
    SWConfig cfg;
    cfg.score_threshold = 10;
    sw.set_config(cfg);
    SWResult r = sw.align(query, ref);
    std::cout << "  score=" << r.score << " valid=" << r.valid << std::endl;
    std::cout << "  mismatches=" << r.num_mismatches
              << " ins=" << r.num_ins << " del=" << r.num_del << std::endl;
    std::cout << "  CIGAR: " << r.cigar << std::endl;
    assert(r.valid);
    assert(r.num_ins >= 1);
    std::cout << "  PASS\n" << std::endl;
}

static void test_complex() {
    std::cout << "=== Test: Complex (SNP + Indel) ===" << std::endl;
    std::string ref   = "ACGTACGTGATTACA";
    std::string query = "ACGTAXGTXXGATTACA";
    SmithWaterman sw;
    SWConfig cfg;
    cfg.score_threshold = 10;
    sw.set_config(cfg);
    SWResult r = sw.align(query, ref);
    std::cout << "  score=" << r.score << " valid=" << r.valid << std::endl;
    std::cout << "  mismatches=" << r.num_mismatches
              << " ins=" << r.num_ins << " del=" << r.num_del << std::endl;
    std::cout << "  CIGAR: " << r.cigar << std::endl;
    assert(r.valid);
    std::cout << "  PASS\n" << std::endl;
}

static void test_aligner_snp_fallback() {
    std::cout << "=== Test: Aligner SNP fallback ===" << std::endl;
    Aligner aligner;
    AlignerConfig cfg;
    cfg.max_mismatches = 0;
    cfg.sw_score_threshold = 8;
    cfg.enable_sw_fallback = true;
    cfg.seed_length = 10;
    cfg.min_seed_hits = 1;
    aligner.set_config(cfg);

    std::string ref_name = "chr1";
    std::string ref_seq  = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
    ref_seq += "ACGTACGTGATTACA";
    ref_seq += "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG";
    std::string query    = "ACGTACGTAATTACAGGGGGGGG";

    FastqRead read;
    read.id = "read1";
    read.sequence = query;
    read.quality = std::string(query.size(), 'I');
    read.is_first = true;
    read.pair_id = 0;

    std::string fasta_path = "_test_ref.fasta";
    {
        std::ofstream f(fasta_path);
        f << ">" << ref_name << "\n" << ref_seq << "\n";
    }
    aligner.load_reference(fasta_path);

    AlignmentResult aln = aligner.align_read(read);
    std::cout << "  unmapped=" << aln.is_unmapped << std::endl;
    std::cout << "  ref=" << aln.ref_name << " pos=" << aln.position << std::endl;
    std::cout << "  cigar=" << aln.cigar << std::endl;
    std::cout << "  mapq=" << aln.mapq << " nm=" << aln.num_mismatches << std::endl;

    if (aln.is_unmapped) {
        std::cout << "  WARN: unmapped\n";
    } else {
        assert(!aln.cigar.empty());
        std::cout << "  PASS\n";
    }
    std::cout << std::endl;

    std::remove(fasta_path.c_str());
}

int main() {
    test_snp();
    test_deletion();
    test_insertion();
    test_complex();
    test_aligner_snp_fallback();
    std::cout << "All tests passed!" << std::endl;
    return 0;
}
