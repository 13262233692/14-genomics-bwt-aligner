#!/usr/bin/env python3
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))

try:
    import bwt_aligner
except ImportError:
    print("Error: bwt_aligner module not found. Please build the C++ extension first.")
    print("Run: mkdir build && cd build && cmake .. && cmake --build . --config Release")
    sys.exit(1)

def test_single_end_alignment():
    print("=" * 60)
    print("Test 1: Single-end Alignment")
    print("=" * 60)

    config = bwt_aligner.AlignerConfig()
    config.num_threads = 4
    config.seed_length = 16
    config.max_mismatches = 2
    config.batch_size = 1000

    aligner = bwt_aligner.Aligner()
    aligner.set_config(config)

    ref_path = os.path.join(os.path.dirname(__file__), "test_data", "reference.fasta")
    fq_path = os.path.join(os.path.dirname(__file__), "test_data", "reads.fastq")
    sam_path = os.path.join(os.path.dirname(__file__), "test_data", "output_single.sam")

    if not os.path.exists(ref_path) or not os.path.exists(fq_path):
        print("Test data not found. Run generate_test_data.py first.")
        return False

    print("Building FM-index...")
    if not aligner.build_index(ref_path):
        print("Failed to build index!")
        return False

    refs = aligner.get_references()
    print(f"Loaded {len(refs)} reference sequences")
    for r in refs:
        print(f"  {r.name}: {r.length} bp")

    print("\nRunning single-end alignment...")
    if not aligner.align_single(fq_path, sam_path):
        print("Alignment failed!")
        return False

    print(f"SAM file written to: {sam_path}")
    print("Single-end alignment PASSED\n")
    return True

def test_paired_end_alignment():
    print("=" * 60)
    print("Test 2: Paired-end Alignment")
    print("=" * 60)

    config = bwt_aligner.AlignerConfig()
    config.num_threads = 4
    config.seed_length = 16
    config.max_mismatches = 2
    config.batch_size = 1000

    aligner = bwt_aligner.Aligner()
    aligner.set_config(config)

    ref_path = os.path.join(os.path.dirname(__file__), "test_data", "reference.fasta")
    fq1_path = os.path.join(os.path.dirname(__file__), "test_data", "reads_R1.fastq")
    fq2_path = os.path.join(os.path.dirname(__file__), "test_data", "reads_R2.fastq")
    sam_path = os.path.join(os.path.dirname(__file__), "test_data", "output_paired.sam")

    if not os.path.exists(ref_path) or not os.path.exists(fq1_path):
        print("Test data not found. Run generate_test_data.py first.")
        return False

    print("Building FM-index...")
    if not aligner.build_index(ref_path):
        print("Failed to build index!")
        return False

    print("Running paired-end alignment...")
    if not aligner.align_paired(fq1_path, fq2_path, sam_path):
        print("Alignment failed!")
        return False

    print(f"SAM file written to: {sam_path}")
    print("Paired-end alignment PASSED\n")
    return True

def test_read_level_alignment():
    print("=" * 60)
    print("Test 3: Low-level FM-index and Search APIs")
    print("=" * 60)

    text = "ACGTACGTGATTACA"
    print(f"Test sequence: {text}")

    fm = bwt_aligner.FMIndex()
    fm.build(text)
    print(f"FM-index built, size: {fm.size()}")

    search = bwt_aligner.BackwardSearch()
    search.set_fm_index(fm)

    patterns = ["ACGT", "TGAT", "TACA", "AAAA"]
    for p in patterns:
        r = search.search(p)
        if r.valid():
            positions = []
            for i in range(r.low, r.high + 1):
                pos = fm.resolve_position(i)
                positions.append(pos)
            print(f"  Pattern '{p}': found at positions {positions}")
        else:
            print(f"  Pattern '{p}': not found")

    print("Low-level API test PASSED\n")
    return True

def test_sequence_utils():
    print("=" * 60)
    print("Test 4: Sequence Utilities")
    print("=" * 60)

    seq = "ACGTAGCT"
    rc = bwt_aligner.reverse_complement(seq)
    print(f"Sequence: {seq}")
    print(f"Reverse complement: {rc}")
    assert rc == "AGCTACGT", f"Expected AGCTACGT, got {rc}"

    for c in "ACGT":
        idx = bwt_aligner.base_to_index(c)
        back = bwt_aligner.index_to_base(idx)
        print(f"  {c} -> idx {idx} -> {back}")
        assert c == back

    print("Sequence utilities PASSED\n")
    return True

def test_fastq_reader():
    print("=" * 60)
    print("Test 5: FASTQ Reader")
    print("=" * 60)

    fq_path = os.path.join(os.path.dirname(__file__), "test_data", "reads.fastq")
    if not os.path.exists(fq_path):
        print("Test data not found. Skipping.")
        return True

    reader = bwt_aligner.FastqReader()
    if not reader.open(fq_path, True):
        print("Failed to open FASTQ!")
        return False

    reads = []
    reader.read_batch(reads, 5)
    print(f"Read {len(reads)} reads:")
    for r in reads[:3]:
        print(f"  {r.id}: {r.sequence[:30]}... (Q={r.quality[:30]}...)")

    print("FASTQ reader PASSED\n")
    reader.close()
    return True

if __name__ == "__main__":
    print("\n" + "=" * 60)
    print("BWT Aligner Python API Test Suite")
    print("=" * 60 + "\n")

    tests = [
        ("Sequence Utilities", test_sequence_utils),
        ("Low-level API", test_read_level_alignment),
        ("FASTQ Reader", test_fastq_reader),
        ("Single-end Alignment", test_single_end_alignment),
        ("Paired-end Alignment", test_paired_end_alignment),
    ]

    passed = 0
    failed = 0

    for name, test_func in tests:
        try:
            if test_func():
                passed += 1
            else:
                failed += 1
                print(f"FAILED: {name}\n")
        except Exception as e:
            failed += 1
            print(f"ERROR in {name}: {e}\n")
            import traceback
            traceback.print_exc()

    print("=" * 60)
    print(f"Results: {passed} passed, {failed} failed")
    print("=" * 60)

    sys.exit(0 if failed == 0 else 1)
