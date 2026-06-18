import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import bwt_aligner

def test_fm_index_debug():
    print("=" * 60)
    print("FM-index Debug Test")
    print("=" * 60)

    text = "ACGTACGTGATTACA"
    print(f"Text: {text}")
    print(f"Text length: {len(text)}")

    fm = bwt_aligner.FMIndex()
    fm.build(text)
    print(f"FM-index size: {fm.size()}")
    print(f"Text stored: {fm.get_text()}")

    sa = fm.get_sa()
    print(f"\nSuffix Array (size {sa.size()}):")
    for i in range(min(sa.size(), 20)):
        sa_val = sa[i]
        suffix = fm.get_text()[sa_val:] if sa_val < len(fm.get_text()) else "$"
        print(f"  SA[{i:2d}] = {sa_val:3d} -> {suffix}")

    bwt_obj = fm.get_bwt()
    bwt_str = bwt_obj.get_bwt()
    print(f"\nBWT string: '{bwt_str}'")
    print(f"BWT primary pos: {bwt_obj.get_primary_pos()}")

    print(f"\nC table:")
    for c in ['$', '#', 'A', 'C', 'G', 'T']:
        print(f"  C['{c}'] = {fm.get_c(c)}")

    print(f"\nOcc table test:")
    for c in ['A', 'C', 'G', 'T', '$']:
        print(f"  Occ('{c}', 0..{len(bwt_str)-1}): ", end="")
        occs = [fm.get_occ(i, c) for i in range(len(bwt_str))]
        print(occs)

    print(f"\nSA samples:")
    sa_size = sa.size()
    sample_rate = 32
    for i in range(0, sa_size, sample_rate):
        print(f"  SA[{i}] = {sa[i]} (sample)")

    print(f"\nSearch test:")
    search = bwt_aligner.BackwardSearch()
    search.set_fm_index(fm)

    patterns = ["ACGT", "TGAT", "TACA", "ATT"]
    for p in patterns:
        r = search.search(p)
        print(f"\nPattern '{p}': SearchRange({r.low}, {r.high}), count={r.count()}")
        if r.valid():
            for i in range(r.low, r.high + 1):
                pos = fm.resolve_position(i)
                sa_val = fm.get_sa_value(i)
                print(f"  row {i}: SA_val={sa_val}, resolved={pos}")

    print(f"\nResolve position test (all rows):")
    bad_count = 0
    for i in range(fm.size()):
        sa_val = fm.get_sa_value(i)
        resolved = fm.resolve_position(i)
        if resolved != sa_val:
            bad_count += 1
            if bad_count <= 5:
                print(f"  MISMATCH row {i}: SA_val={sa_val}, resolved={resolved}")
    if bad_count == 0:
        print("  All positions resolved correctly!")
    else:
        print(f"  Total mismatches: {bad_count}/{fm.size()}")


def count_sam_mapped(sam_path):
    print("\n" + "=" * 60)
    print(f"SAM Analysis: {sam_path}")
    print("=" * 60)

    if not os.path.exists(sam_path):
        print("SAM file not found")
        return

    total = 0
    mapped = 0
    unmapped = 0

    with open(sam_path, 'r') as f:
        for line in f:
            if line.startswith('@'):
                continue
            total += 1
            parts = line.split('\t')
            if len(parts) >= 5:
                flag = int(parts[1])
                pos = int(parts[3])
                if flag & 4:
                    unmapped += 1
                else:
                    mapped += 1
                    if mapped <= 10:
                        print(f"  Mapped read {parts[0]}: chr={parts[2]}, pos={pos}, flag={flag}, mapq={parts[4]}")

    print(f"\nSummary: Total={total}, Mapped={mapped}, Unmapped={unmapped}")
    if total > 0:
        print(f"Mapped rate: {100.0*mapped/total:.2f}%")


if __name__ == "__main__":
    test_fm_index_debug()

    sam_dir = os.path.join(os.path.dirname(__file__), "test_data")
    if os.path.exists(sam_dir):
        count_sam_mapped(os.path.join(sam_dir, "output_single.sam"))
