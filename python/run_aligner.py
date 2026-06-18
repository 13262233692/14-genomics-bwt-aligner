#!/usr/bin/env python3
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))

import bwt_aligner

def main():
    import argparse
    parser = argparse.ArgumentParser(description="BWT-based genomic sequence aligner")
    parser.add_argument("reference", help="Reference genome in FASTA format")
    parser.add_argument("-1", "--read1", help="FASTQ file with read 1 (single-end or paired R1)")
    parser.add_argument("-2", "--read2", help="FASTQ file with read 2 (paired R2)")
    parser.add_argument("-o", "--output", required=True, help="Output SAM file")
    parser.add_argument("-t", "--threads", type=int, default=4, help="Number of threads")
    parser.add_argument("--seed-length", type=int, default=16, help="Seed length")
    parser.add_argument("--max-mismatches", type=int, default=2, help="Max mismatches allowed")
    parser.add_argument("--batch-size", type=int, default=10000, help="Batch size per thread")

    args = parser.parse_args()

    if not args.read1:
        parser.error("At least --read1 is required")

    config = bwt_aligner.AlignerConfig()
    config.num_threads = args.threads
    config.seed_length = args.seed_length
    config.max_mismatches = args.max_mismatches
    config.batch_size = args.batch_size

    aligner = bwt_aligner.Aligner()
    aligner.set_config(config)

    print(f"Building index from {args.reference}...")
    if not aligner.build_index(args.reference):
        print("ERROR: Failed to build index")
        sys.exit(1)

    refs = aligner.get_references()
    total_bp = sum(r.length for r in refs)
    print(f"Loaded {len(refs)} sequences, {total_bp} total bp")

    if args.read2:
        print(f"Aligning paired-end reads: {args.read1} + {args.read2}")
        print(f"Using {args.threads} threads...")
        if not aligner.align_paired(args.read1, args.read2, args.output):
            print("ERROR: Paired-end alignment failed")
            sys.exit(1)
    else:
        print(f"Aligning single-end reads: {args.read1}")
        print(f"Using {args.threads} threads...")
        if not aligner.align_single(args.read1, args.output):
            print("ERROR: Single-end alignment failed")
            sys.exit(1)

    print(f"\nAlignment complete! Output: {args.output}")

if __name__ == "__main__":
    main()
