#!/usr/bin/env python3
import os
import random
import sys

def generate_reference_genome(output_path, num_chromosomes=3, chrom_length=50000):
    with open(output_path, 'w') as f:
        for i in range(num_chromosomes):
            f.write(f">chr{i+1}\n")
            seq = []
            for _ in range(chrom_length):
                seq.append(random.choice(['A', 'C', 'G', 'T']))
            seq_str = ''.join(seq)
            for j in range(0, len(seq_str), 80):
                f.write(seq_str[j:j+80] + '\n')
    print(f"Generated reference genome: {output_path}")

def generate_fastq_reads(ref_path, output_path, num_reads=5000, read_length=100):
    ref_seqs = {}
    current_chrom = None
    current_seq = []

    with open(ref_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_chrom:
                    ref_seqs[current_chrom] = ''.join(current_seq)
                current_chrom = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        if current_chrom:
            ref_seqs[current_chrom] = ''.join(current_seq)

    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}

    def reverse_complement(s):
        return ''.join(complement[c] for c in reversed(s))

    with open(output_path, 'w') as f:
        for i in range(num_reads):
            chrom = random.choice(list(ref_seqs.keys()))
            seq = ref_seqs[chrom]
            pos = random.randint(0, len(seq) - read_length)
            read_seq = seq[pos:pos+read_length]

            for _ in range(random.randint(0, 2)):
                mut_pos = random.randint(0, read_length - 1)
                bases = ['A', 'C', 'G', 'T']
                bases.remove(read_seq[mut_pos])
                read_seq = read_seq[:mut_pos] + random.choice(bases) + read_seq[mut_pos+1:]

            is_reverse = random.random() < 0.5
            if is_reverse:
                read_seq = reverse_complement(read_seq)

            qual = ''.join(chr(random.randint(30, 40)) for _ in range(read_length))

            f.write(f"@read_{i+1} {chrom}:{pos+1}\n")
            f.write(read_seq + '\n')
            f.write("+\n")
            f.write(qual + '\n')

    print(f"Generated FASTQ reads: {output_path} ({num_reads} reads)")

def generate_paired_fastq(ref_path, out1, out2, num_pairs=5000, read_length=100, insert_size=350):
    ref_seqs = {}
    current_chrom = None
    current_seq = []

    with open(ref_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_chrom:
                    ref_seqs[current_chrom] = ''.join(current_seq)
                current_chrom = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        if current_chrom:
            ref_seqs[current_chrom] = ''.join(current_seq)

    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}

    def reverse_complement(s):
        return ''.join(complement[c] for c in reversed(s))

    with open(out1, 'w') as f1, open(out2, 'w') as f2:
        for i in range(num_pairs):
            chrom = random.choice(list(ref_seqs.keys()))
            seq = ref_seqs[chrom]
            start = random.randint(0, len(seq) - insert_size)
            template = seq[start:start+insert_size]

            r1_seq = template[:read_length]
            r2_seq = reverse_complement(template[-read_length:])

            qual1 = ''.join(chr(random.randint(30, 40)) for _ in range(read_length))
            qual2 = ''.join(chr(random.randint(30, 40)) for _ in range(read_length))

            f1.write(f"@read_{i+1}/1 {chrom}:{start+1}\n")
            f1.write(r1_seq + '\n')
            f1.write("+\n")
            f1.write(qual1 + '\n')

            f2.write(f"@read_{i+1}/2 {chrom}:{start+1}\n")
            f2.write(r2_seq + '\n')
            f2.write("+\n")
            f2.write(qual2 + '\n')

    print(f"Generated paired-end FASTQ: {out1}, {out2} ({num_pairs} pairs)")

if __name__ == "__main__":
    random.seed(42)
    os.makedirs("test_data", exist_ok=True)

    ref_path = os.path.join("test_data", "reference.fasta")
    fq_path = os.path.join("test_data", "reads.fastq")
    fq1_path = os.path.join("test_data", "reads_R1.fastq")
    fq2_path = os.path.join("test_data", "reads_R2.fastq")

    generate_reference_genome(ref_path)
    generate_fastq_reads(ref_path, fq_path)
    generate_paired_fastq(ref_path, fq1_path, fq2_path)
    print("Done!")
