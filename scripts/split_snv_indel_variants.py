#!/usr/bin/env python3
import sys, os, re, subprocess
from pathlib import Path

def run(argv):
    """Run a command passed as a list of args."""
    subprocess.run(argv, check=True, text=True)

def run_bash(pipeline):
    """Run a shell pipeline with bash + pipefail."""
    subprocess.run(['bash', '-o', 'pipefail', '-c', pipeline], check=True, text=True)

def check_or_build_indexes(ref):
    ref = Path(ref)
    # samtools faidx
    if not Path(str(ref) + ".fai").exists():
        print(f"[i] Indexing reference with samtools faidx: {ref}")
        run(["samtools", "faidx", str(ref)])
    # BWA index (check for one file)
    if not Path(str(ref) + ".bwt").exists():
        print(f"[i] Indexing reference with bwa index: {ref}")
        run(["bwa", "index", str(ref)])

def infer_sample_name(fq1):
    name = Path(fq1).name
    patterns = [
        r"_R1\.f(ast)?q\.gz$", r"_1\.f(ast)?q\.gz$", r"\.R1\.f(ast)?q\.gz$",
        r"-R1\.f(ast)?q\.gz$", r"\.1\.f(ast)?q\.gz$"
    ]
    for pat in patterns:
        if re.search(pat, name):
            return re.sub(pat, "", name)
    return re.sub(r"\.f(ast)?q\.gz$", "", name)

def count_records(vcf_path):
    """Count non-header lines (#...) in .vcf or .vcf.gz without shell."""
    p = Path(vcf_path)
    if not p.exists():
        return 0
    try:
        if p.suffix == ".gz":
            import gzip
            with gzip.open(p, "rt", encoding="utf-8", errors="ignore") as fh:
                return sum(1 for line in fh if line and line[0] != '#')
        else:
            with open(p, "r", encoding="utf-8", errors="ignore") as fh:
                return sum(1 for line in fh if line and line[0] != '#')
    except Exception:
        return 0

def main():
    import argparse
    ap = argparse.ArgumentParser(
        description="Align, call variants, split SNVs/INDELs, emit MultiQC custom content (portable)."
    )
    ap.add_argument("fq1")
    ap.add_argument("fq2")
    ap.add_argument("ref")
    ap.add_argument("--ploidy", default="2", help="Ploidy for bcftools call (default: 2)")
    ap.add_argument("--threads", type=int, default=4, help="CPUs for bwa/samtools/bcftools (threads used where supported)")
    ap.add_argument("--outdir", default="results", help="Output directory (default: results)")
    args = ap.parse_args()

    fq1, fq2, ref = args.fq1, args.fq2, args.ref
    if not Path(fq1).exists() or not Path(fq2).exists():
        sys.exit("Input FASTQs not found.")
    if not Path(ref).exists():
        sys.exit("Reference FASTA not found.")

    sample = infer_sample_name(fq1)
    wd = Path(os.getcwd())
    results_dir = wd / args.outdir
    bam_dir = results_dir / "bam"
    vcf_dir = results_dir / "vcf"
    snv_dir = vcf_dir / "snvs"
    indel_dir = vcf_dir / "indels"
    for d in (bam_dir, vcf_dir, snv_dir, indel_dir):
        d.mkdir(parents=True, exist_ok=True)

    sorted_bam = bam_dir / f"{sample}.aligned.sorted.bam"
    final_vcf_gz = vcf_dir / f"{sample}_final_variants.vcf.gz"
    snv_gz = snv_dir / f"{sample}_snvs.vcf.gz"
    indel_gz = indel_dir / f"{sample}_indels.vcf.gz"

    check_or_build_indexes(ref)

    print(f"[+] Starting variant pipeline for: {sample}")

    # 1) Align + sort (bash + pipefail). Threads supported by bwa and samtools.
    cmd_align = (
        f"bwa mem -t {args.threads} '{ref}' '{fq1}' '{fq2}' "
        f"| samtools sort -@ {args.threads} -o '{sorted_bam}' -"
    )
    print(f"[run] {cmd_align}")
    run_bash(cmd_align)

    # 2) Index BAM (threads may be ignored by older samtools index; harmless if so)
    print(f"[run] samtools index -@ {args.threads} {sorted_bam}")
    try:
        run(["samtools", "index", "-@", str(args.threads), str(sorted_bam)])
    except subprocess.CalledProcessError:
        # retry without threads flag for compatibility
        print("[warn] 'samtools index -@' not supported; retrying without threads.")
        run(["samtools", "index", str(sorted_bam)])

    # 3) mpileup -> 4) call (portable: do NOT pass '-@' which older bcftools doesn't support)
    cmd_call = (
        f"bcftools mpileup -Ou -f '{ref}' -a FORMAT/DP -q 10 -Q 13 '{sorted_bam}' "
        f"| bcftools call -mv --ploidy {args.ploidy} -Oz -o '{final_vcf_gz}'"
    )
    print(f"[run] {cmd_call}")
    run_bash(cmd_call)

    # 5) Index VCF
    print(f"[run] bcftools index -f {final_vcf_gz}")
    try:
        run(["bcftools", "index", "-f", str(final_vcf_gz)])
    except subprocess.CalledProcessError:
        print("[warn] 'bcftools index -f' not supported; retrying without -f.")
        run(["bcftools", "index", str(final_vcf_gz)])

    # 6) Split SNVs / INDELs
    print(f"[run] bcftools view -v snps -Oz -o {snv_gz} {final_vcf_gz}")
    run(["bcftools", "view", "-v", "snps", "-Oz", "-o", str(snv_gz), str(final_vcf_gz)])
    print(f"[run] bcftools index -f {snv_gz}")
    try:
        run(["bcftools", "index", "-f", str(snv_gz)])
    except subprocess.CalledProcessError:
        print("[warn] 'bcftools index -f' not supported; retrying without -f.")
        run(["bcftools", "index", str(snv_gz)])

    print(f"[run] bcftools view -v indels -Oz -o {indel_gz} {final_vcf_gz}")
    run(["bcftools", "view", "-v", "indels", "-Oz", "-o", str(indel_gz), str(final_vcf_gz)])
    print(f"[run] bcftools index -f {indel_gz}")
    try:
        run(["bcftools", "index", "-f", str(indel_gz)])
    except subprocess.CalledProcessError:
        print("[warn] 'bcftools index -f' not supported; retrying without -f.")
        run(["bcftools", "index", str(indel_gz)])

    # 7) MultiQC custom content
    mqc_dir = results_dir / "multiqc_cc"
    mqc_dir.mkdir(parents=True, exist_ok=True)
    mqc_tsv = mqc_dir / f"{sample}_r1_stats_mqc.tsv"
    snv_count = count_records(snv_gz)
    indel_count = count_records(indel_gz)

    with open(mqc_tsv, "w", encoding="utf-8") as fh:
        fh.write("# id: r1_variant_splitter\n")
        fh.write("# section_name: R1 Variant Splitter — variant summary\n")
        fh.write("# description: SNV/INDEL counts per sample (from bcftools outputs)\n")
        fh.write("# plot_type: table\n")
        fh.write("# file_format: tsv\n")
        fh.write("Sample\tsnvs\tindels\n")
        fh.write(f"{sample}\t{snv_count}\t{indel_count}\n")

    print(f"[✓] Pipeline completed successfully for: {sample}")
    print(f"[i] Outputs:\n  BAM: {sorted_bam}\n  VCF: {final_vcf_gz}\n  SNVs: {snv_gz}\n  INDELs: {indel_gz}\n  MultiQC CC: {mqc_tsv}")

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: split_snv_indel_variants.py <R1.fastq.gz> <R2.fastq.gz> <ref_genome.fa> [--ploidy N] [--threads N] [--outdir PATH]")
        sys.exit(1)
    main()
