# R1 Batch Variant Splitter (SNV/INDEL) - Nextflow + Seqera Wave

A compact, teaching-friendly pipeline that aligns paired-end FASTQs, calls variants, and **splits SNVs and INDELs** into dedicated VCFs.
Designed for **reproducibility** (containers via Seqera Wave), **fault tolerance** (skips incomplete pairs), and **auditability** (per-sample logs + run reports + MultiQC).

**Highlights**
- **Batch-safe**: skips samples missing an R2 mate; logs it, keeps running.
- **End-to-end**: BWA-MEM → SAMtools sort/index → BCFtools mpileup/call → split SNVs/INDELs.
- **MultiQC ready**: emits **custom content** tables per sample and a complete **HTML report**.
- **Portable by default**: **no local Conda solve**; images are built/resolved by **Seqera Wave**, executed locally via **Apptainer/Singularity** (or Docker/Podman if preferred).
- **Provenance out-of-the-box**: Nextflow timeline, trace, report, DAG.

---

## Quick start

> **Prereqs**
> - **Java 17–21** (`java -version`)
> - **Nextflow ≥ 25.x** (`nextflow -version`)
> - A container runtime:
>   - **Apptainer/Singularity** (`singularity --version`) **recommended**, or
>   - Docker/Podman (Wave can target these too, see “Runtimes” below)

### Step 0 - Clone the repository

```bash
git clone https://github.com/samson-olofinsae/R1_batch_variant_splitter_nextflow.git
cd R1_batch_variant_splitter_nextflow
```

### Step 1 - Initialise user workspace

Before running on your own data, initialise empty directories for your FASTQs and reference files:

```bash
./scripts/init_user_workspace.sh
# or manually:
# mkdir -p data/user_fastqs ref/user_ref
```

> This script simply creates `data/user_fastqs/` and `ref/user_ref/` if they do not exist,
> drops `.gitkeep` placeholders, and prints short usage hints.
> It ensures that fresh clones run cleanly without missing-folder errors.

### Step 2 - Run the bundled **demo** dataset (fast!)

```bash
nextflow run main.nf -profile wave,demo --ploidy 2 --max_cpus 2
```

**Outputs** land in `results_demo/`:
```
results_demo/
├── bam/                  # *.aligned.sorted.bam (+ .bai)
├── vcf/                  # *_final_variants.vcf.gz (+ .csi/.tbi)
│   ├── snvs/             # *_snvs.vcf.gz (+ index)
│   └── indels/           # *_indels.vcf.gz (+ index)
├── multiqc_cc/           # *_r1_stats_mqc.tsv + _summary.tsv
├── multiqc_report/       # multiqc_report.html + data/
└── logs/                 # run_batch.log + per-sample .log/.err + .status
```

### Step 3 - Run **your data**

Put your files under `data/user_fastqs/` and your reference under `ref/user_ref/` (or point to your own paths), then:

```bash
nextflow run main.nf -profile wave,user --reads 'data/user_fastqs/*_{R1,R2}.fastq.gz' --ref 'ref/user_ref/<your_ref>.fa' --outdir results_user --ploidy 2 --max_cpus 4
```

> You can resume a previous run safely with `-resume`.

---

## The `init_user_workspace.sh` script

Located under `scripts/`, this helper script prepares your workspace for new data.

**Source:**

```bash
#!/usr/bin/env bash
#!/usr/bin/env bash
set -euo pipefail

# Initialise scaffold directories for user data/references in a fresh clone.
# Creates folders and .gitkeep placeholders so Git tracks the structure.

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "$script_dir/.." && pwd)"

mkdir -p "$repo_root/data/user_fastqs" "$repo_root/ref/user_ref"
touch "$repo_root/data/user_fastqs/.gitkeep" "$repo_root/ref/user_ref/.gitkeep"

printf 'Created: %s and %s\n' "data/user_fastqs/" "ref/user_ref/"
printf 'Drop your FASTQs in data/user_fastqs/\n'
printf 'Put your reference FASTA in ref/user_ref/ (indexes optional; auto-indexed if missing)\n'
```

**Purpose:** Ensures reproducible directory structure for users and contributors cloning the repository.

---

## What the pipeline does (per sample)

1. **Align & sort**  
   `bwa mem` → `samtools sort` → `results/*/bam/{sample}.aligned.sorted.bam`
2. **Index BAM**  
   `samtools index`
3. **Pileup & call**  
   `bcftools mpileup` → `bcftools call -m -v` → `*_final_variants.vcf.gz`
4. **Split by type**  
   `bcftools view -v snps` → `vcf/snvs/{sample}_snvs.vcf.gz`  
   `bcftools view -v indels` → `vcf/indels/{sample}_indels.vcf.gz`
5. **Emit MultiQC custom content**  
   `multiqc_cc/{sample}_r1_stats_mqc.tsv` (SNV/INDEL counts)
6. **Aggregate & report**  
   `_summary.tsv` across samples, **MultiQC HTML** in `multiqc_report/`.

> The reference FASTA is auto-indexed if needed (`samtools faidx`, `bwa index`).

---

## Input conventions

- **FASTQs** must be paired as: `NAME_R1.fastq.gz` and `NAME_R2.fastq.gz`.  
  The pipeline discovers pairs via the glob you pass to `--reads`, e.g. `data/.../*_{R1,R2}.fastq.gz`.
- **Reference FASTA**: any readable path, e.g. `ref/hg19_chr8.fa`. A `.fai` and BWA index are created if missing.

---

## Minimal commands you’ll use most

**Demo:**
```bash
nextflow run main.nf -profile wave,demo --ploidy 2 --max_cpus 2
```

**Your data:**
```bash
nextflow run main.nf -profile wave,user \
  --reads 'data/user_fastqs/*_{R1,R2}.fastq.gz' \
  --ref 'ref/user_ref/<your_ref>.fa' \
  --outdir results_user \
  --ploidy 2 --max_cpus 4 -resume
```

---

## Parameters (CLI `--flag value`)

- `--reads` : Glob for paired FASTQs, e.g. `data/trio/*_{R1,R2}.fastq.gz`
- `--ref` : Reference FASTA path
- `--outdir` : Output directory (default `results`)
- `--ploidy` : Ploidy for `bcftools call` (default `2`)
- `--max_cpus` : Threads handed to tools (default `2`)
- `--do_filter` : Toggle optional VCF normalization/filtering (default `false`)
- `--filter_expr` : BCFtools filter expression if `--do_filter` (default `QUAL>=30 && DP>=5`)

---

## Profiles & runtimes

We ship a small set of **Nextflow profiles**:

- `wave` - **Base** runtime. Tells Nextflow to:
  - build/resolve containers from Conda recipes via **Seqera Wave**  
  - execute them locally via **Apptainer/Singularity** (configurable to Docker/Podman)
  - **no local Conda creation** at all.
- `demo` - Provides demo dataset paths (compose as `-profile wave,demo`).
- `user` - No preset paths (compose as `-profile wave,user`), you **must** pass your own `--reads` and `--ref`.

> If you prefer Docker/Podman: switch the Wave engine in `nextflow.config` (e.g., `wave.engine='docker'`) and enable your runtime.

---

## Reproducibility

- **Containers from recipes**: The pipeline points Wave at `env_tools.yml` (bwa, samtools, bcftools, python) and `env_multiqc.yml` (multiqc). Wave **builds** and **pins** those environments into container images, eliminating local Conda solver variability.
Containers frozen with wave.freezeContainers=true for bitwise reproducibility.
- **Provenance artefacts** (written automatically at the repo root):
  - `pipeline_timeline.html` - interactive run timeline
  - `pipeline_trace.txt` - task-by-task resource stats
  - `pipeline_report.html` - aggregated run metadata
  - `pipeline_dag.svg` - process dependency graph

For stricter locking, you can enable in `nextflow.config`:
```groovy
// inside the 'wave' profile
// wave.freezeContainers = true    // pin built images by digest
// wave.strictPull       = true    // fail if an image can’t be resolved
```

---

## Repository layout

```
project_root/
├── main.nf
├── nextflow.config
├── env_tools.yml
├── env_multiqc.yml
├── scripts/
│   └── split_snv_indel_variants.py
├── data/
│   ├── demo_fastqs/   # bundled tiny paired-end FASTQs
│   └── user_fastqs/   # put your own FASTQs here (optional)
├── ref/
│   ├── demo_ref/      # bundled small reference (demo.fa)
│   └── user_ref/      # put your own reference here (optional)
└── results_*/         # created at runtime
```

---

## Interpreting outputs

- **Per-sample BAM**: `bam/{sample}.aligned.sorted.bam` (+ `.bai`)
- **Per-sample VCF**: `vcf/{sample}_final_variants.vcf.gz` (+ index)
- **Split VCFs**: `vcf/snvs/{sample}_snvs.vcf.gz`, `vcf/indels/{sample}_indels.vcf.gz`
- **MultiQC custom content**: `multiqc_cc/{sample}_r1_stats_mqc.tsv` and combined `_summary.tsv`
- **MultiQC HTML**: `multiqc_report/multiqc_report.html`
- **Logs**:  
  - Global: `logs/run_batch.log` (start/finish + “Processing…/Missing R2…/Completed …”)  
  - Per-sample: `logs/{sample}.log`, `logs/{sample}.err`, and `{sample}.status`

---

## Tuning & tips

- **Threads**: `--max_cpus N` controls `bwa -t`, `samtools sort -@`, and pipeline parallelism.
- **Ploidy**: change `--ploidy` for organism/chromosome context (e.g., `1` for haploid).
- **Optional filtering**:  
  Enable with `--do_filter --filter_expr 'QUAL>=30 && DP>=10'` to produce `*.norm.filt.vcf.gz`.
- **Resume**: safe to use `-resume` to continue partial runs.
- **Disk cleanup**: remove `work/` and prior `results_*` to reclaim space after validation.

---

## Troubleshooting

- **“Missing R2 pair … skipping.”**  
  Check file naming: `*_R1.fastq.gz` and `*_R2.fastq.gz` must both exist for each sample.
- **Reference index errors**  
  The script auto-runs `samtools faidx` and `bwa index` if needed; ensure the reference is writable (or pre-index it in place).
- **Singularity/Apptainer not found**  
  Install Apptainer/Singularity or switch Wave to `docker`/`podman` and enable that runtime.
- **/usr/bin/env: ‘python3\r’: No such file or directory**
  If you see /usr/bin/env: ‘python3\r’: No such file or directory, your checkout has Windows line endings. Fix with git config core.autocrlf input and re-clone, or run sed -i 's/\r$//' scripts/*.py.
---

## Citation

If you use this in a talk, class, or application dossier, please cite:

- **Nextflow** - Di Tommaso *et al.*, *Nat Biotechnol* (2017), doi:10.1038/nbt.3820  
- **MultiQC** - Ewels *et al.*, *Bioinformatics* (2016), btw354


---

## License

MIT (c) 2025 Samson Olofinsae. See `LICENSE`.
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](./LICENSE)

---

## Acknowledgements

- **Seqera** for Wave and containers guidance.
- Phil for steering us toward Seqera Containers and training best practices.
- Everyone testing the demo on WSL2/macOS/Linux - your feedback shaped the ergonomics.

---

### Maintainer notes (optional to keep at the bottom)

- Profiles shipped: `wave` (base), `demo`, `user`.  
  We intentionally **don’t** prefill `user` paths to keep it portable to any dataset layout.  
- The demo FASTQs include a mix of SNVs/INDELs so MultiQC tables aren’t degenerate.  
- For WSL2 users, consider a `.wslconfig` with at least 8 GB RAM and 8 GB swap for smoother local caching.


