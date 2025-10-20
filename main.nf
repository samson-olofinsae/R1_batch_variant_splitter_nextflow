#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Params come from nextflow.config (reads, ref, outdir, max_cpus, ploidy, do_filter, filter_expr)

// -----------------------------
// Pairing preflight (logs) — no conda needed
// -----------------------------
process log_pairing_warnings {
  tag "pairing_check"
  publishDir "${params.outdir}", mode: 'copy', overwrite: true
  conda = null

  input:
    val read_glob

  output:
    path "logs/run_batch.log", emit: batch_log

  cpus 1

  script:
  """
  # drop -u to avoid unbound var on empty shell vars
  set -eo pipefail
  mkdir -p logs

  # Inject the literal glob from Nextflow, convert {R1,R2} -> R1
  r1_glob='${read_glob.replace("{R1,R2}", "R1")}'

  {
    echo "Batch processing started: \$(date -Iseconds)"
    shopt -s nullglob
    found=0
    for r1 in \${r1_glob}; do
      found=1
      base=\$(basename "\$r1")
      id="\${base%_R1.fastq.gz}"
      mate="\${r1/_R1.fastq.gz/_R2.fastq.gz}"
      if [[ -f "\$mate" ]]; then
        echo "Processing \$id..."
      else
        echo "Missing R2 pair for \$id, skipping."
      fi
    done
    if [[ \$found -eq 0 ]]; then
      echo "No FASTQ pairs matched pattern: ${read_glob}"
    fi
  } > logs/run_batch.log
  """
}

// -----------------------------
// Main per-sample calling step (uses conda env)
// -----------------------------
process r1_variant_splitter {
  tag "$sample_id"

  // Publish into the user-requested directory, strip literal 'out/' prefix
  publishDir "${params.outdir}", mode: 'copy', overwrite: true, saveAs: { fname ->
    def s = fname.toString()
    // ignore publishing the top-level 'out' dir itself
    if (s == 'out' || s == './out') return null
    // safely strip the 'out/' prefix if present
    return s.replaceFirst('^out/', '')
  }

  input:
    tuple val(sample_id), path(fq1), path(fq2), path(ref_fa)

  output:
    path "out/**", emit: out_tree
    path "logs/${sample_id}.log", emit: sample_log
    path "logs/${sample_id}.err", emit: sample_err
    path "logs/${sample_id}.status", emit: status
    tuple val(sample_id), path("out/vcf/${sample_id}_final_variants.vcf.gz"), emit: vcf

  cpus { params.max_cpus as int }

  script:
  """
  set -euo pipefail
  mkdir -p logs

  {
    echo "[\$(date -Iseconds)] START sample=${sample_id}"
    split_snv_indel_variants.py \\
      "${fq1}" \\
      "${fq2}" \\
      "${ref_fa}" \\
      --threads ${task.cpus} \\
      --ploidy ${params.ploidy} \\
      --outdir out
    echo "[\$(date -Iseconds)] DONE  sample=${sample_id}"
  } > >(tee logs/${sample_id}.log) 2> >(tee logs/${sample_id}.err >&2)

  echo "Completed ${sample_id} at \$(date -Iseconds)" > "logs/${sample_id}.status"
  """
}

// -----------------------------
// (Optional) normalize & filter (uses conda env)
// -----------------------------
process vcf_norm_filter {
  tag "$sample_id"
  publishDir "${params.outdir}", mode: 'copy', overwrite: true

  input:
    tuple val(sample_id), path(vcf_gz), path(ref_fa)

  output:
    tuple val(sample_id), path("${sample_id}.norm.filt.vcf.gz")

  script:
  """
  set -euo pipefail
  bcftools norm -f "${ref_fa}" -m -any -Oz -o ${sample_id}.norm.vcf.gz "${vcf_gz}"
  bcftools index -f ${sample_id}.norm.vcf.gz
  bcftools filter -i '${params.filter_expr ?: "QUAL>=30 && DP>=5"}' -Oz -o ${sample_id}.norm.filt.vcf.gz ${sample_id}.norm.vcf.gz
  bcftools index -f ${sample_id}.norm.filt.vcf.gz
  """
}

// -----------------------------
// Aggregate counts (tsv table) — no conda needed
// -----------------------------
process summarize_counts {
  tag "mqc_summary"
  publishDir "${params.outdir}/multiqc_cc", mode: 'copy', overwrite: true
  conda = null

  input:
    path mqc_dir

  output:
    path "_summary.tsv"

  script:
  """
  # don't use -u here; it makes 'mqc_dir' fatal if not staged
  set -eo pipefail

  # If the directory isn't there yet (or empty), emit an empty header and exit 0
  if [[ ! -d "${mqc_dir}" ]]; then
    echo -e "Sample\\tsnvs\\tindels" > _summary.tsv
    exit 0
  fi

  {
    echo -e "Sample\\tsnvs\\tindels"
    shopt -s nullglob
    for f in "${mqc_dir}"/*_r1_stats_mqc.tsv; do
      [[ -f "\$f" ]] || continue
      tail -n1 "\$f"
    done | sort -u
  } > _summary.tsv
  """
}

// -----------------------------
// MultiQC HTML report (uses conda env)
// -----------------------------

process multiqc {
  tag "multiqc"
  publishDir "${params.outdir}", mode: 'copy', overwrite: true, saveAs: { f ->
    def s = f.toString()
    return s.startsWith('multiqc_report/') ? s : 'multiqc_report/' + s
  }

  // NEW: reuse tools env instead of a second tiny env
  conda "$projectDir/env_tools.yml"

  input:
    path results_dir

  output:
    path "multiqc_report/**"

  script:
  """
  set -euo pipefail
  mkdir -p multiqc_report
  multiqc -f -o multiqc_report "${results_dir}"
  """
}

// -----------------------------
// Merge batch log (initial + statuses) — no conda needed
// -----------------------------
process batch_log_finalize {
  tag "batch_log"
  publishDir "${params.outdir}", mode: 'copy', overwrite: true
  conda = null

  input:
    path initial_log
    val  done_token

  output:
    path "logs/run_batch.log"

  script:
  """
  set -euo pipefail
  mkdir -p logs
  cp "${initial_log}" logs/run_batch.log
  for s in ${params.outdir}/logs/*.status; do
    [[ -f "\$s" ]] && cat "\$s" >> logs/run_batch.log || true
  done
  echo "Batch finished: \$(date -Iseconds)" >> logs/run_batch.log
  """
}

// -----------------------------
// Workflow wiring
// -----------------------------

workflow {
  log.info "reads=${params.reads}  ref=${params.ref}  outdir=${params.outdir}  cpus=${params.max_cpus}  ploidy=${params.ploidy}"

  // 0) Pairing preflight (shell-only)
  def pairing = log_pairing_warnings(params.reads)
  def batch_log_ch = pairing.batch_log

  // 1) Build R1/R2 pairs
  def pairs_ch = Channel
      .fromFilePairs(params.reads, size: 2, checkIfExists: true)
      .ifEmpty { error "No FASTQ pairs matched: ${params.reads}" }

  // 2) Destructure into (id, R1, R2, ref)
  def tuples_ch = pairs_ch.map { sid, files ->
    tuple( sid as String, files[0], files[1], file(params.ref) )
  }

  // 3) Variant calling per-sample (conda env)
  def call_out = r1_variant_splitter(tuples_ch)

  // 4) Optional normalization / filtering (conda env)
  if (params.do_filter) {
    def vcf_in = call_out.vcf.map { sid, vcf -> tuple(sid, vcf, file(params.ref)) }
    vcf_norm_filter(vcf_in)
  }

  // Gate downstream steps on completion of variant calling
  def done_token = call_out.vcf.collect().map { true }

  // 5) Aggregate counts after samples finished (shell-only)
  def mqc_dir_gate = done_token.map { file("${params.outdir}/multiqc_cc") }
  summarize_counts(mqc_dir_gate)

  // 6) Finalize batch log after all samples finished (shell-only)
  batch_log_finalize(batch_log_ch, done_token)

  // 7) MultiQC HTML after everything is published (conda env)
  def outdir_gate = done_token.map { file(params.outdir) }
  multiqc(outdir_gate)
}
