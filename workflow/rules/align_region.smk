#
# Stage 3 -- alignment (per strain x bubble).
#
# Verbatim transplant of scripts/panfreebayes/align_region.sh's own commands:
#
#   minimap2 -ax map-pb -t <threads> <ref.fasta> <fastq> | samtools sort -o <bam>
#   samtools index <bam>
#
# No flags beyond what's shown are added, matching align_region.sh exactly
# (confirmed by the supervisor that none else were used -- in particular,
# `samtools sort` gets no -@).
#
# `threads:` is sourced from config["minimap2_threads"] and Snakemake's own
# --cores-aware scheduling, rather than align_region.sh's `nproc`-or-4
# fallback -- this workflow's own choice (Snakemake already owns the thread
# budget across concurrently-running rule instances, so re-deriving nproc per
# rule invocation would just fight that).
#
# shell.prefix("set -euo pipefail; ") in the Snakefile makes a minimap2
# failure fail this rule even though it's the left side of a pipe into
# samtools sort -- the same guarantee align_region.sh's own `set -o
# pipefail` provides.
#
rule align_region:
    input:
        ref=f"{RESULTS}/regions/{{bubble}}.fasta",
        fastq=lambda wc: config["strains"][wc.strain],
    output:
        bam=f"{RESULTS}/alignments/{{strain}}/{{bubble}}.bam",
        bai=f"{RESULTS}/alignments/{{strain}}/{{bubble}}.bam.bai",
    threads: config["minimap2_threads"]
    log:
        f"{RESULTS}/logs/align_region/{{strain}}/{{bubble}}.log",
    shell:
        r"""
        minimap2 -ax map-pb -t {threads} {input.ref} {input.fastq} 2> {log} \
            | samtools sort -o {output.bam} 2>> {log}
        samtools index {output.bam} 2>> {log}
        """
