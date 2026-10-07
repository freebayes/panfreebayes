#
# Stage 3 -- alignment (per strain x bubble).
#
#   minimap2 -ax <preset> -t <threads> <ref.fasta> <fastq> | samtools sort -o <bam>
#   samtools index <bam>
#
# No flags beyond what's shown are added (confirmed by the supervisor that
# none else were used -- in particular, `samtools sort` gets no -@).
#
# <preset> is config["minimap2_preset"] (default "map-pb", the PacBio CLR
# preset this workflow was first validated against) -- NOT hardcoded, so a
# different sequencing technology (map-hifi, map-ont, sr, ...) is a config
# change, not a rule edit. This is the one place the workflow used to assume
# a specific read technology; see panfreebayes_milestone6_snakemake_progress.md
# for why it was pulled out.
#
# `threads:` is sourced from config["minimap2_threads"] and Snakemake's own
# --cores-aware scheduling, rather than re-deriving `nproc` per invocation --
# this workflow's own choice (Snakemake already owns the thread budget across
# concurrently-running rule instances, so a per-rule nproc call would just
# fight that).
#
# shell.prefix("set -euo pipefail; ") in the Snakefile makes a minimap2
# failure fail this rule even though it's the left side of a pipe into
# samtools sort -- the same guarantee align_region.sh's own (now-removed)
# `set -o pipefail` provided.
#
rule align_region:
    input:
        ref=f"{RESULTS}/regions/{{bubble}}.fasta",
        fastq=lambda wc: config["strains"][wc.strain],
    output:
        bam=f"{RESULTS}/alignments/{{strain}}/{{bubble}}.bam",
        bai=f"{RESULTS}/alignments/{{strain}}/{{bubble}}.bam.bai",
    params:
        preset=config.get("minimap2_preset", "map-pb"),
    threads: config["minimap2_threads"]
    log:
        f"{RESULTS}/logs/align_region/{{strain}}/{{bubble}}.log",
    shell:
        r"""
        minimap2 -ax {params.preset} -t {threads} {input.ref} {input.fastq} 2> {log} \
            | samtools sort -o {output.bam} 2>> {log}
        samtools index {output.bam} 2>> {log}
        """
