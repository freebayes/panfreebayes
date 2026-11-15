#
# Stage 3 -- alignment (per strain x bubble).
#
#   minimap2 -ax <preset> -t <threads> -R <read group> <ref.fasta> <fastq> [<fastq2>]
#       | samtools sort -@ <threads> -o <bam>
#   samtools index <bam>
#
# The original manual pipeline used no -R and no `samtools sort -@`; both
# were added for multi-sample/joint calling (haplotype_mode):
#   - -R '@RG\tID:<strain>\tSM:<strain>[\tPL:ILLUMINA]' is always set (PL only
#     for preset sr), so each sample keeps its own name through a merged BAM
#     and gets its own VCF column. See common.smk's read_group().
#   - a strains entry may be {r1: ..., r2: ...}; both files go to minimap2
#     as a pair. See common.smk's strain_reads().
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
        fastq=lambda wc: strain_reads(wc.strain),
    output:
        bam=f"{RESULTS}/alignments/{{strain}}/{{bubble}}.bam",
        bai=f"{RESULTS}/alignments/{{strain}}/{{bubble}}.bam.bai",
    params:
        preset=config.get("minimap2_preset", "map-pb"),
        rg=lambda wc: read_group(wc.strain),
    threads: config["minimap2_threads"]
    log:
        f"{RESULTS}/logs/align_region/{{strain}}/{{bubble}}.log",
    shell:
        r"""
        minimap2 -ax {params.preset} -t {threads} -R '{params.rg}' \
            {input.ref} {input.fastq} 2> {log} \
            | samtools sort -@ {threads} -o {output.bam} 2>> {log}
        samtools index {output.bam} 2>> {log}
        """
