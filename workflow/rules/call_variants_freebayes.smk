#
# Equivalence check -- stock freebayes on the SAME extracted region + BAM
# that call_variants (panfreebayes) already used.
#
# This is NOT part of the normal pipeline/`rule all` -- it exists purely to
# answer one question: does panfreebayes agree with stock freebayes on
# IDENTICAL input? This is the same question test/panfreebayes/
# acceptance_check.sh's "equivalence" mode answers, just reached through the
# Snakemake workflow's own already-produced region FASTA/BAM instead of a
# separately-run acceptance_check.sh invocation -- deliberately avoids
# re-running extract_region/align_region a second time (both rules' outputs
# are reused as-is; Snakemake will skip re-running them since they already
# exist and are up to date).
#
# Why this matters right now: comparing a fresh pipeline run against the
# committed baselines (test/panfreebayes/baselines/) conflates two different
# questions -- "does panfreebayes agree with freebayes" and "does our
# extraction/padding match whatever produced the baseline's reference FASTA
# (unknown -- see panfreebayes_milestone6_snakemake_progress.md)". Running
# stock freebayes against the exact same FASTA/BAM panfreebayes just used
# isolates the first question cleanly, with no baseline-provenance ambiguity
# at all.
#
# freebayes -f <fasta> <calling_flags...> <bam>  -- matches buildArgv()'s
# own construction order in src/panfreebayes/panfreebayes_core.cpp (-f first,
# flags, bam last) exactly.
#
rule call_variants_freebayes:
    input:
        ref=f"{RESULTS}/regions/{{bubble}}.fasta",
        bam=f"{RESULTS}/alignments/{{strain}}/{{bubble}}.bam",
        bai=f"{RESULTS}/alignments/{{strain}}/{{bubble}}.bam.bai",
    output:
        vcf=f"{RESULTS}/calls_freebayes/{{strain}}/{{bubble}}.vcf",
    params:
        fb=config["freebayes_bin"],
        flags=" ".join(config["calling_flags"]),
    log:
        f"{RESULTS}/logs/call_variants_freebayes/{{strain}}/{{bubble}}.log",
    shell:
        r"""
        {params.fb} -f {input.ref} {params.flags} {input.bam} \
            > {output.vcf} 2> {log}
        """
