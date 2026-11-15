#
# Stage 4 -- variant calling (per strain x bubble).
#
# Invokes the compiled panfreebayes CLI exactly as the manual pipeline did:
#
#   <panfreebayes_bin> call --ref <fasta> --bam <bam> -- <calling_flags...>
#
# including the "--" separator before forwarded calling flags. The default
# calling_flags in config/config.yaml are the brief's own baseline
# combination (--pooled-continuous --min-alternate-count 2
# --min-alternate-fraction 0.2 --limit-coverage 200), the exact flags
# test/panfreebayes/baselines/ were produced with.
#
# panfreebayes_bin is assumed pre-built (config-provided path) -- this
# workflow does not add a rule that builds it via ninja; see
# panfreebayes_milestone6_snakemake_progress.md for why (this workflow's own
# choice, not instructed).
#
# Note calling_flags is a YAML list, not a single string, specifically to
# avoid a second layer of shell-quoting/splitting ambiguity when Snakemake's
# shell: formats it into a command line.
#
# The VCF is written to <vcf>.tmp and only moved into place when
# panfreebayes exits 0, so a crashed/killed call never leaves a partial
# VCF that looks like a finished output. A non-zero exit fails the rule with
# panfreebayes's own exit status. PFB_CALL_SHELL is shared with
# call_variants_joint.smk.
#
PFB_CALL_SHELL = r"""
        tmp={output.vcf}.tmp
        rc=0
        {params.pfb} call --ref {input.ref} --bam {input.bam} -- {params.flags} \
            > "$tmp" 2> {log} || rc=$?
        if [ "$rc" -ne 0 ]; then
          rm -f "$tmp"
          echo "panfreebayes call failed with exit status $rc (see {log})" >&2
          exit "$rc"
        fi
        mv "$tmp" {output.vcf}
        """


rule call_variants:
    input:
        ref=f"{RESULTS}/regions/{{bubble}}.fasta",
        bam=f"{RESULTS}/alignments/{{strain}}/{{bubble}}.bam",
        bai=f"{RESULTS}/alignments/{{strain}}/{{bubble}}.bam.bai",
    output:
        vcf=f"{RESULTS}/calls/{{strain}}/{{bubble}}.vcf",
    params:
        pfb=config["panfreebayes_bin"],
        flags=" ".join(config["calling_flags"]),
    log:
        f"{RESULTS}/logs/call_variants/{{strain}}/{{bubble}}.log",
    shell:
        PFB_CALL_SHELL
