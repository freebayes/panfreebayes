#
# Joint calling (optional, config["joint"]: true) -- one VCF with one sample
# column per strain, instead of one VCF per strain.
#
#   samtools merge <every strain's BAM for this region> -> alignments/joint/<region>.bam
#   panfreebayes call --ref <region.fasta> --bam <merged.bam> -> calls/joint/<region>.vcf
#
# Each strain's reads already carry their own read group (ID/SM = strain
# name, set by align_region), and `samtools merge` keeps those @RG lines and
# RG tags, so the engine sees each strain as its own sample. This is what
# lets a tumour/normal pair be called together with one column each.
#
# "joint" is reserved as a directory name here; the Snakefile's strain
# wildcard constraint excludes it, so these outputs can't be confused with
# align_region's alignments/<strain>/... or call_variants' calls/<strain>/...
#
# Only haplotype_mode requests these outputs; bubble_mode/path_mode are
# unchanged.
#
rule merge_bams:
    input:
        bams=lambda wc: expand(
            f"{RESULTS}/alignments/{{strain}}/{wc.bubble}.bam",
            strain=list(config["strains"].keys()),
        ),
        bais=lambda wc: expand(
            f"{RESULTS}/alignments/{{strain}}/{wc.bubble}.bam.bai",
            strain=list(config["strains"].keys()),
        ),
    output:
        bam=f"{RESULTS}/alignments/joint/{{bubble}}.bam",
        bai=f"{RESULTS}/alignments/joint/{{bubble}}.bam.bai",
    threads: config["minimap2_threads"]
    log:
        f"{RESULTS}/logs/merge_bams/{{bubble}}.log",
    shell:
        r"""
        samtools merge -@ {threads} -f {output.bam} {input.bams} 2> {log}
        samtools index {output.bam} 2>> {log}
        """


rule call_variants_joint:
    input:
        ref=f"{RESULTS}/regions/{{bubble}}.fasta",
        bam=f"{RESULTS}/alignments/joint/{{bubble}}.bam",
        bai=f"{RESULTS}/alignments/joint/{{bubble}}.bam.bai",
    output:
        vcf=f"{RESULTS}/calls/joint/{{bubble}}.vcf",
    params:
        pfb=config["panfreebayes_bin"],
        flags=" ".join(config["calling_flags"]),
    log:
        f"{RESULTS}/logs/call_variants_joint/{{bubble}}.log",
    shell:
        PFB_CALL_SHELL
