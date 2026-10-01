#
# Stage 2 -- region extraction (per bubble).
#
# Matches the manual pipeline's own extraction step exactly:
#
#   odgi extract -i <graph.og> -r "<path>:<start>-<end>" -c <padding> -o <bubble>.og
#   odgi paths -i <bubble>.og -f > <bubble>.fasta
#
# {wildcards.bubble} is exactly the output-basename convention that step used
# (see common.smk's _bubble_id): the PanSN path with every '#' replaced by
# '_', plus "_<start>-<end>". params.path/start/end are looked up from the
# discover_bubbles checkpoint's TSV via bubble_row() (common.smk).
#
# --graph must already be odgi-native (.og); no .gfa->.og conversion is done
# here.
#
rule extract_region:
    input:
        graph=config["graph"],
        tsv=lambda wc: checkpoints.discover_bubbles.get().output.tsv,
    output:
        og=f"{RESULTS}/regions/{{bubble}}.og",
        fasta=f"{RESULTS}/regions/{{bubble}}.fasta",
    params:
        path=lambda wc: bubble_row(wc)[0],
        start=lambda wc: bubble_row(wc)[1],
        end=lambda wc: bubble_row(wc)[2],
        padding=config["padding"],
    log:
        f"{RESULTS}/logs/extract_region/{{bubble}}.log",
    shell:
        r"""
        odgi extract -i {input.graph} -r "{params.path}:{params.start}-{params.end}" \
            -c {params.padding} -o {output.og} > {log} 2>&1
        odgi paths -i {output.og} -f > {output.fasta} 2>> {log}
        """
