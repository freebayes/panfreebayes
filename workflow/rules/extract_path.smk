#
# Whole-path extraction -- an alternative entry point to extract_region.smk's
# bubble-scoped extraction.
#
# extract_region.smk takes a locally-extracted, padded BUBBLE region (found
# via discover_bubbles scanning a vg-deconstruct VCF). This rule instead
# takes one ENTIRE named path from the graph (e.g. one strain's whole
# chromosome) and uses it as a complete, standard reference -- no bubble
# discovery, no size-diff threshold, no padding, nothing to pad around when
# you're already taking the whole thing.
#
# This matches a second, simpler framing of what "panfreebayes" should do:
# pick one path from the pangenome graph, use it as the reference (replacing
# what would normally be a single assembled genome), remap every strain's
# raw reads against it, call variants with the unmodified engine -- as
# opposed to the bubble-local approach extract_region.smk implements. Both
# now coexist; see panfreebayes_milestone6_snakemake_progress.md for why.
#
# Deliberately reuses the exact same odgi mechanism already validated for
# bubble extraction (odgi extract -r "<path>:<start>-<end>" + odgi paths -f),
# confirmed on a real graph to correctly isolate to just the one named path,
# not every path in the subgraph -- the only difference here is the region
# span covers the path's entire length (0 to the path's own length) instead
# of a small bubble window.
#
# NOT YET VALIDATED against a real odgi binary: the `odgi paths -i ... -L`
# invocation below, to list path names and lengths so the full-length -r
# span can be computed, is written from best understanding of odgi's CLI,
# not confirmed against real odgi --help/output. Same honesty standard as
# everything else in this project that could only be checked with real
# tools on vesuvio -- verify the exact flag and output format there before
# trusting this blindly, and adjust the awk column/field below if the real
# output format differs.
#
# Deliberately produces outputs under the SAME regions/{bubble}.{og,fasta}
# path pattern extract_region.smk uses (the wildcard is named {path_id} here,
# but that's just this rule's own internal name for it -- what matters for
# Snakemake's cross-rule matching is the literal resolved file path, not the
# wildcard variable name). align_region/call_variants downstream need ZERO
# changes to support this new extraction mode.
#
rule extract_path:
    input:
        graph=config["graph"],
    output:
        og=f"{RESULTS}/regions/{{path_id}}.og",
        fasta=f"{RESULTS}/regions/{{path_id}}.fasta",
    # Rule-level override: the exact negation of extract_region's bubble
    # constraint (must NOT end in "_<digits>-<digits>"), so the two rules'
    # output patterns are provably non-overlapping -- see extract_region.smk's
    # own comment on this; without it, Snakemake can't tell which rule should
    # produce a given regions/*.fasta path once both exist.
    wildcard_constraints:
        path_id=r"(?!.*_[0-9]+-[0-9]+$).*",
    params:
        path_name=lambda wc: path_name_for(wc.path_id),
    log:
        f"{RESULTS}/logs/extract_path/{{path_id}}.log",
    shell:
        r"""
        length=$(odgi paths -i {input.graph} -L | awk -F'\t' -v p="{params.path_name}" '$1==p {{print $2; exit}}')
        if [ -z "$length" ]; then
          echo "extract_path: path not found in graph (or unexpected 'odgi paths -L' output format -- verify against real odgi): {params.path_name}" >&2
          exit 1
        fi
        odgi extract -i {input.graph} -r "{params.path_name}:0-${{length}}" -o {output.og} > {log} 2>&1
        odgi paths -i {output.og} -f > {output.fasta} 2>> {log}
        """
