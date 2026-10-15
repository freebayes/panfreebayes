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
# REVISED after a real run on vesuvio: the first version of this rule used
# `odgi paths -i ... -L` to list path names AND lengths, to compute a
# full-length `-r "<path>:0-<length>"` span for `odgi extract`. That was
# wrong -- confirmed directly, not assumed: `odgi paths -L` only lists path
# NAMES, one per line, no lengths at all. Rather than guess a second time at
# which odgi flag actually reports lengths, this version avoids needing a
# length (or `odgi extract`) at all: `odgi paths -i <full graph> -f` dumps
# FASTA for every path in the graph directly (the exact same mechanism
# already proven to work, just applied to the whole graph instead of a
# pre-extracted subgraph), and a plain `awk` multi-FASTA-record filter keeps
# only the one path we actually want. Simpler than the first version, and
# removes the one place this rule was previously guessing at odgi's CLI.
#
# No `.og` output in this version -- there's no subgraph-extraction step to
# produce one from, and nothing downstream (align_region/call_variants) ever
# reads it; the earlier `.og` output was just extract_region's own byproduct
# carried over by habit, not something actually needed here.
#
# Deliberately produces its FASTA under the SAME regions/{bubble}.fasta path
# pattern extract_region.smk uses (the wildcard is named {path_id} here, but
# that's just this rule's own internal name for it -- what matters for
# Snakemake's cross-rule matching is the literal resolved file path, not the
# wildcard variable name). align_region/call_variants downstream need ZERO
# changes to support this extraction mode.
#
# STILL NOT YET VALIDATED end-to-end against real data (this fix hasn't been
# re-run on vesuvio yet) -- but the one piece of odgi behavior it depends on
# (`odgi paths -f` producing real per-path FASTA records) has now actually
# been exercised successfully on this real graph, unlike the length-lookup
# approach this replaces.
#
rule extract_path:
    input:
        graph=config["graph"],
    output:
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
        odgi paths -i {input.graph} -f 2> {log} | awk -v name=">{params.path_name}" '
          /^>/ {{ keep = ($0 == name) }}
          keep
        ' > {output.fasta}
        if [ ! -s {output.fasta} ]; then
          echo "extract_path: path not found in graph: {params.path_name}" >&2
          exit 1
        fi
        """
