#
# Shared helpers: bubble-id naming and checkpoint-TSV lookups.
#
# "checkpoint" (see discover_bubbles.smk) is Snakemake's mechanism for a DAG
# branch whose shape isn't known until a step actually runs -- here, the
# number and identity of qualifying bubbles is only knowable after the awk
# filter has scanned the real deconstruct VCF, not from config.yaml. Every
# downstream rule that needs the bubble list calls
# checkpoints.discover_bubbles.get(), which blocks/re-triggers DAG evaluation
# for that branch until the checkpoint has actually executed, then reads its
# real output.
#


def _bubble_id(path, start, end):
    # Verbatim transplant of extract_region.sh's own output-basename
    # convention: PanSN path with every '#' replaced by '_', plus
    # "_<start>-<end>". Keeping this identical means Snakemake's {bubble}
    # wildcard values look exactly like the filenames the bash script already
    # produces, e.g. DL238_1_chrII_JAFETN010000011.1_1664999-1748623.
    return f"{path.replace('#', '_')}_{start}-{end}"


def _load_bubbles(tsv_path):
    # bubbles.tsv columns (see discover_bubbles.smk):
    #   path  start  end  ref_len  alt_len  size_diff  ref  alt
    rows = {}
    with open(tsv_path) as f:
        next(f)  # header
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            path, start, end = line.split("\t")[:3]
            rows[_bubble_id(path, start, end)] = (path, start, end)
    return rows


def bubble_row(wildcards):
    """Look up (path, start, end) for wildcards.bubble from the checkpoint's TSV."""
    tsv = checkpoints.discover_bubbles.get().output.tsv
    rows = _load_bubbles(tsv)
    return rows[wildcards.bubble]


def all_bubble_ids(wildcards):
    """All bubble ids discovered by the checkpoint -- drives rule all's fan-out."""
    tsv = checkpoints.discover_bubbles.get().output.tsv
    return list(_load_bubbles(tsv).keys())
