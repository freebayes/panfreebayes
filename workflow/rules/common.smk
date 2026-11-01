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
    """All bubble ids discovered by the checkpoint -- drives rule bubble_mode's fan-out."""
    tsv = checkpoints.discover_bubbles.get().output.tsv
    return list(_load_bubbles(tsv).keys())


# --- Whole-path extraction (extract_path.smk) -------------------------------
#
# A second, simpler entry point alongside bubble extraction: take one whole
# named path from the graph (e.g. one strain's entire chromosome) and use it
# as a complete reference, with no bubble discovery, no size-diff filtering,
# no padding. See extract_path.smk's own header comment for why this exists
# as a separate mode rather than folded into bubble extraction.


import re

# A whole_paths entry is either a bare path name (the whole path) or
# "<path>:<start>-<end>" (1-based, inclusive -- samtools region syntax) to
# keep only that sub-range of the path.
_RANGE_RE = re.compile(r"^(.*):([0-9]+)-([0-9]+)$")


def _parse_path_entry(entry):
    """'<path>[:<start>-<end>]' -> (path, start, end); start/end None for a whole path."""
    m = _RANGE_RE.match(entry)
    if m:
        return m.group(1), int(m.group(2)), int(m.group(3))
    return entry, None, None


def _path_id(entry):
    # Same '#' -> '_' convention as _bubble_id(). A sub-range gets a
    # "_sub<start>-<end>" suffix -- NOT bubble's "_<start>-<end>", which
    # extract_region's wildcard constraint would claim instead.
    path, start, end = _parse_path_entry(entry)
    pid = path.replace("#", "_")
    if start is not None:
        pid += f"_sub{start}-{end}"
    return pid


def path_entry_for(path_id):
    """Reverse-lookup: path_id -> (path, start, end), from config["whole_paths"]."""
    for entry in config.get("whole_paths", []):
        if _path_id(entry) == path_id:
            return _parse_path_entry(entry)
    raise ValueError(f"no configured whole_paths entry matches path_id={path_id}")


def all_path_ids(wildcards):
    """All configured whole-path ids -- drives rule path_mode's fan-out."""
    return [_path_id(p) for p in config.get("whole_paths", [])]
