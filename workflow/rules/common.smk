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


def require_config(key, mode):
    """config[key], or a clear error naming the mode that needs it.

    Used from input/params functions so the check only fires when a rule that
    reads `key` is actually in the requested DAG (e.g. deconstruct_vcf is only
    needed by bubble_mode, not path_mode/haplotype_mode).
    """
    if key not in config or config[key] in (None, ""):
        raise ValueError(f"config key '{key}' is required for {mode}")
    return config[key]


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


# --- Haplotype extraction (extract_haplotype.smk) ---------------------------
#
# A haplotypes entry is a PanSN "<sample>#<hap>" prefix (e.g. "HG00423#1");
# every graph path named "<sample>#<hap>#..." goes into ONE multi-record
# reference FASTA, called as one BAM / one VCF.


def _hap_id(prefix):
    # Same '#' -> '_' convention as _path_id(): "HG00423#1" -> "HG00423_1".
    return prefix.replace("#", "_")


def hap_prefix_for(hap_id):
    """Reverse-lookup: hap_id -> the configured "<sample>#<hap>" prefix."""
    for prefix in config.get("haplotypes") or []:
        if _hap_id(prefix) == hap_id:
            return prefix
    raise ValueError(f"no configured haplotypes entry matches hap_id={hap_id}")


def all_hap_ids(wildcards):
    """All configured haplotype ids -- drives rule haplotype_mode's fan-out."""
    return [_hap_id(p) for p in config.get("haplotypes") or []]


# extract_path and extract_haplotype both write regions/<id>.fasta, so each
# rule's wildcard is constrained to exactly its own configured ids (and,
# like before, never to a bubble-shaped "..._<start>-<end>" id, which
# extract_region owns). An id configured for both would be ambiguous.
_NOT_BUBBLE = r"(?!.*_[0-9]+-[0-9]+$)"


def _ids_regex(ids):
    # "(?!)" never matches: no configured ids means the rule never applies.
    return "(?:" + "|".join(re.escape(i) for i in ids) + ")" if ids else "(?!)"


def _check_region_ids():
    for prefix in config.get("haplotypes") or []:
        if not prefix or prefix.endswith("#") or prefix.count("#") != 1:
            raise ValueError(
                f"haplotypes entry {prefix!r} must be '<sample>#<hap>', e.g. 'HG00423#1'"
            )
    path_ids = {_path_id(p) for p in config.get("whole_paths") or []}
    hap_ids = {_hap_id(p) for p in config.get("haplotypes") or []}
    both = path_ids & hap_ids
    if both:
        raise ValueError(
            "same id configured in both whole_paths and haplotypes: " + ", ".join(sorted(both))
        )


_check_region_ids()
PATH_ID_REGEX = _NOT_BUBBLE + _ids_regex([_path_id(p) for p in config.get("whole_paths") or []])
HAP_ID_REGEX = _NOT_BUBBLE + _ids_regex([_hap_id(p) for p in config.get("haplotypes") or []])


# --- Reads per strain (align_region.smk) ------------------------------------


def strain_reads(strain):
    """The FASTQ(s) for one strain: [fastq] or, for paired-end, [r1, r2].

    A strains entry is either a single FASTQ path (single-end / long reads,
    the original form) or a mapping {r1: ..., r2: ...}.
    """
    entry = config["strains"][strain]
    if isinstance(entry, dict):
        missing = [k for k in ("r1", "r2") if not entry.get(k)]
        if missing or set(entry) - {"r1", "r2"}:
            raise ValueError(
                f"strains entry '{strain}' must be a FASTQ path or {{r1: ..., r2: ...}}"
            )
        return [entry["r1"], entry["r2"]]
    return [entry]


def read_group(strain):
    # minimap2 -R value. \t stays a literal backslash-t here; minimap2
    # expands it. PL only for short reads (preset sr), the one case where
    # the platform is known from the preset alone.
    rg = rf"@RG\tID:{strain}\tSM:{strain}"
    if config.get("minimap2_preset", "map-pb") == "sr":
        rg += r"\tPL:ILLUMINA"
    return rg
