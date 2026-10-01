#!/usr/bin/env bash
#
# panfreebayes Snakemake workflow -- real-tool validation (Milestone 6).
#
# Runs the ACTUAL workflow/Snakefile against real inputs with real odgi,
# minimap2, samtools, and a real built panfreebayes -- no fakes, no stubs.
# This is the real-data validation that panfreebayes_milestone6_snakemake_progress.md
# flags as still outstanding; it supersedes the fake-stub wiring test that
# used to live in this directory (removed once real tools/data became
# available here, since maintaining both would be redundant).
#
# What it does:
#   1. Materialises just the discover_bubbles checkpoint output (bubbles.tsv)
#      against your real --deconstruct-vcf.
#   2. Greps it for the known chrII contig (JAFETN010000011.1) to find
#      whichever bubble id(s) the real discovery step actually produced near
#      the project's nominal bubble coordinates (1664999-1748623) -- this is
#      deliberately NOT hardcoded, because of the known, still-unresolved
#      34bp discrepancy between that nominal range and the committed
#      baseline FASTA's own coordinates (see panfreebayes_milestone4_5_progress.md).
#      If more than one bubble matches, you choose which to validate.
#   3. Builds results/calls/DL238/<bubble>.vcf and results/calls/MY2693/<bubble>.vcf
#      ONLY for that one bubble (not the whole graph -- `snakemake all` would
#      fan out over every bubble the real deconstruct VCF contains, which is
#      not what a quick validation run wants).
#   4. Diffs both against test/panfreebayes/baselines/*.vcf.gz, same
#      pass/fail convention as acceptance_check.sh.
#
# Usage:
#   test/panfreebayes/snakemake/run_real_validation.sh \
#     --graph /path/to/real/graph.og \
#     --deconstruct-vcf /path/to/real/deconstruct.vcf.gz \
#     --dl238-fastq /path/to/DL238.fastq \
#     --my2693-fastq /path/to/MY2693.fastq \
#     [--panfreebayes build/panfreebayes] [--results-dir results_real_validation] \
#     [--bubble <explicit bubble id, skips the contig search/prompt>]
#
# Must run where `snakemake`, `odgi`, `minimap2`, and `samtools` are all on
# PATH (add snakemake to whatever guix shell recipe you already use for
# panfreebayes/odgi/minimap2 -- see the progress doc's Guix section).
#
set -euo pipefail

here=$(cd "$(dirname "$0")" && pwd)
root=$(cd "$here/../../.." && pwd)

GRAPH="" DECONSTRUCT_VCF="" DL238_FASTQ="" MY2693_FASTQ=""
PFB="$root/build/panfreebayes"
RESULTS_DIR="$root/results_real_validation"
EXPLICIT_BUBBLE=""
CONTIG="JAFETN010000011.1"

while [ $# -gt 0 ]; do
  case "$1" in
    --graph)            GRAPH=$2; shift 2 ;;
    --deconstruct-vcf)  DECONSTRUCT_VCF=$2; shift 2 ;;
    --dl238-fastq)      DL238_FASTQ=$2; shift 2 ;;
    --my2693-fastq)     MY2693_FASTQ=$2; shift 2 ;;
    --panfreebayes)     PFB=$2; shift 2 ;;
    --results-dir)      RESULTS_DIR=$2; shift 2 ;;
    --bubble)           EXPLICIT_BUBBLE=$2; shift 2 ;;
    -h|--help)          sed -n '2,40p' "$0"; exit 0 ;;
    *) echo "unknown arg: $1" >&2; exit 2 ;;
  esac
done

for tool in snakemake odgi minimap2 samtools; do
  command -v "$tool" >/dev/null 2>&1 || { echo "required tool not on PATH: $tool" >&2; exit 2; }
done

[ -n "$GRAPH" ] && [ -n "$DECONSTRUCT_VCF" ] && [ -n "$DL238_FASTQ" ] && [ -n "$MY2693_FASTQ" ] \
  || { echo "need --graph, --deconstruct-vcf, --dl238-fastq and --my2693-fastq" >&2; exit 2; }
[ -f "$GRAPH" ]            || { echo "no such graph: $GRAPH" >&2; exit 2; }
[ -f "$DECONSTRUCT_VCF" ]  || { echo "no such deconstruct VCF: $DECONSTRUCT_VCF" >&2; exit 2; }
[ -f "$DL238_FASTQ" ]      || { echo "no such fastq: $DL238_FASTQ" >&2; exit 2; }
[ -f "$MY2693_FASTQ" ]     || { echo "no such fastq: $MY2693_FASTQ" >&2; exit 2; }
[ -x "$PFB" ]              || { echo "panfreebayes not executable: $PFB (build it first)" >&2; exit 2; }

mkdir -p "$RESULTS_DIR"
cfg="$RESULTS_DIR/config.generated.yaml"

# A self-contained config -- not an edited copy of config/config.yaml, so this
# script has no fragile dependency on that file's exact formatting. The
# non-path values below are the project's confirmed defaults (see
# config/config.yaml and panfreebayes_milestone4_5_progress.md); edit here if
# you deliberately want to validate against different tuning.
cat > "$cfg" <<EOF
graph: "$GRAPH"
deconstruct_vcf: "$DECONSTRUCT_VCF"
strains:
  DL238: "$DL238_FASTQ"
  MY2693: "$MY2693_FASTQ"
min_size_diff: 500
bubble_limit:
padding: 5000
minimap2_threads: 4
panfreebayes_bin: "$PFB"
calling_flags:
  - "--pooled-continuous"
  - "--min-alternate-count"
  - "2"
  - "--min-alternate-fraction"
  - "0.2"
  - "--limit-coverage"
  - "200"
results_dir: "$RESULTS_DIR"
EOF

cd "$root"

echo "=== [1/3] discovering bubbles (real deconstruct VCF) ==="
snakemake --cores 1 --snakefile workflow/Snakefile --configfile "$cfg" \
  "$RESULTS_DIR/bubbles/bubbles.tsv"

tsv="$RESULTS_DIR/bubbles/bubbles.tsv"
[ -f "$tsv" ] || { echo "discover_bubbles did not produce $tsv" >&2; exit 1; }

echo
echo "=== [2/3] locating the $CONTIG bubble ==="
matches=$(awk -F'\t' -v contig="$CONTIG" 'NR>1 && $1 ~ contig' "$tsv")
if [ -z "$matches" ]; then
  echo "FAILED: no bubble matching contig '$CONTIG' found in $tsv" >&2
  echo "  (check --deconstruct-vcf actually covers this contig, or that" >&2
  echo "  --min-size-diff in this script's generated config isn't excluding it)" >&2
  exit 1
fi

bubble_id_of() {
  # path start end ... -> path with '#'->'_' plus _<start>-<end>, matching
  # common.smk's _bubble_id() exactly.
  awk -F'\t' '{ p=$1; gsub(/#/,"_",p); print p "_" $2 "-" $3 }' <<<"$1"
}

if [ -n "$EXPLICIT_BUBBLE" ]; then
  bubble="$EXPLICIT_BUBBLE"
  echo "  using explicitly given --bubble $bubble"
else
  n=$(wc -l <<<"$matches" | tr -d ' ')
  if [ "$n" -gt 1 ]; then
    echo "multiple bubbles matched '$CONTIG' -- re-run with one of these via --bubble:" >&2
    while IFS= read -r row; do
      echo "  $(bubble_id_of "$row")   (start=$(cut -f2 <<<"$row") end=$(cut -f3 <<<"$row") size_diff=$(cut -f6 <<<"$row"))" >&2
    done <<<"$matches"
    exit 2
  fi
  bubble=$(bubble_id_of "$matches")
  echo "  found bubble: $bubble"
  echo "  (start=$(cut -f2 <<<"$matches") end=$(cut -f3 <<<"$matches") -- compare against the" \
       "project's nominal 1664999-1748623 and the committed baseline FASTA's own" \
       "1664999-1748589; a mismatch here is the known, already-flagged discrepancy," \
       "not a new bug)"
fi

echo
echo "=== [3/3] extracting, aligning, calling for $bubble (DL238, MY2693) ==="
snakemake --cores "${SNAKEMAKE_CORES:-4}" --snakefile workflow/Snakefile --configfile "$cfg" \
  "$RESULTS_DIR/calls/DL238/$bubble.vcf" "$RESULTS_DIR/calls/MY2693/$bubble.vcf"

body() { grep -v '^##' "$1"; }
read_baseline() { gzip -dc "$1"; }

FAIL=0
check() {
  local strain=$1 baseline=$2
  local vcf="$RESULTS_DIR/calls/$strain/$bubble.vcf"
  local diff_out="$RESULTS_DIR/$strain.$bubble.diff"
  echo
  echo "--- $strain ---"
  echo "  calls: $(grep -vc '^#' "$vcf" || true)"
  if body "$vcf" | diff - <(read_baseline "$baseline") > "$diff_out"; then
    echo "  REGRESSION: $strain == baseline  (PASS)"
  else
    echo "  REGRESSION: $strain != baseline  ($(grep -c '^[<>]' "$diff_out") lines -> $diff_out)  (FAIL)"
    FAIL=1
  fi
}

check DL238  "$root/test/panfreebayes/baselines/DL238_JAFETN010000011.1.vcf.gz"
check MY2693 "$root/test/panfreebayes/baselines/MY2693_JAFETN010000011.1.vcf.gz"

echo
if [ "$FAIL" = 0 ]; then
  echo "REAL-TOOL VALIDATION PASSED -- Snakemake workflow output matches both baselines"
else
  echo "REAL-TOOL VALIDATION FAILED -- see diffs above. Likely suspects, in order:" >&2
  echo "  1. the extracted bubble's actual start/end not matching what produced the" >&2
  echo "     baseline FASTA (the known 34bp discrepancy noted above)" >&2
  echo "  2. --padding (5000 here) not matching what the baseline was actually built with" >&2
  echo "  3. a different underlying graph than whatever produced the original BAMs" >&2
  exit 1
fi
