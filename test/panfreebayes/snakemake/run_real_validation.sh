#!/usr/bin/env bash
#
# panfreebayes Snakemake workflow -- real-tool validation (Milestone 6).
#
# Runs the ACTUAL workflow/Snakefile against real inputs with real odgi,
# minimap2, samtools, and a real built panfreebayes -- no fakes, no stubs.
# Organism/strain-agnostic: takes an arbitrary number of strains and an
# optional contig filter, so this isn't tied to any one species or
# strain-pair (originally built and exercised against two C. elegans
# strains, DL238/MY2693 -- see panfreebayes_milestone6_snakemake_progress.md
# for that history and for the exact flags that reproduce it under this
# generalized CLI).
#
# What it does:
#   1. Materialises just the discover_bubbles checkpoint output (bubbles.tsv)
#      against your real --deconstruct-vcf.
#   2. Locates a bubble to validate: if --contig is given, narrows to rows
#      matching that substring; if omitted, considers every bubble in the
#      file. Either way, --bubble picks a specific one explicitly and skips
#      the search; otherwise a single match is used automatically, and
#      multiple matches are listed for you to choose from via --bubble.
#   3. Builds results/calls/<strain>/<bubble>.vcf for every --strain given,
#      for that one bubble only (not the whole graph -- `snakemake all`
#      would fan out over every bubble the real deconstruct VCF contains,
#      which is not what a quick validation run wants).
#   4. For every --strain that also has a matching --baseline, diffs its
#      VCF against that baseline (same pass/fail convention as
#      acceptance_check.sh). Strains with no baseline given just get their
#      call count reported -- useful when validating a new organism/strain
#      with no existing ground truth yet.
#
# Usage:
#   test/panfreebayes/snakemake/run_real_validation.sh \
#     --graph /path/to/real/graph.og \
#     --deconstruct-vcf /path/to/real/deconstruct.vcf.gz \
#     --strain NAME=/path/to/NAME.fastq [--strain NAME2=/path/to/NAME2.fastq ...] \
#     [--contig SUBSTRING] \
#     [--baseline NAME=/path/to/baseline.vcf.gz ...] \
#     [--panfreebayes build/panfreebayes] [--results-dir results_real_validation] \
#     [--bubble <explicit bubble id, skips the contig search/prompt>]
#
# --strain is repeatable, at least one required. --baseline is repeatable
# and optional; a --baseline NAME only applies to a --strain of the same
# NAME. --contig is optional -- omit it to search every discovered bubble,
# not just ones on a contig you already have to know the name of.
#
# The exact flags that reproduce this script's original two-strain
# C. elegans invocation:
#   --strain DL238=<path> --strain MY2693=<path> --contig JAFETN010000011.1 \
#   --baseline DL238=test/panfreebayes/baselines/DL238_JAFETN010000011.1.vcf.gz \
#   --baseline MY2693=test/panfreebayes/baselines/MY2693_JAFETN010000011.1.vcf.gz
#
# Must run where `snakemake`, `odgi`, `minimap2`, and `samtools` are all on
# PATH (add snakemake to whatever guix shell recipe you already use for
# panfreebayes/odgi/minimap2 -- see the progress doc's Guix section).
#
set -euo pipefail

here=$(cd "$(dirname "$0")" && pwd)
root=$(cd "$here/../../.." && pwd)

GRAPH="" DECONSTRUCT_VCF=""
PFB="$root/build/panfreebayes"
RESULTS_DIR="$root/results_real_validation"
EXPLICIT_BUBBLE=""
CONTIG=""
STRAIN_NAMES=() STRAIN_FASTQS=()
BASELINE_NAMES=() BASELINE_PATHS=()

while [ $# -gt 0 ]; do
  case "$1" in
    --graph)            GRAPH=$2; shift 2 ;;
    --deconstruct-vcf)  DECONSTRUCT_VCF=$2; shift 2 ;;
    --strain)
      STRAIN_NAMES+=("${2%%=*}")
      STRAIN_FASTQS+=("${2#*=}")
      shift 2 ;;
    --baseline)
      BASELINE_NAMES+=("${2%%=*}")
      BASELINE_PATHS+=("${2#*=}")
      shift 2 ;;
    --contig)           CONTIG=$2; shift 2 ;;
    --panfreebayes)     PFB=$2; shift 2 ;;
    --results-dir)      RESULTS_DIR=$2; shift 2 ;;
    --bubble)           EXPLICIT_BUBBLE=$2; shift 2 ;;
    -h|--help)          sed -n '2,56p' "$0"; exit 0 ;;
    *) echo "unknown arg: $1" >&2; exit 2 ;;
  esac
done

for tool in snakemake odgi minimap2 samtools; do
  command -v "$tool" >/dev/null 2>&1 || { echo "required tool not on PATH: $tool" >&2; exit 2; }
done

[ -n "$GRAPH" ] && [ -n "$DECONSTRUCT_VCF" ] \
  || { echo "need --graph and --deconstruct-vcf" >&2; exit 2; }
[ "${#STRAIN_NAMES[@]}" -gt 0 ] \
  || { echo "need at least one --strain NAME=FASTQ" >&2; exit 2; }
[ -f "$GRAPH" ]            || { echo "no such graph: $GRAPH" >&2; exit 2; }
[ -f "$DECONSTRUCT_VCF" ]  || { echo "no such deconstruct VCF: $DECONSTRUCT_VCF" >&2; exit 2; }
# C-style indexing throughout (not "${!arr[@]}"/"${arr[@]}") -- these arrays
# can legitimately be empty (no --baseline given at all), and old bash
# (3.2, macOS's system default) crashes with "unbound variable" under set -u
# when @-expanding a declared-but-empty array; "${#arr[@]}" + C-style ((...))
# avoids that entirely, portable on any bash.
for ((i = 0; i < ${#STRAIN_NAMES[@]}; i++)); do
  [ -f "${STRAIN_FASTQS[$i]}" ] || { echo "no such fastq for ${STRAIN_NAMES[$i]}: ${STRAIN_FASTQS[$i]}" >&2; exit 2; }
done
for ((i = 0; i < ${#BASELINE_PATHS[@]}; i++)); do
  [ -f "${BASELINE_PATHS[$i]}" ] || { echo "no such baseline for ${BASELINE_NAMES[$i]}: ${BASELINE_PATHS[$i]}" >&2; exit 2; }
done
[ -x "$PFB" ] || { echo "panfreebayes not executable: $PFB (build it first)" >&2; exit 2; }

baseline_for() {
  local name=$1
  # Note: deliberately `if`, not a bare `[ ... ] && { ...; }` statement --
  # under set -e, the latter aborts the whole script on the first
  # NON-matching entry (a bare "test && action" with no else branch has the
  # test's own failure as its exit status, which set -e does not exempt the
  # way it exempts an if/while condition).
  for ((i = 0; i < ${#BASELINE_NAMES[@]}; i++)); do
    if [ "${BASELINE_NAMES[$i]}" = "$name" ]; then
      echo "${BASELINE_PATHS[$i]}"
      return
    fi
  done
}

mkdir -p "$RESULTS_DIR"
cfg="$RESULTS_DIR/config.generated.yaml"

# A self-contained config -- not an edited copy of config/config.yaml, so this
# script has no fragile dependency on that file's exact formatting. The
# non-path values below are the project's confirmed defaults (see
# config/config.yaml and panfreebayes_milestone4_5_progress.md); edit here if
# you deliberately want to validate against different tuning.
{
  echo "graph: \"$GRAPH\""
  echo "deconstruct_vcf: \"$DECONSTRUCT_VCF\""
  echo "strains:"
  for ((i = 0; i < ${#STRAIN_NAMES[@]}; i++)); do
    echo "  ${STRAIN_NAMES[$i]}: \"${STRAIN_FASTQS[$i]}\""
  done
  cat <<EOF
min_size_diff: 500
bubble_limit:
padding: 5000
minimap2_threads: 4
minimap2_preset: "map-pb"
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
} > "$cfg"

cd "$root"

echo "=== [1/3] discovering bubbles (real deconstruct VCF) ==="
# Target(s) given BEFORE --configfile deliberately: Snakemake 8.x made
# --configfile nargs='+' (to support layering multiple config files), so it
# now greedily swallows whatever token follows it on the command line -- if
# the target path came after --configfile here, Snakemake would try to
# open() the (not-yet-created) target as a second config file and fail with
# a confusing FileNotFoundError, instead of treating it as a build target.
snakemake "$RESULTS_DIR/bubbles/bubbles.tsv" \
  --cores 1 --snakefile workflow/Snakefile --configfile "$cfg"

tsv="$RESULTS_DIR/bubbles/bubbles.tsv"
[ -f "$tsv" ] || { echo "discover_bubbles did not produce $tsv" >&2; exit 1; }

echo
if [ -n "$CONTIG" ]; then
  echo "=== [2/3] locating the $CONTIG bubble ==="
else
  echo "=== [2/3] locating a bubble (no --contig given, considering all) ==="
fi
matches=$(awk -F'\t' -v contig="$CONTIG" 'NR>1 && (contig=="" || $1 ~ contig)' "$tsv")
if [ -z "$matches" ]; then
  echo "FAILED: no bubble found${CONTIG:+ matching contig '$CONTIG'} in $tsv" >&2
  echo "  (check --deconstruct-vcf actually covers what you expect, or that" >&2
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
    echo "multiple bubbles matched -- re-run with one of these via --bubble:" >&2
    while IFS= read -r row; do
      echo "  $(bubble_id_of "$row")   (start=$(cut -f2 <<<"$row") end=$(cut -f3 <<<"$row") size_diff=$(cut -f6 <<<"$row"))" >&2
    done <<<"$matches"
    exit 2
  fi
  bubble=$(bubble_id_of "$matches")
  echo "  found bubble: $bubble"
  echo "  (start=$(cut -f2 <<<"$matches") end=$(cut -f3 <<<"$matches") -- if this doesn't match" \
       "coordinates you expected from elsewhere, that's worth checking: the real" \
       "discovered bubble end doesn't always match a nominal/expected one exactly)"
fi

echo
echo "=== [3/3] extracting, aligning, calling for $bubble (${STRAIN_NAMES[*]}) ==="
targets=()
for name in "${STRAIN_NAMES[@]}"; do
  targets+=("$RESULTS_DIR/calls/$name/$bubble.vcf")
done
# Same Snakemake 8.x --configfile-is-greedy ordering fix as above.
snakemake "${targets[@]}" \
  --cores "${SNAKEMAKE_CORES:-4}" --snakefile workflow/Snakefile --configfile "$cfg"

body() { grep -v '^##' "$1"; }
read_baseline() { gzip -dc "$1"; }

FAIL=0
check() {
  local strain=$1 baseline=$2
  local vcf="$RESULTS_DIR/calls/$strain/$bubble.vcf"
  echo
  echo "--- $strain ---"
  echo "  calls: $(grep -vc '^#' "$vcf" || true)"
  if [ -z "$baseline" ]; then
    echo "  (no --baseline given for $strain -- call count only, not diffed)"
    return
  fi
  local diff_out="$RESULTS_DIR/$strain.$bubble.diff"
  if body "$vcf" | diff - <(read_baseline "$baseline") > "$diff_out"; then
    echo "  REGRESSION: $strain == baseline  (PASS)"
  else
    echo "  REGRESSION: $strain != baseline  ($(grep -c '^[<>]' "$diff_out") lines -> $diff_out)  (FAIL)"
    FAIL=1
  fi
}

for name in "${STRAIN_NAMES[@]}"; do
  check "$name" "$(baseline_for "$name")"
done

echo
if [ "$FAIL" = 0 ]; then
  echo "REAL-TOOL VALIDATION PASSED -- Snakemake workflow output matches every given baseline"
else
  echo "REAL-TOOL VALIDATION FAILED -- see diffs above. Likely suspects, in order:" >&2
  echo "  1. the extracted bubble's actual start/end not matching what produced the" >&2
  echo "     baseline (extraction padding/coordinates differ from whatever made it)" >&2
  echo "  2. --padding (5000 here) not matching what the baseline was actually built with" >&2
  echo "  3. a different underlying graph than whatever produced the original BAMs" >&2
  exit 1
fi
