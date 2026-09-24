#!/usr/bin/env bash
#
# panFreebayes regression / equivalence harness (Milestone 2).
#
# Milestone 1 acceptance criterion: panfreebayes must produce EXACTLY the same
# VCF as stock single-process freebayes on identical inputs. (The historical
# "1,413 / 2,451" figures are superseded: they predate a July-2026 re-alignment
# of the BAMs. Ground truth is now "whatever current stock freebayes produces".)
#
# Three modes:
#
#   1. equivalence (default) -- run panfreebayes AND freebayes, diff their VCF
#      bodies. Definitive but slow (freebayes can take hours on a coverage spike).
#
#        acceptance_check.sh --ref R.fasta --bam A.bam
#
#   2. regression -- run only panfreebayes, diff against a stored baseline VCF
#      (plain or .gz). Fast; this is the CI / repeatable-check mode.
#
#        acceptance_check.sh --ref R.fasta --bam A.bam --baseline base.vcf.gz
#
#   3. make-baseline -- run stock freebayes and save its VCF as the baseline.
#      Run once to (re)establish ground truth.
#
#        acceptance_check.sh --ref R.fasta --bam A.bam --save-baseline base.vcf.gz
#
# Nothing is hard-coded to a machine layout. Must run inside the guix build env
# so the binaries find their shared libs:
#
#   guix shell --pure bash grep gzip coreutils diffutils gcc-toolchain \
#     zlib xz bzip2 htslib simde tabixpp fastahack smithwaterman intervaltree \
#     -- bash test/panfreebayes/acceptance_check.sh --ref ... --bam ... --baseline ...
#
# Default calling params (the empirically-needed long-read set); override with
# `-- <args>` (which REPLACES them entirely):
#   --pooled-continuous --min-alternate-count 2 --min-alternate-fraction 0.2 --limit-coverage 200
#
set -euo pipefail

here=$(cd "$(dirname "$0")" && pwd)
build_default=$(cd "$here/../.." && pwd)/build

REF="" BAM="" OUTDIR="." BASELINE="" SAVE_BASELINE=""
FB="$build_default/freebayes"
PFB="$build_default/panfreebayes"
EXTRA=()

while [ $# -gt 0 ]; do
  case "$1" in
    --ref)           REF=$2; shift 2 ;;
    --bam)           BAM=$2; shift 2 ;;
    --baseline)      BASELINE=$2; shift 2 ;;
    --save-baseline) SAVE_BASELINE=$2; shift 2 ;;
    --outdir)        OUTDIR=$2; shift 2 ;;
    --freebayes)     FB=$2; shift 2 ;;
    --panfreebayes)  PFB=$2; shift 2 ;;
    --)              shift; EXTRA=("$@"); break ;;
    -h|--help)       sed -n '2,45p' "$0"; exit 0 ;;
    *) echo "unknown arg: $1" >&2; exit 2 ;;
  esac
done

[ -n "$REF" ] && [ -n "$BAM" ] || { echo "need --ref and --bam" >&2; exit 2; }
[ -x "$PFB" ] || { echo "panfreebayes not executable: $PFB" >&2; exit 2; }
[ -f "$REF" ] || { echo "no such ref: $REF" >&2; exit 2; }
[ -f "$BAM" ] || { echo "no such bam: $BAM" >&2; exit 2; }
mkdir -p "$OUTDIR"

if [ "${#EXTRA[@]}" -eq 0 ]; then
  EXTRA=(--pooled-continuous --min-alternate-count 2 --min-alternate-fraction 0.2 --limit-coverage 200)
fi

run_fb="no"
[ -z "$BASELINE" ] && run_fb="yes"
[ -n "$SAVE_BASELINE" ] && run_fb="yes"
if [ "$run_fb" = "yes" ] && [ ! -x "$FB" ]; then
  echo "freebayes needed here but not executable: $FB" >&2; exit 2
fi

base=$(basename "$BAM" .bam)
pfb_vcf="$OUTDIR/${base}.panfreebayes.vcf"
fb_vcf="$OUTDIR/${base}.freebayes.vcf"
diff_out="$OUTDIR/${base}.diff"

body() { grep -v '^##' "$1"; }          # keep the #CHROM line, drop ## meta
read_baseline() {
  case "$1" in
    *.gz) gzip -dc "$1" ;;
    *)    cat "$1" ;;
  esac
}

echo "ref:   $REF"
echo "bam:   $BAM  ($(stat -c '%s bytes, mtime %y' "$BAM" 2>/dev/null || echo '?'))"
echo "args:  ${EXTRA[*]}"
echo

echo "running panfreebayes ..."
time "$PFB" call --ref "$REF" --bam "$BAM" "${EXTRA[@]}" > "$pfb_vcf"
echo "  panfreebayes calls: $(grep -vc '^#' "$pfb_vcf" || true)"
echo

if [ "$run_fb" = "yes" ]; then
  echo "running freebayes ..."
  time "$FB" -f "$REF" "$BAM" "${EXTRA[@]}" > "$fb_vcf"
  echo "  freebayes calls: $(grep -vc '^#' "$fb_vcf" || true)"
  echo
fi

if [ -n "$SAVE_BASELINE" ]; then
  case "$SAVE_BASELINE" in
    *.gz) body "$fb_vcf" | gzip -c > "$SAVE_BASELINE" ;;
    *)    body "$fb_vcf" > "$SAVE_BASELINE" ;;
  esac
  echo "baseline saved: $SAVE_BASELINE  ($(grep -vc '^#' "$fb_vcf") calls, from stock freebayes)"
fi

rc=0
if [ -n "$BASELINE" ]; then
  [ -f "$BASELINE" ] || { echo "no such baseline: $BASELINE" >&2; exit 2; }
  if body "$pfb_vcf" | diff - <(read_baseline "$BASELINE") > "$diff_out"; then
    echo "REGRESSION: panfreebayes == baseline  (PASS)"
  else
    echo "REGRESSION: panfreebayes != baseline  ($(grep -c '^[<>]' "$diff_out") lines -> $diff_out)  (FAIL)"
    rc=1
  fi
elif [ "$run_fb" = "yes" ]; then
  if body "$pfb_vcf" | diff - <(body "$fb_vcf") > "$diff_out"; then
    echo "EQUIVALENCE: panfreebayes == stock freebayes  (PASS)"
  else
    echo "EQUIVALENCE: panfreebayes != stock freebayes  ($(grep -c '^[<>]' "$diff_out") lines -> $diff_out)  (FAIL)"
    rc=1
  fi
fi

exit $rc
