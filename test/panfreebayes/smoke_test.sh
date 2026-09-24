#!/usr/bin/env bash
#
# panfreebayes Milestone 1 smoke test.
#
# Not a correctness test against real data. Confirms:
#   1. panfreebayes and freebayes both run to completion on the same inputs
#   2. their VCF output is identical after stripping volatile header lines
#      (i.e. the extracted core reproduces the reference implementation)
#   3. a sharply-localised coverage spike (Milestone 1 risk: buildHaplotypeAlleles
#      window-growth) does not hang, crash, or blow up memory
#
# Requires (add to the guix shell): samtools, python; plus a built ./build/ tree.
#
#   guix shell --pure gcc-toolchain coreutils meson ninja pkg-config cmake perl \
#     samtools python zlib xz bzip2 htslib simde tabixpp fastahack smithwaterman intervaltree
#   ninja -C build
#   test/panfreebayes/smoke_test.sh
#
set -euo pipefail

here=$(cd "$(dirname "$0")" && pwd)
root=$(cd "$here/../.." && pwd)
build=${BUILD_DIR:-$root/build}
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

REF=${REF:-$root/validation_data/DL238_1_chrII_JAFETN010000011.1_1664999-1748589.fasta}
FB=$build/freebayes
PFB=$build/panfreebayes

[ -x "$FB" ]  || { echo "missing $FB — run: ninja -C $build" >&2; exit 1; }
[ -x "$PFB" ] || { echo "missing $PFB — run: ninja -C $build" >&2; exit 1; }
[ -f "$REF" ] || { echo "missing reference $REF" >&2; exit 1; }

cp "$REF" "$work/ref.fasta"; REF="$work/ref.fasta"

TIME=""
command -v /usr/bin/time >/dev/null 2>&1 && TIME="/usr/bin/time -v"

strip() { grep -Ev '^##(commandline|fileDate|reference|source)='; }

FAIL=0
run_case() {
  local tag=$1 bam=$2; shift 2
  echo "=== case: $tag   (args: $*) ==="

  "$FB" -f "$REF" "$@" "$bam" | strip > "$work/$tag.fb.vcf"

  if ! $TIME "$PFB" call --ref "$REF" --bam "$bam" "$@" 2> "$work/$tag.err" \
        | strip > "$work/$tag.pfb.vcf"; then
    echo "  panfreebayes FAILED"; sed -n '1,40p' "$work/$tag.err"; FAIL=1; return
  fi

  local fb pfb
  fb=$(grep -vc '^#' "$work/$tag.fb.vcf"  || true)
  pfb=$(grep -vc '^#' "$work/$tag.pfb.vcf" || true)
  echo "  freebayes calls:    $fb"
  echo "  panfreebayes calls: $pfb"

  if diff -u "$work/$tag.fb.vcf" "$work/$tag.pfb.vcf" > "$work/$tag.diff"; then
    echo "  VCF: IDENTICAL"
  else
    echo "  VCF: DIFFERS  (first 30 diff lines)"; sed -n '1,30p' "$work/$tag.diff"; FAIL=1
  fi
  grep -E 'Maximum resident set size|Elapsed \(wall clock\)' "$work/$tag.err" 2>/dev/null | sed 's/^/  /' || true
}

ACC_ARGS=(--pooled-continuous --min-alternate-count 2 --min-alternate-fraction 0.2 --limit-coverage 200)

python3 "$here/make_synthetic_bam.py" --ref "$REF" --out "$work/flat.bam" \
  --depth 30 --readlen 150 --sample DL238
run_case flat "$work/flat.bam" "${ACC_ARGS[@]}"

python3 "$here/make_synthetic_bam.py" --ref "$REF" --out "$work/spike.bam" \
  --depth 30 --readlen 150 --sample DL238 --spike-at 0.5 --spike-depth 23000 --spike-width 1
run_case spike        "$work/spike.bam" "${ACC_ARGS[@]}"
run_case spike_nolimit "$work/spike.bam" --pooled-continuous --min-alternate-count 2 --min-alternate-fraction 0.2

echo "=== case: region_flag_rejected ==="
if "$PFB" call --ref "$REF" --bam "$work/flat.bam" -- -r foo:1-2 \
      > "$work/region.out" 2> "$work/region.err"; then
  echo "  FAILED: panfreebayes accepted -r instead of rejecting it"; FAIL=1
else
  if grep -q "not supported" "$work/region.err"; then
    echo "  correctly rejected -r (region/target flags must stay unreachable)"
  else
    echo "  FAILED: rejected -r for the wrong reason:"; cat "$work/region.err"; FAIL=1
  fi
fi

echo
if [ "$FAIL" = 0 ]; then
  echo "SMOKE TEST PASSED — extracted core matches freebayes on all cases"
else
  echo "SMOKE TEST FAILED — see diffs above"; exit 1
fi
