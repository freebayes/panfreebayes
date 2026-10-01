#!/usr/bin/env bash
#
# Milestone 4, part 1 -- bubble discovery.
#
# FIRST VERSION, NOT FINAL -- see panfreebayes_milestone4_5_progress.md.
# This formalizes the manual pipeline's own awk filter: sites where a
# structural event's REF and ALT alleles differ in length by MORE than
# --min-size-diff (strict, matching the original `diff > 500`), sorted by
# that size difference, largest first. 500 is the confirmed threshold from
# that filter, not a guess.
#
# Reads `vg deconstruct` output (a VCF, plain or .gz). A vg-deconstruct
# VCF's CHROM is the reference PATH it was deconstructed against (e.g. a
# PanSN path like "DL238#1#chrII_JAFETN010000011.1", already correctly
# qualified upstream of this script -- no header rewriting is done here).
#
# Multi-allelic handling (pick the ALT with the largest size difference
# from REF; skip symbolic ALTs like <INS> rather than guess their length)
# is THIS SCRIPT'S OWN design choice, not ported from the original -- the
# manual filter never actually encountered multi-allelic sites in practice.
# See panfreebayes_milestone4_5_progress.md.
#
# Usage:
#   discover_bubbles.sh <deconstruct.vcf[.gz]> [--min-size-diff 500] [--limit N]
#
# Output (TSV to stdout), one row per qualifying site, largest event first:
#   path  start  end  ref_len  alt_len  size_diff  ref  alt
#
# start/end are 1-based inclusive (VCF convention), matching what
# extract_region.sh's -r "<path>:<start>-<end>" expects directly.
#
set -euo pipefail

usage() {
  cat <<'EOF'
usage: discover_bubbles.sh <deconstruct.vcf[.gz]> [--min-size-diff 500] [--limit N]

  --min-size-diff N   report sites where |REF_len - ALT_len| is STRICTLY
                      greater than N (default: 500, the confirmed manual
                      pipeline threshold: diff > 500)
  --limit N           report at most N sites (largest size-diff first)
EOF
}

min_size_diff=500
limit=
vcf=

while [ $# -gt 0 ]; do
  case "$1" in
    --min-size-diff) min_size_diff=$2; shift 2 ;;
    --limit)         limit=$2; shift 2 ;;
    -h|--help)       usage; exit 0 ;;
    -*)              echo "discover_bubbles.sh: unknown flag: $1" >&2; usage >&2; exit 2 ;;
    *)               vcf=$1; shift ;;
  esac
done

[ -n "$vcf" ]   || { echo "discover_bubbles.sh: missing <deconstruct.vcf[.gz]>" >&2; usage >&2; exit 2; }
[ -f "$vcf" ]   || { echo "discover_bubbles.sh: no such file: $vcf" >&2; exit 2; }

reader=(cat -- "$vcf")
case "$vcf" in
  *.gz) reader=(gzip -dc -- "$vcf") ;;
esac

tmp=$(mktemp)
trap 'rm -f "$tmp" "$tmp.limited"' EXIT

"${reader[@]}" | awk -v min_size_diff="$min_size_diff" '
  BEGIN { FS = "\t"; OFS = "\t" }
  /^#/  { next }
  NF < 5 { next }
  {
    chrom = $1; pos = $2 + 0; ref = $4; alt = $5
    ref_len = length(ref)

    n = split(alt, alts, ",")
    best_diff = -1; best_alt = ""; best_alt_len = 0
    for (i = 1; i <= n; i++) {
      a = alts[i]
      if (a == "" || substr(a, 1, 1) == "<") continue   # symbolic ALT -- skip, do not guess a length
      alt_len = length(a)
      diff = ref_len - alt_len
      if (diff < 0) diff = -diff
      if (diff > best_diff) { best_diff = diff; best_alt = a; best_alt_len = alt_len }
    }
    if (best_diff < 0) next                # no usable (non-symbolic) ALT
    if (best_diff <= min_size_diff) next    # strict: matches the original `diff > 500`

    start = pos
    end   = pos + ref_len - 1
    print chrom, start, end, ref_len, best_alt_len, best_diff, ref, best_alt
  }
' | sort -t "$(printf '\t')" -k6,6nr > "$tmp"

if [ -n "$limit" ]; then
  head -n "$limit" "$tmp" > "$tmp.limited"
  mv "$tmp.limited" "$tmp"
fi

printf 'path\tstart\tend\tref_len\talt_len\tsize_diff\tref\talt\n'
cat "$tmp"

echo "$(wc -l < "$tmp" | tr -d ' ') candidate site(s) with size-diff > ${min_size_diff}" >&2
