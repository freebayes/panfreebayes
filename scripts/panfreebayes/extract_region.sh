#!/usr/bin/env bash
#
# Milestone 4, part 2 -- region extraction (odgi extract + odgi paths -f).
#
# FIRST VERSION, NOT FINAL -- see panfreebayes_milestone4_5_progress.md.
# Matches 13_extract_bubble_region.sh exactly:
#
#   odgi extract -i <graph.og> -r "<path>:<start>-<end>" -c <padding> -o <base>.og
#   odgi paths -i <base>.og -f > <base>.fasta
#
# No PanSN header rewriting is done here, by design -- that happens
# upstream of this script, before pggb ever runs. --path is expected to
# already be a correctly-qualified PanSN path string (exactly what
# discover_bubbles.sh's "path" column provides); this script does not
# inspect or rewrite it.
#
# --graph must already be in odgi's native format (.og) -- `odgi extract -i`
# takes .og, not a raw .gfa. Converting a .gfa (`odgi build`) is not done
# here; not part of the script this was matched against.
#
# Usage:
#   extract_region.sh --graph graph.og --path "DL238#1#chrII_JAFETN010000011.1" \
#       --start 1664999 --end 1748623 [--padding 5000] [--outdir .]
#
# Output naming (matches the manual pipeline exactly): the PanSN path with
# every '#' replaced by '_', plus "_<start>-<end>", e.g.:
#   DL238#1#chrII_JAFETN010000011.1 : 1664999-1748623
#     -> DL238_1_chrII_JAFETN010000011.1_1664999-1748623.{og,fasta}
#
# Prints the resulting .fasta path to stdout on success; the .og path to
# stderr (informational).
#
set -euo pipefail

usage() {
  cat <<'EOF'
usage: extract_region.sh --graph <graph.og> --path <PanSN-path> --start N --end N
                          [--padding 5000] [--outdir .]
EOF
}

graph=
region_path=
start=
end=
padding=5000
outdir=.

while [ $# -gt 0 ]; do
  case "$1" in
    --graph)   graph=$2; shift 2 ;;
    --path)    region_path=$2; shift 2 ;;
    --start)   start=$2; shift 2 ;;
    --end)     end=$2; shift 2 ;;
    --padding) padding=$2; shift 2 ;;
    --outdir)  outdir=$2; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *)         echo "extract_region.sh: unknown argument: $1" >&2; usage >&2; exit 2 ;;
  esac
done

if [ -z "$graph" ] || [ -z "$region_path" ] || [ -z "$start" ] || [ -z "$end" ]; then
  echo "extract_region.sh: --graph, --path, --start and --end are all required" >&2
  usage >&2
  exit 2
fi

case "$start" in ''|*[!0-9]*) echo "extract_region.sh: --start must be a non-negative integer, got: $start" >&2; exit 2 ;; esac
case "$end"   in ''|*[!0-9]*) echo "extract_region.sh: --end must be a non-negative integer, got: $end" >&2; exit 2 ;; esac

if [ "$start" -ge "$end" ]; then
  echo "extract_region.sh: --start ($start) must be < --end ($end)" >&2
  exit 2
fi

command -v odgi >/dev/null 2>&1 || { echo "extract_region.sh: required tool 'odgi' not found on PATH" >&2; exit 2; }
[ -f "$graph" ] || { echo "extract_region.sh: no such graph file: $graph" >&2; exit 2; }

mkdir -p -- "$outdir"

base="${region_path//#/_}_${start}-${end}"
og="$outdir/$base.og"
fasta="$outdir/$base.fasta"
region="${region_path}:${start}-${end}"

odgi extract -i "$graph" -r "$region" -c "$padding" -o "$og"
odgi paths -i "$og" -f > "$fasta"

echo "$og" >&2
echo "$fasta"
