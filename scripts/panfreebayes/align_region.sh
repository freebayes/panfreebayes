#!/usr/bin/env bash
#
# Milestone 5 -- alignment wrapper (minimap2 + samtools sort/index).
#
# FIRST VERSION, NOT FINAL -- see panfreebayes_milestone4_5_progress.md.
# Matches 14_align_reads_to_region.sh exactly:
#
#   minimap2 -ax map-pb -t <threads> <ref.fasta> <fastq> | samtools sort -o <out>.bam
#   samtools index <out>.bam
#
# No flags beyond what's shown were used in the original -- none are added
# here. The only change from the original: thread count is a --threads flag
# (default: nproc, falling back to 4 if nproc isn't available) instead of
# the hardcoded 16. Threading applies to minimap2 only (-t); `samtools sort`
# is invoked exactly as in the original, with no -@.
#
# set -o pipefail (below) means a minimap2 failure fails this script even
# though it's the left side of a pipe into samtools sort -- this is the
# bash-native equivalent of checking both processes' exit codes.
#
# Usage:
#   align_region.sh --ref region.fasta --fastq strain.fastq --out aligned.bam [--threads N]
#
set -euo pipefail

usage() {
  cat <<'EOF'
usage: align_region.sh --ref <region.fasta> --fastq <strain.fastq> --out <out.bam> [--threads N]
EOF
}

ref=
fastq=
out=
threads=

while [ $# -gt 0 ]; do
  case "$1" in
    --ref)     ref=$2; shift 2 ;;
    --fastq)   fastq=$2; shift 2 ;;
    --out)     out=$2; shift 2 ;;
    --threads) threads=$2; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *)         echo "align_region.sh: unknown argument: $1" >&2; usage >&2; exit 2 ;;
  esac
done

if [ -z "$ref" ] || [ -z "$fastq" ] || [ -z "$out" ]; then
  echo "align_region.sh: --ref, --fastq and --out are all required" >&2
  usage >&2
  exit 2
fi

command -v minimap2 >/dev/null 2>&1 || { echo "align_region.sh: required tool 'minimap2' not found on PATH" >&2; exit 2; }
command -v samtools  >/dev/null 2>&1 || { echo "align_region.sh: required tool 'samtools' not found on PATH" >&2; exit 2; }
[ -f "$ref" ]   || { echo "align_region.sh: no such reference: $ref" >&2; exit 2; }
[ -f "$fastq" ] || { echo "align_region.sh: no such fastq: $fastq" >&2; exit 2; }

if [ -z "$threads" ]; then
  threads=$(nproc 2>/dev/null || echo 4)
fi

minimap2 -ax map-pb -t "$threads" "$ref" "$fastq" | samtools sort -o "$out"
samtools index "$out"

echo "$out"
