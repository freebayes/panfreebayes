#!/usr/bin/env bash
#
# panfreebayes pipeline -- chain M4 (extract) -> M5 (align) -> M3 (call).
#
# FIRST VERSION, NOT FINAL -- see panfreebayes_milestone4_5_progress.md.
#
#   bubble region on a graph -> odgi extract + odgi paths -f (extract_region.sh)
#                             -> minimap2 + samtools sort/index (align_region.sh)
#                             -> panfreebayes call (the compiled M3 CLI)
#                             -> VCF on stdout
#
# Usage:
#   pipeline.sh --graph graph.og \
#       --region "DL238#1#chrII_JAFETN010000011.1:1664999-1748623" \
#       --fastq strain.fastq \
#       [--padding 5000] [--threads N] [--workdir .] [--panfreebayes PATH] \
#       [-- <flags forwarded to 'panfreebayes call', e.g. --pooled-continuous \
#            --min-alternate-count 2 --min-alternate-fraction 0.2 --limit-coverage 200>] \
#       > out.vcf
#
# Intermediate files (<region-basename>.og / .fasta / .bam[.bai]) are
# written to --workdir (default: current directory) and left in place --
# not cleaned up, matching how the manual pipeline's stages leave their
# outputs around.
#
# Bubble *discovery* (discover_bubbles.sh) is deliberately not part of this
# command -- it's a separate, one-time step over a whole graph that
# produces a list of candidate (path, start, end) triples; this script
# consumes one already-chosen region at a time (one independent invocation
# per bubble per strain -- no chromosome-chunking, by design).
#
set -euo pipefail

here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

usage() {
  cat <<'EOF'
usage: pipeline.sh --graph <graph.og> --region "<path>:<start>-<end>" --fastq <reads.fastq>
                    [--padding 5000] [--threads N] [--workdir .] [--panfreebayes PATH]
                    [-- <flags forwarded to 'panfreebayes call'>]
                    > out.vcf
EOF
}

graph=
region=
fastq=
padding=5000
threads=
workdir=.
panfreebayes=
call_args=()

while [ $# -gt 0 ]; do
  case "$1" in
    --graph)        graph=$2; shift 2 ;;
    --region)       region=$2; shift 2 ;;
    --fastq)        fastq=$2; shift 2 ;;
    --padding)      padding=$2; shift 2 ;;
    --threads)      threads=$2; shift 2 ;;
    --workdir)      workdir=$2; shift 2 ;;
    --panfreebayes) panfreebayes=$2; shift 2 ;;
    -h|--help)      usage; exit 0 ;;
    --)             shift; call_args=("$@"); break ;;
    *)              echo "pipeline.sh: unknown argument: $1" >&2; usage >&2; exit 2 ;;
  esac
done

if [ -z "$graph" ] || [ -z "$region" ] || [ -z "$fastq" ]; then
  echo "pipeline.sh: --graph, --region and --fastq are all required" >&2
  usage >&2
  exit 2
fi

# --region "<path>:<start>-<end>" -- path itself may contain '#' (PanSN) but
# not ':', so split on the LAST ':' and the FIRST '-' after it.
case "$region" in
  *:*-*) ;;
  *) echo "pipeline.sh: --region must look like '<path>:<start>-<end>', got: $region" >&2; exit 2 ;;
esac
region_path="${region%:*}"
coords="${region##*:}"
start="${coords%%-*}"
end="${coords#*-}"
case "$start" in ''|*[!0-9]*) echo "pipeline.sh: could not parse --start from --region: $region" >&2; exit 2 ;; esac
case "$end"   in ''|*[!0-9]*) echo "pipeline.sh: could not parse --end from --region: $region" >&2; exit 2 ;; esac

mkdir -p -- "$workdir"

echo "=== [1/3] extracting region: $region (padding ${padding}bp) ===" >&2
fasta=$("$here/extract_region.sh" --graph "$graph" --path "$region_path" \
            --start "$start" --end "$end" --padding "$padding" --outdir "$workdir")
echo "extracted: $fasta" >&2

base="${region_path//#/_}_${start}-${end}"
bam="$workdir/$base.bam"

echo "=== [2/3] aligning reads: $fastq -> $bam ===" >&2
align_opts=()
[ -n "$threads" ] && align_opts+=(--threads "$threads")
"$here/align_region.sh" --ref "$fasta" --fastq "$fastq" --out "$bam" \
    "${align_opts[@]+"${align_opts[@]}"}" >&2
echo "aligned: $bam" >&2

panfreebayes="${panfreebayes:-$here/../../build/panfreebayes}"
if [ ! -x "$panfreebayes" ]; then
  echo "pipeline.sh: panfreebayes binary not found or not executable: $panfreebayes" >&2
  echo "  (build it, or pass --panfreebayes /path/to/panfreebayes)" >&2
  exit 2
fi

echo "=== [3/3] calling variants ===" >&2
echo "+ $panfreebayes call --ref $fasta --bam $bam ${call_args[*]+"${call_args[*]}"}" >&2
exec "$panfreebayes" call --ref "$fasta" --bam "$bam" "${call_args[@]+"${call_args[@]}"}"
