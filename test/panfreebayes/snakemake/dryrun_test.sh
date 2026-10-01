#!/usr/bin/env bash
#
# panfreebayes Snakemake workflow -- wiring test.
#
# Not a correctness test against real data (same caveat as
# test/panfreebayes/smoke_test.sh and the bash pipeline's own fake-stub
# tests). Confirms, by substituting logging fake odgi/minimap2/samtools/
# panfreebayes executables on PATH:
#   1. `snakemake -n` resolves the pre-checkpoint DAG correctly against the
#      synthetic fixtures (discover_bubbles + all only -- Snakemake cannot
#      show the downstream extract/align/call fan-out in a dry run, since it
#      depends on bubbles.tsv, which only exists after the checkpoint
#      actually runs; this is expected checkpoint behaviour, not a gap)
#   2. a real `snakemake --cores N` run against those fakes actually fans out
#      over every (strain, bubble) pair (2 bubbles x 2 strains here) and
#      produces a VCF for each
#   3. the exact command construction is correct: odgi extract's -r/-c/-o,
#      odgi paths -f, minimap2's -ax map-pb -t, samtools sort -o / index,
#      and panfreebayes call's --ref/--bam/-- <calling flags>
#
# Requires `snakemake` on PATH (add it to the guix shell package list used
# for panfreebayes work -- see panfreebayes_milestone6_snakemake_progress.md).
# Does NOT require odgi, minimap2, samtools, or a real panfreebayes build --
# those are all faked. Real-data validation remains a separate, outstanding
# step on vesuvio (see the progress doc).
#
set -euo pipefail

here=$(cd "$(dirname "$0")" && pwd)
root=$(cd "$here/../../.." && pwd)
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

command -v snakemake >/dev/null 2>&1 || { echo "missing snakemake on PATH -- see this script's header comment" >&2; exit 1; }

export PATH="$here/fakebin:$PATH"
export FAKE_LOG="$work/fake_invocations.log"
touch "$FAKE_LOG"

cd "$root"

FAIL=0

echo "=== dry run ==="
snakemake all -n --cores 1 --snakefile workflow/Snakefile \
  --configfile "$here/config.yaml" \
  --config results_dir="$work/results" | tee "$work/dryrun.out"

# Snakemake dry-runs a checkpoint's DAG only up to the checkpoint itself --
# it cannot know the downstream extract/align/call fan-out until
# discover_bubbles has actually run and produced bubbles.tsv, so the correct
# expectation here is exactly "checkpoint discover_bubbles" + "all", plus
# Snakemake's own note that the DAG will grow after the checkpoint runs.
if grep -q 'checkpoint discover_bubbles' "$work/dryrun.out" \
   && grep -q 'localrule all' "$work/dryrun.out" \
   && grep -q 'alteration of the DAG of jobs' "$work/dryrun.out"; then
  echo "  dry run correctly shows only the checkpoint + all pre-expansion"
else
  echo "  FAILED: dry run output didn't look like the expected pre-checkpoint DAG"; FAIL=1
fi

echo
echo "=== real run against fake odgi/minimap2/samtools/panfreebayes ==="
if ! snakemake all --cores 2 --snakefile workflow/Snakefile \
      --configfile "$here/config.yaml" \
      --config results_dir="$work/results" > "$work/run.out" 2> "$work/run.err"; then
  echo "  FAILED: snakemake run did not complete"; sed -n '1,60p' "$work/run.err"; FAIL=1
fi

BUBBLES=(
  "DL238_1_chrII_JAFETN010000011.1_200-200"
  "DL238_1_chrII_JAFETN010000011.1_500-527"
)
STRAINS=(DL238 MY2693)

for strain in "${STRAINS[@]}"; do
  for bubble in "${BUBBLES[@]}"; do
    vcf="$work/results/calls/$strain/$bubble.vcf"
    if [ -s "$vcf" ]; then
      echo "  VCF present: calls/$strain/$bubble.vcf"
    else
      echo "  FAILED: missing or empty $vcf"; FAIL=1
    fi
  done
done

echo
echo "=== checking exact command construction in $FAKE_LOG ==="
check() {
  local desc=$1 pattern=$2
  if grep -qE "$pattern" "$FAKE_LOG"; then
    echo "  ok: $desc"
  else
    echo "  FAILED: expected invocation not found: $desc"; FAIL=1
  fi
}

check "odgi extract with -r <path:start-end> and -c <padding>" \
  '^odgi extract -i .* -r DL238#1#chrII_JAFETN010000011\.1:200-200 -c 100 -o .*\.og$'
check "odgi paths -f" \
  '^odgi paths -i .*\.og -f$'
check "minimap2 -ax map-pb -t <threads>" \
  '^minimap2 -ax map-pb -t 1 .*\.fasta .*\.fastq$'
check "samtools sort -o" \
  '^samtools sort -o .*\.bam$'
check "samtools index" \
  '^samtools index .*\.bam$'
check "panfreebayes call with -- and the four calling flags" \
  '^panfreebayes call --ref .*\.fasta --bam .*\.bam -- --pooled-continuous --min-alternate-count 2 --min-alternate-fraction 0\.2 --limit-coverage 200$'

echo
if [ "$FAIL" = 0 ]; then
  echo "WIRING TEST PASSED -- DAG, fan-out, and exact command construction all match"
else
  echo "WIRING TEST FAILED -- see above"; exit 1
fi
