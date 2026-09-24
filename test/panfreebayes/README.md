# panFreebayes test harness

## The CLI (Milestone 3)

```
panfreebayes call --ref <ref.fasta> --bam <aln.bam> [flags] > out.vcf
panfreebayes version
panfreebayes help
```

`call` analyses the **whole** reference/BAM as one region (no `-r`/`-t`). Flags:
`--pooled-continuous`, `--min-alternate-count N`, `--min-alternate-fraction X`,
`--limit-coverage N`; `-- <args>` passes extra tokens to the engine verbatim
(region/target/stdin flags are rejected). Output is bit-for-bit identical to
stock single-process `freebayes` with the equivalent arguments.

Pre-flight checks before calling: `--ref`/`--bam` readable; BAM opens and has a
header; BAM is `SO:coordinate`; at least one BAM `@SQ` name is present in the
FASTA (else "BAM and reference do not match"). Exit codes: `0` ok, `1` bad
command, `2` bad arguments / failed pre-flight.

---

## Test scripts

Both meant to run **inside the guix build environment** (so the built binaries
find their shared libraries):

```sh
guix shell --pure bash grep gzip coreutils diffutils gcc-toolchain \
  meson ninja pkg-config cmake perl samtools python \
  zlib xz bzip2 htslib simde tabixpp fastahack smithwaterman intervaltree
ninja -C build
```

---

## `smoke_test.sh` — synthetic, fast, no external data

Generates synthetic BAMs (flat coverage + a ~23,000× spike) against the bundled
`validation_data/*.fasta`, runs both `freebayes` and `panfreebayes` with the
standard long-read params, and asserts their VCFs are byte-identical. Also checks
the spike case doesn't hang/OOM and that `-r`/`-t` flags are rejected.

```sh
test/panfreebayes/smoke_test.sh
```

Perfect synthetic reads produce ~0 variant calls, so this proves the *plumbing*
(argv, headers, sample handling, no crash) — not the calling path. Use it as the
CI check and the post-build sanity check.

---

## `acceptance_check.sh` — real data, equivalence + regression

Milestone 1 acceptance criterion: **panfreebayes must produce exactly the same
VCF as stock single-process `freebayes`** on identical inputs.

> The historical "1,413 / 2,451" call counts are **superseded**. They predate a
> 2026-07-17 re-alignment of the strain BAMs (minimap2 2.28, `-ax map-pb`).
> Ground truth is now "whatever current stock `freebayes` produces" — which is
> **4,102** for DL238 on the current BAM.

### 1. Establish / refresh the baseline (slow, run once)

```sh
test/panfreebayes/acceptance_check.sh \
  --ref  /path/to/DL238_..._1664999-1748589.fasta \
  --bam  /path/to/DL238_DL238_..._1664999-1748589.bam \
  --save-baseline test/panfreebayes/baselines/DL238_JAFETN010000011.1.vcf.gz
```

Runs `freebayes`, saves its VCF body (gzip). `freebayes` can take **hours** on a
region with a coverage spike — this is why the baseline exists.

Record what the baseline corresponds to alongside it:

```sh
{ echo "bam_md5=$(md5sum <bam> | cut -d' ' -f1)"
  echo "made=$(date -u +%FT%TZ)  freebayes=$(freebayes --version)"
  samtools view -H <bam> | grep '^@PG'
} > test/panfreebayes/baselines/DL238_JAFETN010000011.1.provenance
```

### 2. Regression check (fast, repeatable)

```sh
test/panfreebayes/acceptance_check.sh \
  --ref  /path/to/....fasta \
  --bam  /path/to/....bam \
  --baseline test/panfreebayes/baselines/DL238_JAFETN010000011.1.vcf.gz
```

Runs only `panfreebayes`, diffs its VCF body against the baseline.
Exit 0 = `PASS`, exit 1 = `FAIL` (diff written to `<bam>.diff` in `--outdir`,
default `.`).

### 3. Full equivalence (no baseline — runs both)

```sh
test/panfreebayes/acceptance_check.sh --ref ....fasta --bam ....bam
```

`-- <args>` overrides the default calling params
(`--pooled-continuous --min-alternate-count 2 --min-alternate-fraction 0.2 --limit-coverage 200`)
entirely, passed verbatim to both tools.

---

## Status

- [x] `panfreebayes` ≡ stock `freebayes`, byte-identical, on synthetic data (smoke test)
- [x] `panfreebayes` ≡ stock `freebayes`, byte-identical, on the real DL238 bubble
      (4,102 calls, 2026-09-08) — baseline committed under `baselines/`
- [x] M3 CLI (`panfreebayes call`) + pre-flight checks; smoke test runs through it
- [ ] MY2693 bubble baseline + regression (run in progress)
- [ ] DL238 regression re-run *through `panfreebayes call`* (belt-and-suspenders;
      the smoke test already exercises `call`, and identical `Options` → identical
      engine invocation, so this only re-confirms the arg parsing)
- [ ] second validation region (`…JAFETN010000040.1:497444-552651`) checked
