# panfreebayes

**panfreebayes** embeds FreeBayes's variant-calling engine as a callable
library + CLI, built for calling variants on individual regions ("bubbles")
extracted from a pangenome graph, one region at a time, for any number of
strains/samples, rather than FreeBayes's usual whole-genome, single-BAM
mode.

The underlying FreeBayes engine itself is not modified; it's extracted and
wrapped, not rewritten. For FreeBayes's own documentation (what it detects,
how it models haplotypes, the full CLI flag reference, citation info), see
[`src/README.md`](src/README.md).

## Repository layout

```
src/panfreebayes/        panfreebayes_core.{h,cpp} (the extracted engine),
                          panfreebayes_main.cpp (the CLI: `panfreebayes call`)
workflow/                 Snakemake workflow, rules/*.smk (see below)
config/config.yaml        workflow configuration (edit this before running)
test/panfreebayes/        regression/validation tests
```

## Building

panfreebayes is built as part of the normal FreeBayes build (meson/ninja),
with one addition: `-Dprefer_system_deps=false`, which switches to this
project's vendored `contrib/vcflib-min` instead of requiring a system-wide
`vcflib` install.

```sh
# a C++17 compiler, meson, ninja, and FreeBayes's usual dependencies
# (htslib, zlib, liblzma, etc.; see src/README.md) need to be available.
meson setup build --buildtype release -Dprefer_system_deps=false
ninja -C build
```

This produces both `build/freebayes` (stock FreeBayes, unmodified) and
`build/panfreebayes` (the CLI wrapper).

## Running panfreebayes directly

```sh
build/panfreebayes call --ref region.fasta --bam region.bam \
  -- --pooled-continuous --min-alternate-count 2 --min-alternate-fraction 0.2 --limit-coverage 200 \
  > calls.vcf
```

Anything after `--` is forwarded verbatim to the underlying FreeBayes
engine: any FreeBayes flag not explicitly modelled by the CLI works this
way. Region/target flags (`-r`/`-t`/`-c`/`-L`/`--stdin`/`--bam-list`) are
rejected deliberately. panfreebayes always analyses its entire given
reference/BAM as one region, by design. Extract the region you want
*before* calling, which is exactly what the workflow below automates.

## Running the Snakemake workflow

The workflow orchestrates the full pipeline (discovering candidate bubble
regions from a `vg deconstruct` VCF, extracting each from the graph,
aligning reads, and calling variants) across any number of strains and any
number of discovered bubbles at once.

**Requirements on `PATH`**: `snakemake`, `odgi`, `minimap2`, `samtools`, and
a built `panfreebayes`/`freebayes` (above). None of these need to come from
the same environment as the build itself; only the already-compiled
binaries matter at run time.

1. **Edit `config/config.yaml`** with your real inputs:
   ```yaml
   graph: "/path/to/graph.og"                     # odgi-native
   deconstruct_vcf: "/path/to/deconstruct.vcf.gz"  # vg deconstruct output
   strains:
     sample1: "/path/to/sample1.fastq"
     sample2: "/path/to/sample2.fastq"
     # any number of strains, any names
   minimap2_preset: "map-pb"   # map-pb (PacBio CLR), map-hifi, map-ont, or sr (short reads)
   panfreebayes_bin: "build/panfreebayes"
   ```
   See the file's own comments for every field (bubble-size threshold,
   extraction padding, calling flags, output directory).

2. **Dry run** to sanity-check before anything executes:
   ```sh
   snakemake all -n --cores 1 --snakefile workflow/Snakefile --configfile config/config.yaml
   ```
   (This can only show the bubble-discovery step up front. The number of
   bubbles found, and everything downstream, isn't known until that step
   actually runs. That's expected, not a bug; see "What each stage does"
   below.)

3. **Real run**:
   ```sh
   snakemake all --cores <N> --snakefile workflow/Snakefile --configfile config/config.yaml
   ```
   Output VCFs land at `<results_dir>/calls/<strain>/<bubble-id>.vcf`, one
   per strain per discovered bubble clearing the configured size threshold.

### What each stage does

The workflow is five Snakemake rules, each in its own file under
`workflow/rules/`, run in this order:

- **`discover_bubbles`** (`discover_bubbles.smk`). Scans the `deconstruct_vcf`
  and, for every site, computes `|REF length - ALT length|`. Sites where
  that exceeds `min_size_diff` (default 500) are kept, sorted largest
  first, into `bubbles.tsv`. This is a Snakemake *checkpoint*, not a plain
  rule, because the number and identity of qualifying bubbles can only be
  known after this step actually runs; nothing else in the workflow can be
  planned ahead of it, which is why a dry run can't show the full DAG
  up front.
- **`extract_region`** (`extract_region.smk`). For one bubble, runs
  `odgi extract` to pull just that region (plus `padding` extra bases on
  each side) out of the graph, then `odgi paths -f` to write it out as a
  plain FASTA. Produces `<bubble>.og` and `<bubble>.fasta`.
- **`align_region`** (`align_region.smk`). For one strain and one bubble,
  aligns that strain's FASTQ reads against the extracted FASTA with
  `minimap2` (preset from `minimap2_preset`), then sorts and indexes the
  result with `samtools`. Produces `<strain>/<bubble>.bam` and its index.
- **`call_variants`** (`call_variants.smk`). For one strain and one bubble,
  runs `panfreebayes call` on the extracted FASTA and aligned BAM, using
  the flags from `calling_flags`. Produces `<strain>/<bubble>.vcf`, the
  workflow's final output.
- **`call_variants_freebayes`** (`call_variants_freebayes.smk`). Optional,
  not part of `rule all`. Runs stock `freebayes` (not `panfreebayes`) on
  the exact same FASTA/BAM a `call_variants` run already produced, purely
  to confirm panfreebayes agrees with the engine it wraps. Only needed when
  checking that extraction itself hasn't changed anything; see "Validating
  against real tools/data" below for how to invoke it.

### Validating against real tools/data

`test/panfreebayes/snakemake/run_real_validation.sh` runs the workflow
end-to-end against real inputs for a single chosen bubble, optionally
diffing the result against a known-good baseline VCF per strain:

```sh
test/panfreebayes/snakemake/run_real_validation.sh \
  --graph /path/to/graph.og \
  --deconstruct-vcf /path/to/deconstruct.vcf.gz \
  --strain sample1=/path/to/sample1.fastq \
  --strain sample2=/path/to/sample2.fastq \
  [--contig SUBSTRING]                      \
  [--baseline sample1=/path/to/baseline.vcf.gz]
```

`--strain` is repeatable (at least one required). `--contig` narrows bubble
discovery to a substring match if given, otherwise considers every
discovered bubble. `--baseline` is repeatable and optional: a strain with
no baseline just gets its call count reported, useful when there's no
existing ground truth yet. Run `--help` for the full flag reference.
