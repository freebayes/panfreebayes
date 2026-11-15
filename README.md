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
a built `panfreebayes`/`freebayes` (above), plus `vg` for haplotype mode on
a `.gbz` graph. None of these need to come from
the same environment as the build itself; only the already-compiled
binaries matter at run time.

**The `pfb` command**: `scripts/pfb` is a short front-end for the workflow. It
fills in `--snakefile` and (unless you pass your own) `--configfile` from
this checkout, so it works from any directory. Put it on `PATH` with a
symlink (not a copy, since it finds the workflow next to its own real location):
```sh
ln -s "$PWD/scripts/pfb" ~/.local/bin/pfb     # or: export PATH="$PWD/scripts:$PATH"
```
`pfb bubble` runs `bubble_mode`, `pfb path` runs `path_mode`, and `pfb bubble path`
runs both. Every other argument (`-n`, `--cores`, `--configfile my.yaml`, ...)
goes straight to snakemake. `results_dir` is relative to the directory you
run `pfb` from. Plain `snakemake <target> --snakefile workflow/Snakefile ...`
still works exactly as before.

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
   pfb bubble -n --cores 1
   ```
   (This can only show the bubble-discovery step up front. The number of
   bubbles found, and everything downstream, isn't known until that step
   actually runs. That's expected, not a bug; see "What each stage does"
   below.)

3. **Real run**:
   ```sh
   pfb bubble --cores <N>
   ```
   Output VCFs land at `<results_dir>/calls/<strain>/<bubble-id>.vcf`, one
   per strain per discovered bubble clearing the configured size threshold.

**Whole-path mode, instead of bubble discovery**: set `whole_paths` in
`config/config.yaml` to one or more path names from the graph, then run
`pfb path --cores <N>` in place of `pfb bubble`. This
extracts each listed path in full and uses it as a complete reference, no
bubble discovery or padding involved. See "What each stage does" below.

### Haplotype mode (e.g. HG008 tumour/normal)

Uses **every path of one haplotype** as a single reference: for each
`"<sample>#<hap>"` prefix in `haplotypes:` (e.g. `"HG00423#1"`), all graph
paths named `<sample>#<hap>#...` are extracted into one multi-contig FASTA,
`<results_dir>/regions/<sample>_<hap>.fasta`. Each strain is aligned to it
once and called once (not once per contig).

```sh
pfb haplotype -n --cores 1 --configfile config/config.hg008.yaml   # dry run
pfb haplotype --cores <N> --configfile config/config.hg008.yaml
```

`config/config.hg008.yaml` is a template (graph, FASTQ paths are
placeholders). What's specific to this mode:
- **Graph**: `.gbz` (extracted with `vg paths`) or `.og` (`odgi paths`),
  picked by file extension. **`vg` must be on `PATH` for `.gbz` graphs.**
  **A `.gbz` must have at least one reference sample set**, otherwise
  `vg paths` refuses the graph ("graph does not contain any paths").
  Add one with:
  ```sh
  vg gbwt -Z in.gbz --set-reference <SAMPLE> [--set-reference ...] --gbz-format -g out.ref.gbz
  ```
  Haplotype paths are selected by name prefix (`vg paths -F -Q '<sample>#<hap>#'`),
  and FASTA headers are kept exactly as `vg` prints them, including the
  trailing phase-block field (e.g. `HG00423#1#JBHIIO010000010.1#0`).
  The extracted contig names and lengths are written to
  `<results_dir>/logs/extract_haplotype/<hap_id>.log`; the step fails if no
  path matches the prefix.
- **Paired-end reads**: a `strains:` entry can be `{r1: ..., r2: ...}` as
  well as a single FASTQ (this works in every mode). Every alignment gets a
  read group with `ID`/`SM` set to the strain name (`PL:ILLUMINA` too with
  `minimap2_preset: "sr"`).
- **Joint calling**: with `joint: true`, all strains' BAMs are merged and
  called together into one VCF with a column per sample,
  `<results_dir>/calls/joint/<hap_id>.vcf`. Without it, you get one VCF per
  strain, as in the other modes. `joint` can't be used as a strain name.

`deconstruct_vcf`, `min_size_diff` and `padding` are only needed for bubble
mode.

### What each stage does

The workflow is six Snakemake rules, each in its own file under
`workflow/rules/`. `discover_bubbles` → `extract_region` is one way to
reach a reference (bubble-local regions found by scanning a graph for
structurally divergent sites); `extract_path` is a second, independent way
(an entire named path from the graph, used whole, no discovery involved).
`extract_haplotype` is a third (all paths of one haplotype as one
multi-contig reference; see "Haplotype mode" above, which also covers the
optional joint-calling rules `merge_bams`/`call_variants_joint`). All of
them feed into the same `align_region`/`call_variants` steps.

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
- **`extract_path`** (`extract_path.smk`). Alternative to the two rules
  above: given one path name listed in `whole_paths` in the config, dumps
  FASTA for every path in the graph with `odgi paths -f` and keeps only the
  one record that matches (no bubble discovery, no padding, no subgraph
  extraction needed). Produces `<path-id>.fasta`. Reached via
  `pfb path ...` instead of `pfb bubble ...`, kept as
  a separate target since the two modes are different things you may or may
  not want in the same run.
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
