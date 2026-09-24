#set document(title: "panFreebayes filtering audit", author: "panFreebayes")
#set page(numbering: "1", margin: (x: 2.4cm, y: 2.6cm))
#set text(font: "New Computer Modern", size: 10.5pt)
#set par(justify: true, leading: 0.62em)
#show heading: set block(above: 1.4em, below: 0.8em)
#set heading(numbering: "1.1")
#show raw.where(block: true): set block(fill: luma(245), inset: 8pt, radius: 3pt, width: 100%)

#align(center)[
  #text(17pt, weight: "bold")[Filtering audit: what FreeBayes filters, and what the extraction keeps]
  #v(2pt)
  #text(10pt)[panFreebayes · answer to Pjotr's Milestone-1 question]
]

#outline(depth: 2, indent: auto)
#v(1em)

= The question and the method

Pjotr's question: _"Is there any existing filtering in FreeBayes that isn't
included in what we extracted? Anything that filters candidate variants or
reads that we may have implicitly dropped by scoping to 'engine only'."_

Method: every `if (…) { skip = true }`, `continue`, `pop_front`, `erase`,
`remove_if`, threshold comparison and `--flag`-gated branch in the
read → candidate-allele → emitted-record path was enumerated from
`freebayes.cpp`, `AlleleParser.cpp`, `Parameters.cpp`, `Allele.cpp`,
`Sample.cpp` and `BedReader.cpp`, and each was classified as *in the extraction*
or *driver-only*, cross-checked against the Milestone-0 "driver vs engine"
split and the actual extraction (`src/panfreebayes/panfreebayes_core.cpp`).

The relevant fact about the extraction: `callVariantsArgv` is a *verbatim* copy
of `freebayes.cpp`'s `main()` loop body (the only diffs are `argv` marshalling,
the output stream, and deleted comments — confirmed by `diff`), and it drives an
*unmodified* `AlleleParser`. So the question reduces to: _what filtering only
runs on a code path that `panfreebayes` cannot reach?_

= Answer

*No candidate-variant or read filter is lost.* Every read-level, base-level,
allele-level, candidate-level and site-level filter in FreeBayes runs
identically in panFreebayes, because the entire `main()` loop and the entire
`AlleleParser` are included unchanged.

The only things scoped out are `-r/--region`, `-t/--targets`, `-c/--stdin` and
`-L/--bam-list`. Of these, _only `-t` has any filtering character_, and:

- `-t` is *positional windowing*, not a candidate-variant filter — it makes the
  driver *visit* only positions inside the BED intervals (by jumping between
  them). `AlleleParser::inTarget()`, the one function that could apply a BED
  test inside the calling loop, *is never called* — it is dead code
  (`AlleleParser.cpp:866`, zero call sites).
- panFreebayes replaces it structurally: "one bubble = one reference FASTA + one
  BAM", analysed whole. There is no sub-region to exclude because the input
  *is* the region.

The single caveat (§4): if a future use wants to *exclude* a sub-region of a
bubble (a known-bad repeat, say), FreeBayes could do that with `-t`; panFreebayes
cannot, and the work-around is to pre-mask the reference or pre-filter the BAM.

= Complete filter inventory (all present in the extraction)

Every entry below runs in `panfreebayes call` exactly as in stock
single-process `freebayes`. Defaults from `Parameters.cpp`; ★ = overridden by
the panFreebayes acceptance command
(`--pooled-continuous --min-alternate-count 2 --min-alternate-fraction 0.2
--limit-coverage 200`).

== Read-level — `updateAlignmentQueue` (`AlleleParser.cpp:1954`)

#table(
  columns: (auto, auto, 1fr),
  stroke: 0.4pt + luma(200), inset: 5pt,
  [*filter*], [*param / default*], [*where*],
  [read-group not in analysed sample set], [—], [`:2010`],
  [duplicate read], [`--use-duplicate-reads` off ⇒ dropped], [`:2016`],
  [unmapped], [—], [`:2022`],
  [zero aligned bases], [—], [`:2028`],
  [secondary alignment], [`SecondaryFlag()`], [`:2034`],
  [coordinate out of order], [—], [`:2040`],
  [mapping quality], [`-m` / `MQL0` = 1], [`:2051`],
  [per-position coverage cap (deterministic)], [`-g` / `--skip-coverage` = 0 (off)], [`:2075`],
  [base-quality cap], [`--base-quality-cap` = 0 (off)], [`:2069`],
  [no alleles produced], [—], [`:2104`],
  [mismatch fraction], [`-z` / `--read-max-mismatch-fraction` = 1.0 (off)], [`:2105`],
  [absolute mismatch count], [`-U` / `--read-mismatch-limit` = 10⁷ (off)], [`:2106`],
  [SNP count], [`--read-snp-limit` = 10⁷ (off)], [`:2107`],
  [indel count], [`--read-indel-limit` = 10⁷ (off)], [`:2108`],
)

== Read-level — `--limit-coverage` random downsample ★

`freebayes.cpp:163–220` (in the lifted loop). Per site, per sample: if coverage
exceeds `--limit-coverage`, alleles are dropped at random (`rand()`, seeded
`srand(13)` at `:109`) until the cap is met. In the extraction. *Non-deterministic
across differently-chunked runs* — see the realignment doc §7 and the M0 report
§7. panFreebayes re-seeds per call and never chunks, so it reproduces
single-process FreeBayes.

== Base / allele-level — `registerAlignment`, `makeAllele`, `getAlleles`

#table(
  columns: (auto, auto, 1fr),
  stroke: 0.4pt + luma(200), inset: 5pt,
  [*filter*], [*param / default*], [*where*],
  [mismatch quality to count as a mismatch], [`-Q` / `BQL2` = 10], [`AlleleParser.cpp:1491`],
  [read base not `ACGT` ⇒ `ALLELE_NULL` (uncallable)], [—], [`:1465`, `:1515`],
  [reference base `N` ⇒ position skipped], [—], [`freebayes.cpp:142`],
  [terminal deletion (first/last CIGAR op) ignored], [—], [`AlleleParser.cpp:1681`],
  [soft clip ⇒ `ALLELE_NULL`], [—], [`:1783`],
  [`--trim-complex-tail` tail split], [0 (off)], [`:1840`],
  [haplotype-basis-allele gating ⇒ demote to reference], [`--haplotype-basis-alleles` (off)], [`:1222`],
  [allele quality], [`-q` / `BQL0` = 0], [`getAlleles`, `:3739`],
  [allele `currentBase == "N"` or empty alt], [—], [`:3739`],
  [allele-type mask], [`--no-snps` / `--no-indels` / `--no-mnps` / `--no-complex` (all on)], [`freebayes.cpp:88`],
)

== Candidate-allele-level — `genotypeAlleles` (`AlleleParser.cpp:3843`) + main loop

#table(
  columns: (auto, auto, 1fr),
  stroke: 0.4pt + luma(200), inset: 5pt,
  [*filter*], [*param / default*], [*where*],
  [group alt not all-`ACGT`], [—], [`:3860`],
  [group size], [`-G` / `--min-alternate-total` / `minAltTotal` = 1], [`:3864`],
  [supporting base-quality sum], [`-3` params / `minSupportingAlleleQualitySum`], [`:3877`],
  [supporting mapping-quality sum], [`minSupportingMappingQualitySum`], [`:3877`],
  [per-sample alt quality sum], [`--min-alternate-qsum` / `minAltQSum` = 0], [`:3927`],
  [per-sample alt count ★], [`-C` / `--min-alternate-count` = 2], [`:3928`],
  [per-sample alt fraction ★], [`-F` / `--min-alternate-fraction` = 0.05], [`:3929`],
  [keep only N best alleles], [`-n` / `--use-best-n-alleles` = 0 (all)], [`:3949`+],
  [`sufficientAlternateObservations` pre-gate], [uses `-C` and `-F`], [`freebayes.cpp:227`],
  [`--min-coverage`], [`-!` / `minCoverage` = 0], [`freebayes.cpp:157`],
  [coverage == 0 ⇒ skip], [—], [`freebayes.cpp:154`],
  [only the reference allele survives ⇒ skip], [—], [`freebayes.cpp:312`],
)

== Site / output-level — `freebayes.cpp`

#table(
  columns: (auto, auto, 1fr),
  stroke: 0.4pt + luma(200), inset: 5pt,
  [*filter*], [*param / default*], [*where*],
  [posterior probability of variation], [`-P` / `--pvar` / `PVL` = 0.0], [`:645`],
  [emit only if `alts` non-empty (unless `PVL == 0`)], [—], [`:645`],
  [`--report-monomorphic`], [off], [`:226`, `:231`],
  [`--report-genotype-likelihood-max`], [off], [`:603`],
  [no sample data likelihoods ⇒ skip], [—], [`:393`],
)

With `--pooled-continuous` and the default `PVL = 0.0`, the site-level gate
reduces to "emit if there is at least one alt allele that cleared the
candidate-level filters" (`freebayes.cpp:574–582`, `:645`).

= Driver-only code, and whether it filters

#table(
  columns: (auto, 1fr, auto),
  stroke: 0.4pt + luma(200), inset: 5pt,
  [*driver feature*], [*what it does*], [*filter gap?*],
  [`-r` / `--region`], [analyse only `chr:S-E`; `loadTargets` (`:722`), window clip in `toNextPosition` (`:2901`)], [No — positional window; panFreebayes uses `[0, contigLen)`],
  [`-t` / `--targets` (BED)], [analyse only positions inside BED intervals, by *jumping* between them (`toNextTarget`); interval tree built at `:792`], [*Partial* — see below],
  [`AlleleParser::inTarget()`], [would test "current position inside a BED target"], [*Dead code* — 0 call sites (`:866`)],
  [BAM-index region seek (`SetRegion`)], [`loadTarget` (`:2732`) — only reached when `targets` non-empty], [No — I/O],
  [`-c` / `--stdin`, `-L` / `--bam-list`], [alternative BAM input plumbing], [No — I/O],
  [`freebayes-parallel` wrapper], [`fasta_generate_regions.py` split + `parallel` + `vcffirstheader | vcfstreamsort | vcfuniq`], [No — the `vcfuniq` step only removes *duplicate* records at chunk boundaries; with a single un-chunked region there are none. See §5],
)

== The one partial gap: `-t` BED as an exclusion mechanism

`-t` is normally used to *include* regions of interest. Because the driver walks
BED intervals and clips `[left, right]` on each (`toNextPosition:2901`,
`loadTarget:2722`), a BED file can also be used to *exclude* a sub-region by
simply not covering it.

panFreebayes has no equivalent: `panfreebayes call` rejects `-t` outright
(`panfreebayes_core.cpp`, `isRegionOrTargetFlag`) and analyses the whole
reference. This is by design — the tool's unit of work is a single bubble whose
reference FASTA *is* the region — but it does mean:

- if a bubble's local reference contains a stretch you want FreeBayes to *not*
  call in (e.g. a collapsed repeat that generates spurious calls, like the
  ~23,000× spike in the DL238 bubble), the FreeBayes way (`-t` a BED that
  excludes it) is unavailable;
- the panFreebayes way is upstream: hard-mask that stretch to `N` in the local
  reference FASTA (`N` reference bases are skipped, `freebayes.cpp:142`), or
  drop those reads from the local BAM before calling.

This is a *usability* difference, not a correctness one: nothing is being
silently filtered *differently*; a capability is simply not exposed.

= `freebayes-parallel` and post-hoc filtering

The manual pipeline's historical numbers may have come through
`scripts/freebayes-parallel`, whose output pipeline is:

```
… | vcffirstheader | vcfstreamsort -w 1000 | vcfuniq
```

None of these is a variant *filter* in the sense of "drops calls that pass
FreeBayes' thresholds":

- `vcffirstheader` — keeps one header;
- `vcfstreamsort` — re-sorts records within a 1000-record window (chunk
  boundaries interleave);
- `vcfuniq` — removes records that are *byte-identical to the previous line*,
  which only happens for a variant reported by two adjacent chunks that both
  overlap it.

A single-region panFreebayes run produces no chunk-boundary duplicates, so
there is nothing for `vcfuniq` to remove. The realignment doc §7 / M0 §7 point
stands separately: chunked and un-chunked runs can *genuinely* disagree on
`--limit-coverage` sites because of the shared PRNG — that is a divergence in
*which reads survive*, upstream of any of these post-filters, and panFreebayes
matches the *single-process* result.

= Recommendation

+ Record in the milestone report that the extraction preserves 100% of
  FreeBayes' read/base/allele/candidate/site filtering, and that the only
  scoped-out capability is `-r`/`-t` positional restriction, intentionally
  replaced by the per-bubble reference model.
+ If sub-bubble exclusion is ever needed, do it by masking the local reference
  (`N`) or pre-filtering the local BAM — not by adding region flags back to the
  core (which would re-introduce the BAM-index dependency the guard exists to
  prevent).
+ `AlleleParser::inTarget()` is dead code in upstream FreeBayes; worth a
  one-line note if the maintainer contact happens (alongside the two
  non-determinism findings).
