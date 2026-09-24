#ifndef PANFREEBAYES_CORE_H
#define PANFREEBAYES_CORE_H

//
// panFreebayes core engine
// ------------------------
//
// A standalone entry point into FreeBayes's realignment + variant-calling
// engine, callable without the `freebayes` CLI binary.
//
// Milestone 1 scope (see panfreebayes_instructions.md):
//   - takes a local reference FASTA + a BAM aligned to it + calling parameters
//   - runs the EXACT per-position engine loop from src/freebayes.cpp
//     (realignment -> allele registration -> haplotype window -> data
//     likelihoods -> Bayesian combo search -> VCF record)
//   - no region (-r) / target (-t) selection, no BAM-index jumping:
//     the whole reference is swept, exactly as `freebayes -f ref.fa aln.bam`
//     with no region flags. This is enforced, not just the default: passing
//     -r/-t/--stdin/-L (directly or via Options::extraArgs) is rejected, so
//     the index-seeking code path can never be reached silently.
//
// Behaviour is intended to be BIT-FOR-BIT identical to running the stock
// `freebayes` binary with the equivalent arguments. In particular the global
// PRNG used by --limit-coverage (srand(13) once, then rand()) is seeded and
// consumed in the same order. See panfreebayes_milestone0_orientation.md sec. 7.
//

#include <string>
#include <vector>
#include <iosfwd>

namespace panfreebayes {

// High-level options. Any field left at its sentinel is simply not passed,
// so FreeBayes's own default applies. `extraArgs` is a raw passthrough for
// any freebayes flag not modelled here (inserted before the BAM argument).
struct Options {
    std::string fasta;                    // required: local reference FASTA
    std::string bam;                      // required: BAM aligned to `fasta`

    // The tuning flags needed empirically for noisy long-read data
    // (sentinel = "not set", use freebayes default):
    bool   pooledContinuous     = false;  // --pooled-continuous
    int    minAlternateCount    = -1;     // --min-alternate-count  (>=0 to set)
    double minAlternateFraction = -1.0;   // --min-alternate-fraction (>=0 to set)
    int    limitCoverage        = -1;     // --limit-coverage (>=0 to set)

    std::vector<std::string> extraArgs;   // raw freebayes flags, verbatim
};

// Build the argv vector that `callVariantsArgv` would run for `opt`
// (argv[0] == "freebayes"). Exposed for logging / the regression harness.
std::vector<std::string> buildArgv(const Options& opt);

// Run the engine for `opt`, writing a VCF stream to `out`.
// Returns 0 on success. NOTE: on fatal input errors the underlying
// FreeBayes code calls exit() (Milestone 1 limitation, documented).
int callVariants(const Options& opt, std::ostream& out);

// argv-level entry point. `argv` is a full freebayes command line with
// argv[0] == "freebayes" (or anything). This is what the CLI uses, and is the
// single place the src/freebayes.cpp main-loop body is reproduced.
int callVariantsArgv(const std::vector<std::string>& argv, std::ostream& out);

// Version string for the CLI / logging (single source of truth).
std::string version();

} // namespace panfreebayes

#endif
