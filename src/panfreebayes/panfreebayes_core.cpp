//
// panFreebayes core engine — implementation
//
// callVariantsArgv() below is a VERBATIM lift of the main() loop body from
// src/freebayes.cpp (as of the build-fix branch). The only differences:
//   - argc/argv come from a std::vector<std::string> instead of the process
//   - the output stream is redirected to the caller's std::ostream
//   - no signal() / segfault handler install (that stays in the CLI wrapper)
//   - returns instead of exit(0), and does not delete-then-continue past errors
//
// Keep this in lock-step with src/freebayes.cpp until the two are refactored
// onto a shared helper (tracked as tech debt for Milestone 1 follow-up).
//

#include "panfreebayes_core.h"

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <list>
#include <algorithm>
#include <iterator>
#include <cmath>
#include <cassert>
#include <cstdlib>

#include "TryCatch.h"
#include "Parameters.h"
#include "Allele.h"
#include "Sample.h"
#include "AlleleParser.h"
#include "Utility.h"

#include <vcflib/multichoose.h>

#include "Genotype.h"
#include "DataLikelihood.h"
#include "Marginals.h"
#include "ResultData.h"

#include "Bias.h"
#include "Contamination.h"
#include "NonCall.h"
#include "Logging.h"
#include "version_git.h"
#include <vcflib/Variant.h>

using namespace std;

namespace panfreebayes {

std::string version() {
    return std::string("panfreebayes (freebayes engine ") + VERSION_GIT + ")";
}

// ---------------------------------------------------------------------------
// Milestone 1, "Step 2" (decision 2).
//
// Decision 2 asked for reads to be pre-registered rather than streamed through
// a bounded queue, with the stated goal of removing the BAM-index /
// region-iteration dependency. Investigating that change surfaced that:
//
//   - AlleleParser::updateAlignmentQueue() only ever pulls reads with
//     POSITION <= currentPosition (AlleleParser.cpp:1978, 2118-2120), so
//     buildHaplotypeAlleles' window-growth loop (AlleleParser.cpp:3266-3280)
//     can only see reads that already ended inside the growing window, never
//     ones that merely *start* inside it. Pre-registering all reads gives it
//     those too, which can grow haplotypeLength further, change `theta`
//     (freebayes.cpp:309), and change which sites get absorbed into a
//     haplotype block -- i.e. it can change the call set. That would mean
//     panfreebayes could no longer reproduce freebayes' 1,413 / 2,451
//     acceptance numbers, and "provably identical before/after" would not
//     hold (flagged and resolved with the supervisor 2026-09-04/05).
//
//   - Separately, and this is what actually closes decision 2's real goal:
//     the SeqLib BAM backend this build uses (HAVE_BAMTOOLS undefined) never
//     calls LocateIndexes() at all (AlleleParser.cpp's non-BAMTOOLS
//     openBams() branch has that call commented out), and
//     AlleleParser::loadTarget()'s bamMultiReader.SetRegion() -- the only
//     index-dependent seek in the engine -- is reachable only when
//     `parameters.targets` is non-empty, i.e. only when -r/--region or
//     -t/--targets is passed. buildArgv() above never emits either flag.
//     So panfreebayes already runs FreeBayes's no-targets, whole-BAM,
//     index-free streaming path *by construction* -- with byte-identical
//     behaviour to stock freebayes, confirmed by the Step 1 smoke test.
//
// So "Step 2" is not an engine change: it is making that already-true
// invariant impossible to break by accident (a caller smuggling -r/-t/--stdin
// into extraArgs would silently re-enable the index-seeking path). Reads
// still stream positionally through the exact unmodified AlleleParser code.
// ---------------------------------------------------------------------------
namespace {
bool isRegionOrTargetFlag(const std::string& tok) {
    static const std::vector<std::string> disallowed = {
        "-r", "--region", "-t", "--targets", "-c", "--stdin", "-L", "--bam-list",
    };
    return std::find(disallowed.begin(), disallowed.end(), tok) != disallowed.end();
}
} // namespace

std::vector<std::string> buildArgv(const Options& opt) {
    std::vector<std::string> a;
    a.push_back("freebayes");
    a.push_back("-f");
    a.push_back(opt.fasta);
    if (opt.pooledContinuous) {
        a.push_back("--pooled-continuous");
    }
    if (opt.minAlternateCount >= 0) {
        a.push_back("--min-alternate-count");
        a.push_back(std::to_string(opt.minAlternateCount));
    }
    if (opt.minAlternateFraction >= 0.0) {
        std::ostringstream s;
        s << opt.minAlternateFraction;
        a.push_back("--min-alternate-fraction");
        a.push_back(s.str());
    }
    if (opt.limitCoverage >= 0) {
        a.push_back("--limit-coverage");
        a.push_back(std::to_string(opt.limitCoverage));
    }
    for (const auto& e : opt.extraArgs) {
        a.push_back(e);
    }
    a.push_back(opt.bam);
    return a;
}

int callVariants(const Options& opt, std::ostream& out) {
    return callVariantsArgv(buildArgv(opt), out);
}

// ---------------------------------------------------------------------------
// The lifted engine loop.
// ---------------------------------------------------------------------------
int callVariantsArgv(const std::vector<std::string>& argvStrings, std::ostream& outStream) {

    // Guarantee the no-targets, index-free, whole-BAM streaming path (see the
    // comment above isRegionOrTargetFlag()): reject -r/-t/--stdin/etc. rather
    // than silently falling onto AlleleParser's index-seeking SetRegion() path.
    for (const auto& tok : argvStrings) {
        if (isRegionOrTargetFlag(tok)) {
            std::cerr << "panfreebayes: '" << tok << "' is not supported -- "
                       << "panfreebayes always analyses the entire supplied "
                       << "reference/BAM as one region (no -r/-t/--stdin/-L). "
                       << "Extract the region into its own reference FASTA + BAM "
                       << "upstream instead." << std::endl;
            return 1;
        }
    }

    // marshal argv
    std::vector<char*> argv;
    argv.reserve(argvStrings.size() + 1);
    for (const auto& s : argvStrings) {
        argv.push_back(const_cast<char*>(s.c_str()));
    }
    argv.push_back(nullptr);
    int argc = static_cast<int>(argvStrings.size());

    AlleleParser* parser = new AlleleParser(argc, argv.data());
    parser->output = &outStream;   // redirect VCF stream to the caller

    Parameters& parameters = parser->parameters;
    list<Allele*> alleles;

    Samples samples;
    NonCalls nonCalls;

    ostream& out = *(parser->output);

    Bias observationBias;
    if (!parameters.alleleObservationBiasFile.empty()) {
        observationBias.open(parameters.alleleObservationBiasFile);
    }

    Contamination contaminationEstimates(0.5+parameters.probContamination, parameters.probContamination);
    if (!parameters.contaminationEstimateFile.empty()) {
        contaminationEstimates.open(parameters.contaminationEstimateFile);
    }

    // this can be uncommented to force operation on a specific set of genotypes
    vector<Allele> allGenotypeAlleles;
    allGenotypeAlleles.push_back(genotypeAllele(ALLELE_GENOTYPE, "A", 1));
    allGenotypeAlleles.push_back(genotypeAllele(ALLELE_GENOTYPE, "T", 1));
    allGenotypeAlleles.push_back(genotypeAllele(ALLELE_GENOTYPE, "G", 1));
    allGenotypeAlleles.push_back(genotypeAllele(ALLELE_GENOTYPE, "C", 1));

    int allowedAlleleTypes = ALLELE_REFERENCE;
    if (parameters.allowSNPs) {
        allowedAlleleTypes |= ALLELE_SNP;
    }
    if (parameters.allowIndels) {
        allowedAlleleTypes |= ALLELE_INSERTION;
        allowedAlleleTypes |= ALLELE_DELETION;
    }
    if (parameters.allowMNPs) {
        allowedAlleleTypes |= ALLELE_MNP;
    }
    if (parameters.allowComplex) {
        allowedAlleleTypes |= ALLELE_COMPLEX;
    }

    // output VCF header
    if (parameters.output == "vcf") {
        out << parser->variantCallFile.header << endl;
    }

    if (0 < parameters.limitCoverage) {
        srand(13);
    }

    Allele nullAllele = genotypeAllele(ALLELE_NULL, "N", 1, "1N");

    unsigned long total_sites = 0;
    unsigned long processed_sites = 0;

    while (parser->getNextAlleles(samples, allowedAlleleTypes)) {

        ++total_sites;

        DEBUG2("at start of main loop");

        // did we switch chromosomes or exceed our gVCF chunk size, or do we not want to use chunks?
        // if so, we may need to output a gVCF record
        Results results;
        if (parameters.gVCFout
               &&  !(nonCalls.empty())
               &&  (  (parameters.gVCFNoChunk)
                   || (nonCalls.begin()->first != parser->currentSequenceName)
                   || (parameters.gVCFchunk
                       && nonCalls.lastPos().second - nonCalls.firstPos().second >= parameters.gVCFchunk
                      )
                  )
            ){
            vcflib::Variant var(parser->variantCallFile);
            out << results.gvcf(var, nonCalls, parser) << endl;
            nonCalls.clear();
        }

        // don't process non-ATGCN's in the reference
        string cb = parser->currentReferenceBaseString();
        if (cb != "A" && cb != "T" && cb != "C" && cb != "G" && cb != "N") {
            DEBUG2("current reference base is not in { A T G C N }");
            continue;
        }

        int coverage = countAlleles(samples);

        DEBUG("position: " << parser->currentSequenceName << ":" << (long unsigned int) parser->currentPosition + 1 << " coverage: " << coverage);

        bool skip = false;
        if (!parser->hasInputVariantAllelesAtCurrentPosition()) {
            // skips 0-coverage regions
            if (coverage == 0) {
                DEBUG("no alleles left at this site after filtering");
                skip = true;
            } else if (coverage < parameters.minCoverage) {
                DEBUG("post-filtering coverage of " << coverage << " is less than --min-coverage of " << parameters.minCoverage);
                skip = true;
            } else if (parameters.onlyUseInputAlleles) {
                DEBUG("no input alleles, but using only input alleles for analysis, skipping position");
                skip = true;
            } else if (0 < parameters.limitCoverage) {
                // go through each sample
                for (Samples::iterator s = samples.begin(); s != samples.end(); ++s) {
                    string sampleName = s->first;
                    Sample& sample = s->second;
                    // get the coverage for this sample
                    int sampleCoverage = 0;
                    for (Sample::iterator sg = sample.begin(); sg != sample.end(); ++sg) {
                        sampleCoverage += sg->second.size();
                    }
                    if (sampleCoverage <= parameters.limitCoverage) {
                        continue;
                    }

                    DEBUG("coverage " << sampleCoverage << " for sample " << sampleName << " was > " << parameters.limitCoverage << ", so we will remove " << (sampleCoverage - parameters.limitCoverage) << " genotypes");
                    vector<string> genotypesToErase;
                    do {
                        double probRemove = (sampleCoverage - parameters.limitCoverage) / (double)sampleCoverage;
                        vector<string> genotypesToErase;
                        // iterate through the genotypes
                        for (Sample::iterator sg = sample.begin(); sg != sample.end(); ++sg) {
                            vector<Allele*> allelesToKeep;
                            // iterate through each allele
                            for (int alleleIndex = 0; alleleIndex < sg->second.size(); alleleIndex++) {
                                // only if we have more alleles to remove
                                if (parameters.limitCoverage < sampleCoverage) {
                                    double r = rand() / (double)RAND_MAX;
                                    if (r < probRemove) { // skip over this allele
                                        sampleCoverage--;
                                        continue;
                                    }
                                }
                                // keep it
                                allelesToKeep.push_back(sg->second[alleleIndex]);
                            }
                            // re-assign the alleles to this genotype
                            if (allelesToKeep.size() < sg->second.size()) {
                                sg->second.assign(allelesToKeep.begin(), allelesToKeep.end());
                            }
                            // if no more alleles for this genotype, remove it later
                            if (sg->second.empty()) {
                                genotypesToErase.push_back(sg->first);
                            }
                        }
                        // remove empty genotypes
                        for (vector<string>::iterator gt = genotypesToErase.begin(); gt != genotypesToErase.end(); ++gt) {
                            sample.erase(*gt);
                        }
                    } while (parameters.limitCoverage < sampleCoverage);
                    sampleCoverage = 0;
                    for (Sample::iterator sg = sample.begin(); sg != sample.end(); ++sg) {
                        sampleCoverage += sg->second.size();
                    }
                    DEBUG("coverage for sample " << sampleName << " is now " << sampleCoverage);
                }
                // update coverage
                coverage = countAlleles(samples);
            }

            DEBUG2("coverage " << parser->currentSequenceName << ":" << parser->currentPosition << " == " << coverage);

            // establish a set of possible alternate alleles to evaluate at this location

            if (!parameters.reportMonomorphic
                && !sufficientAlternateObservations(samples, parameters.minAltCount, parameters.minAltFraction)) {
                DEBUG("insufficient alternate observations");
                skip = true;
            }
            if (parameters.reportMonomorphic) {
                DEBUG("calling at site even though there are no alternate observations");
            }
        }

        if (skip) {
            // record data for gVCF
            if (parameters.gVCFout) {
                nonCalls.record(parser->currentSequenceName, parser->currentPosition, samples);
            }
            // and step ahead
            continue;
        }

        // to ensure proper ordering of output stream
        vector<string> sampleListPlusRef;

        for (vector<string>::iterator s = parser->sampleList.begin(); s != parser->sampleList.end(); ++s) {
            sampleListPlusRef.push_back(*s);
        }
        if (parameters.useRefAllele) {
            sampleListPlusRef.push_back(parser->currentSequenceName);
        }

        // establish genotype alleles using input filters
        map<string, vector<Allele*> > alleleGroups;
        groupAlleles(samples, alleleGroups);
        DEBUG2("grouped alleles by equivalence");

        vector<Allele> genotypeAlleles = parser->genotypeAlleles(alleleGroups, samples, parameters.onlyUseInputAlleles);

        // always include the reference allele as a possible genotype, even when we don't include it by default
        if (!parameters.useRefAllele) {
            vector<Allele> refAlleleVector;
            refAlleleVector.push_back(genotypeAllele(ALLELE_REFERENCE, string(1, parser->currentReferenceBase), 1, "1M"));
            genotypeAlleles = alleleUnion(genotypeAlleles, refAlleleVector);
        }

        map<string, vector<Allele*> > partialObservationGroups;
        map<Allele*, set<Allele*> > partialObservationSupport;

        // build haplotype alleles matching the current longest allele (often will do nothing)
        // this will adjust genotypeAlleles if changes are made
        DEBUG("building haplotype alleles, currently there are " << genotypeAlleles.size() << " genotype alleles");
        DEBUG(genotypeAlleles);
        parser->buildHaplotypeAlleles(genotypeAlleles,
                                      samples,
                                      alleleGroups,
                                      partialObservationGroups,
                                      partialObservationSupport,
                                      allowedAlleleTypes);
        DEBUG("built haplotype alleles, now there are " << genotypeAlleles.size() << " genotype alleles");
        DEBUG(genotypeAlleles);

        string referenceBase = parser->currentReferenceHaplotype();

        // re-calculate coverage, as this could change now that we've built haplotype alleles
        coverage = countAlleles(samples);

        // estimate theta using the haplotype length
        long double theta = parameters.TH * parser->lastHaplotypeLength;

        // if we have only one viable allele, we don't have evidence for variation at this site
        if (!parser->hasInputVariantAllelesAtCurrentPosition() && !parameters.reportMonomorphic && genotypeAlleles.size() <= 1 && genotypeAlleles.front().isReference()) {
            DEBUG("no alternate genotype alleles passed filters at " << parser->currentSequenceName << ":" << parser->currentPosition);
            continue;
        }
        DEBUG("genotype alleles: " << genotypeAlleles);

        // add the null genotype
        bool usingNull = false;
        if (parameters.excludeUnobservedGenotypes && genotypeAlleles.size() > 2) {
            genotypeAlleles.push_back(nullAllele);
            usingNull = true;
        }

        ++processed_sites;

        // generate possible genotypes

        // for each possible ploidy in the dataset, generate all possible genotypes
        vector<int> ploidies = parser->currentPloidies(samples);
        map<int, vector<Genotype> > genotypesByPloidy = getGenotypesByPloidy(ploidies, genotypeAlleles);
        int numCopiesOfLocus = parser->copiesOfLocus(samples);

        DEBUG2("generated all possible genotypes:");
        if (parameters.debug2) {
            for (map<int, vector<Genotype> >::iterator s = genotypesByPloidy.begin(); s != genotypesByPloidy.end(); ++s) {
                vector<Genotype>& genotypes = s->second;
                for (vector<Genotype>::iterator g = genotypes.begin(); g != genotypes.end(); ++g) {
                    DEBUG2(*g);
                }
            }
        }

        // get estimated allele frequencies using sum of estimated qualities
        map<string, double> estimatedAlleleFrequencies = samples.estimatedAlleleFrequencies();
        double estimatedMaxAlleleFrequency = 0;
        double estimatedMaxAlleleCount = 0;
        double estimatedMajorFrequency = estimatedAlleleFrequencies[referenceBase];
        if (estimatedMajorFrequency < 0.5) estimatedMajorFrequency = 1-estimatedMajorFrequency;
        double estimatedMinorFrequency = 1-estimatedMajorFrequency;
        int estimatedMinorAllelesAtLocus = max(1, (int) ceil((double) numCopiesOfLocus * estimatedMinorFrequency));

        map<string, vector<vector<SampleDataLikelihood> > > sampleDataLikelihoodsByPopulation;
        map<string, vector<vector<SampleDataLikelihood> > > variantSampleDataLikelihoodsByPopulation;
        map<string, vector<vector<SampleDataLikelihood> > > invariantSampleDataLikelihoodsByPopulation;

        map<string, int> inputAlleleCounts;
        int inputLikelihoodCount = 0;

        DEBUG2("calculating data likelihoods");
        calculateSampleDataLikelihoods(
            samples,
            results,
            parser,
            genotypesByPloidy,
            parameters,
            usingNull,
            observationBias,
            genotypeAlleles,
            contaminationEstimates,
            estimatedAlleleFrequencies,
            sampleDataLikelihoodsByPopulation,
            variantSampleDataLikelihoodsByPopulation,
            invariantSampleDataLikelihoodsByPopulation);

        DEBUG2("finished calculating data likelihoods");

        // if somehow we get here without any possible sample genotype likelihoods, bail out
        bool hasSampleLikelihoods = false;
        for (map<string, vector<vector<SampleDataLikelihood> > >::iterator s = sampleDataLikelihoodsByPopulation.begin();
             s != sampleDataLikelihoodsByPopulation.end(); ++s) {
            if (!s->second.empty()) {
                hasSampleLikelihoods = true;
                break;
            }
        }
        if (!hasSampleLikelihoods) {
            continue;
        }

        DEBUG2("calulating combo posteriors over " << parser->populationSamples.size() << " populations");

        BigFloat pVar = 1.0;
        BigFloat pHom = 0.0;

        long double bestComboOddsRatio = 0;

        bool bestOverallComboIsHet = false;
        GenotypeCombo bestCombo;

        GenotypeCombo bestGenotypeComboByMarginals;
        vector<vector<SampleDataLikelihood> > allSampleDataLikelihoods;

        DEBUG("searching genotype space");

        map<string, list<GenotypeCombo> > genotypeCombosByPopulation;
        int genotypingTotalIterations = 0; // tally total iterations required to reach convergence
        map<string, list<GenotypeCombo> > glMaxCombos;

        for (map<string, SampleDataLikelihoods>::iterator p = sampleDataLikelihoodsByPopulation.begin(); p != sampleDataLikelihoodsByPopulation.end(); ++p) {

            const string& population = p->first;
            SampleDataLikelihoods& sampleDataLikelihoods = p->second;
            list<GenotypeCombo>& populationGenotypeCombos = genotypeCombosByPopulation[population];

            DEBUG2("genqerating banded genotype combinations from " << sampleDataLikelihoods.size() << " sample genotypes in population " << population);

            // cap the number of iterations at 2 x the number of alternate alleles
            // max it at parameters.genotypingMaxIterations iterations, min at 10
            int itermax = min(max(10, 2 * estimatedMinorAllelesAtLocus), parameters.genotypingMaxIterations);

            // XXX HACK
            // passing 0 for bandwidth and banddepth means "exhaustive local search"
            // this produces properly normalized GQ's at polyallelic sites
            int adjustedBandwidth = 0;
            int adjustedBanddepth = 0;
            // however, this can lead to huge performance problems at complex sites,
            // so we implement this hack...
            if (parameters.genotypingMaxBandDepth > 0 &&
                genotypeAlleles.size() > parameters.genotypingMaxBandDepth) {
                adjustedBandwidth = 1;
                adjustedBanddepth = parameters.genotypingMaxBandDepth;
            }

            GenotypeCombo nullCombo;
            SampleDataLikelihoods nullSampleDataLikelihoods;

            // this is the genotype-likelihood maximum
            if (parameters.reportGenotypeLikelihoodMax) {
                GenotypeCombo comboKing;
                vector<int> initialPosition;
                initialPosition.assign(sampleDataLikelihoods.size(), 0);
                SampleDataLikelihoods nullDataLikelihoods; // dummy variable
                makeComboByDatalLikelihoodRank(comboKing,
                                               initialPosition,
                                               sampleDataLikelihoods,
                                               nullDataLikelihoods,
                                               inputAlleleCounts,
                                               theta,
                                               parameters.pooledDiscrete,
                                               parameters.ewensPriors,
                                               parameters.permute,
                                               parameters.hwePriors,
                                               parameters.obsBinomialPriors,
                                               parameters.alleleBalancePriors,
                                               parameters.diffusionPriorScalar);

                glMaxCombos[population].push_back(comboKing);
            }

            // search much longer for convergence
            convergentGenotypeComboSearch(
                populationGenotypeCombos,
                nullCombo,
                sampleDataLikelihoods, // vary everything
                sampleDataLikelihoods,
                nullSampleDataLikelihoods,
                samples,
                genotypeAlleles,
                inputAlleleCounts,
                adjustedBandwidth,
                adjustedBanddepth,
                theta,
                parameters.pooledDiscrete,
                parameters.ewensPriors,
                parameters.permute,
                parameters.hwePriors,
                parameters.obsBinomialPriors,
                parameters.alleleBalancePriors,
                parameters.diffusionPriorScalar,
                itermax,
                genotypingTotalIterations,
                true); // add homozygous combos
        }

        // generate the GL max combo
        GenotypeCombo glMax;
        if (parameters.reportGenotypeLikelihoodMax) {
            list<GenotypeCombo> glMaxGenotypeCombos;
            combinePopulationCombos(glMaxGenotypeCombos, glMaxCombos);
            glMax = glMaxGenotypeCombos.front();
        }

        // accumulate combos from independently-calculated populations into the list of combos
        list<GenotypeCombo> genotypeCombos; // build new combos into this list
        combinePopulationCombos(genotypeCombos, genotypeCombosByPopulation);

        // re-get posterior normalizer
        vector<long double> comboProbs;
        for (list<GenotypeCombo>::iterator gc = genotypeCombos.begin(); gc != genotypeCombos.end(); ++gc) {
            comboProbs.push_back(gc->posteriorProb);
        }
        long double posteriorNormalizer = logsumexp_probs(comboProbs);

        // recalculate posterior normalizer
        pVar = 1.0;
        pHom = 0.0;
        // calculates pvar and gets the best het combo
        list<GenotypeCombo>::iterator gc = genotypeCombos.begin();
        bestCombo = *gc;
        for ( ; gc != genotypeCombos.end(); ++gc) {
            if (gc->isHomozygous() && gc->alleles().front() == referenceBase) {
                pVar -= big_exp(gc->posteriorProb - posteriorNormalizer);
                pHom += big_exp(gc->posteriorProb - posteriorNormalizer);
            } else if (gc == genotypeCombos.begin()) {
                bestOverallComboIsHet = true;
            }
        }

        // odds ratio between the first and second-best combinations
        if (genotypeCombos.size() > 1) {
            bestComboOddsRatio = genotypeCombos.front().posteriorProb - (++genotypeCombos.begin())->posteriorProb;
        }

        if (parameters.calculateMarginals) {
            // make a combined, all-populations sample data likelihoods vector to accumulate marginals
            SampleDataLikelihoods allSampleDataLikelihoods;
            for (map<string, SampleDataLikelihoods>::iterator p = sampleDataLikelihoodsByPopulation.begin(); p != sampleDataLikelihoodsByPopulation.end(); ++p) {
                SampleDataLikelihoods& sdls = p->second;
                allSampleDataLikelihoods.reserve(allSampleDataLikelihoods.size() + distance(sdls.begin(), sdls.end()));
                allSampleDataLikelihoods.insert(allSampleDataLikelihoods.end(), sdls.begin(), sdls.end());
            }
            // calculate the marginal likelihoods for this population
            marginalGenotypeLikelihoods(genotypeCombos, allSampleDataLikelihoods);
            // store the marginal data likelihoods in the results, for easy parsing
            results.update(allSampleDataLikelihoods);
        }

        map<string, int> repeats;
        if (parameters.showReferenceRepeats) {
            repeats = parser->repeatCounts(parser->currentSequencePosition(), parser->currentSequence, 12);
        }

        vector<Allele> alts;
        if (parameters.onlyUseInputAlleles
            || parameters.reportAllHaplotypeAlleles
            || parameters.pooledContinuous) {
            for (vector<Allele>::iterator a = genotypeAlleles.begin(); a != genotypeAlleles.end(); ++a) {
                if (!a->isReference()) {
                    alts.push_back(*a);
                }
            }
        } else {
            // get the unique alternate alleles in this combo, sorted by frequency in the combo
            vector<pair<Allele, int> > alternates = alternateAlleles(bestCombo, referenceBase);
            for (vector<pair<Allele, int> >::iterator a = alternates.begin(); a != alternates.end(); ++a) {
                Allele& alt = a->first;
                if (!alt.isNull() && !alt.isReference())
                    alts.push_back(alt);
            }
            // if there are no alternate alleles in the best combo, use the genotype alleles
            if (alts.empty()) {
                for (vector<Allele>::iterator a = genotypeAlleles.begin(); a != genotypeAlleles.end(); ++a) {
                    if (!a->isReference()) {
                        alts.push_back(*a);
                    }
                }
            }
        }

        // reporting the GL maximum *over all alleles*
        if (parameters.reportGenotypeLikelihoodMax) {
            bestCombo = glMax;
        } else {
            // the default behavior is to report the GL maximum genotyping over the alleles in the best posterior genotyping
            vector<Allele> alleles = alts;
            for (vector<Allele>::iterator a = genotypeAlleles.begin(); a != genotypeAlleles.end(); ++a) {
                if (a->isReference()) {
                    alleles.push_back(*a);
                }
            }
            map<string, list<GenotypeCombo> > glMaxComboBasedOnAltsByPop;
            for (map<string, SampleDataLikelihoods>::iterator p = sampleDataLikelihoodsByPopulation.begin(); p != sampleDataLikelihoodsByPopulation.end(); ++p) {
                const string& population = p->first;
                SampleDataLikelihoods& sampleDataLikelihoods = p->second;
                GenotypeCombo glMaxBasedOnAlts;
                for (SampleDataLikelihoods::iterator v = sampleDataLikelihoods.begin(); v != sampleDataLikelihoods.end(); ++v) {
                    SampleDataLikelihood* m = NULL;
                    for (vector<SampleDataLikelihood>::iterator d = v->begin(); d != v->end(); ++d) {
                        if (d->genotype->matchesAlleles(alleles)) {
                            m = &*d;
                            break;
                        }
                    }
                    assert(m != NULL);
                    glMaxBasedOnAlts.push_back(m);
                }
                glMaxComboBasedOnAltsByPop[population].push_back(glMaxBasedOnAlts);
            }
            list<GenotypeCombo> glMaxBasedOnAltsGenotypeCombos; // build new combos into this list
            combinePopulationCombos(glMaxBasedOnAltsGenotypeCombos, glMaxComboBasedOnAltsByPop);
            bestCombo = glMaxBasedOnAltsGenotypeCombos.front();
        }

        DEBUG("best combo: " << bestCombo);

        // output

        if ((!alts.empty() && (1 - pHom.ToDouble()) >= parameters.PVL) || parameters.PVL == 0){

            // write the last gVCF record(s)
            if (parameters.gVCFout && !nonCalls.empty()) {
                vcflib::Variant var(parser->variantCallFile);
                out << results.gvcf(var, nonCalls, parser) << endl;
                nonCalls.clear();
            }

            vcflib::Variant var(parser->variantCallFile);

            out << results.vcf(
                var,
                pHom,
                bestComboOddsRatio,
                samples,
                referenceBase,
                alts,
                repeats,
                genotypingTotalIterations,
                parser->sampleList,
                coverage,
                bestCombo,
                alleleGroups,
                partialObservationGroups,
                partialObservationSupport,
                genotypesByPloidy,
                parser->sequencingTechnologies,
                parser)
                << endl;

        } else if (parameters.gVCFout) {
            // record statistics for gVCF output
            nonCalls.record(parser->currentSequenceName, parser->currentPosition, samples);
        }
        DEBUG2("finished position");

    }

    // write the last gVCF record
    if (parameters.gVCFout && !nonCalls.empty() && !parameters.gVCFNoChunk) {
        Results results;
        vcflib::Variant var(parser->variantCallFile);
        out << results.gvcf(var, nonCalls, parser) << endl;
        nonCalls.clear();
    }

    DEBUG("total sites: " << total_sites << endl
          << "processed sites: " << processed_sites << endl
          << "ratio: " << (float) processed_sites / (float) total_sites);

    delete parser;

    return 0;
}

} // namespace panfreebayes
