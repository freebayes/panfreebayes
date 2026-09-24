//
// panfreebayes CLI (Milestone 3)
//
//   panfreebayes call --ref <ref.fasta> --bam <aln.bam> [flags] [-- <extra engine args>] > out.vcf
//   panfreebayes version
//   panfreebayes help
//
// `call` runs the extracted FreeBayes engine over the ENTIRE reference / BAM as a
// single region (no -r/-t). Output is bit-for-bit identical to stock
// single-process `freebayes` run with the equivalent arguments.
//

#include "panfreebayes_core.h"

#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <set>
#include <sstream>
#include <string>
#include <vector>
#include <csignal>

#include "SeqLib/BamReader.h"
#include "SeqLib/BamHeader.h"
#include "SegfaultHandler.h"

namespace {

const char* USAGE =
    "usage: panfreebayes <command> [options]\n"
    "\n"
    "commands:\n"
    "  call      call variants on a local reference + an aligned BAM\n"
    "  version   print version\n"
    "  help      print this message\n"
    "\n"
    "panfreebayes call --ref <ref.fasta> --bam <aln.bam> [options] > out.vcf\n"
    "\n"
    "  --ref FILE                   local reference FASTA (required)\n"
    "  --bam FILE                   coordinate-sorted BAM of reads aligned to --ref (required)\n"
    "  --pooled-continuous          treat the sample as a pool; report continuous allele frequencies\n"
    "  --min-alternate-count N      min observations supporting an alternate allele\n"
    "  --min-alternate-fraction X   min fraction of observations supporting an alternate allele\n"
    "  --limit-coverage N           downsample to N-fold coverage (seeded srand(13), deterministic)\n"
    "  -- <args...>                  pass <args> verbatim to the engine (advanced;\n"
    "                               region / target / stdin flags are rejected)\n"
    "\n"
    "The whole reference is analysed as one region. To call on a sub-region,\n"
    "extract it into its own reference FASTA + BAM first.\n";

bool fileReadable(const std::string& p) {
    std::ifstream f(p.c_str());
    return f.good();
}

int toInt(const std::string& s, const char* flag) {
    try {
        size_t pos = 0;
        int v = std::stoi(s, &pos);
        if (pos == s.size()) return v;
    } catch (...) {}
    std::cerr << "panfreebayes call: " << flag << " expects an integer, got '" << s << "'\n";
    std::exit(2);
}

double toDouble(const std::string& s, const char* flag) {
    try {
        size_t pos = 0;
        double v = std::stod(s, &pos);
        if (pos == s.size()) return v;
    } catch (...) {}
    std::cerr << "panfreebayes call: " << flag << " expects a number, got '" << s << "'\n";
    std::exit(2);
}

// reference sequence names: prefer the .fai, else scan '>' lines
std::set<std::string> fastaSeqNames(const std::string& ref) {
    std::set<std::string> names;
    std::ifstream fai((ref + ".fai").c_str());
    if (fai) {
        std::string line;
        while (std::getline(fai, line)) {
            if (!line.empty()) names.insert(line.substr(0, line.find('\t')));
        }
        return names;
    }
    std::ifstream fa(ref.c_str());
    std::string line;
    while (std::getline(fa, line)) {
        if (!line.empty() && line[0] == '>') {
            std::string n = line.substr(1);
            size_t ws = n.find_first_of(" \t\r");
            if (ws != std::string::npos) n.resize(ws);
            names.insert(n);
        }
    }
    return names;
}

// 0 = ok; non-zero (with a message on stderr) = do not proceed
int preflight(const std::string& ref, const std::string& bam) {
    if (ref.empty() || bam.empty()) {
        std::cerr << "panfreebayes call: --ref and --bam are both required\n\n" << USAGE;
        return 2;
    }
    if (!fileReadable(ref)) {
        std::cerr << "panfreebayes call: cannot read reference FASTA: " << ref << "\n";
        return 2;
    }
    if (!fileReadable(bam)) {
        std::cerr << "panfreebayes call: cannot read BAM: " << bam << "\n";
        return 2;
    }

    SeqLib::BamReader reader;
    if (!reader.Open(bam)) {
        std::cerr << "panfreebayes call: could not open BAM (corrupt, truncated, or not a BAM?): " << bam << "\n";
        return 2;
    }
    SeqLib::BamHeader hdr = reader.Header();
    if (hdr.isEmpty()) {
        std::cerr << "panfreebayes call: BAM has no readable header: " << bam << "\n";
        return 2;
    }

    if (hdr.AsString().find("SO:coordinate") == std::string::npos) {
        std::cerr << "panfreebayes call: BAM is not coordinate-sorted (need @HD ... SO:coordinate).\n"
                  << "  sort it first:  samtools sort -o sorted.bam " << bam << "\n";
        return 2;
    }

    std::set<std::string> refNames = fastaSeqNames(ref);
    if (refNames.empty()) {
        std::cerr << "panfreebayes call: no sequences found in reference FASTA: " << ref << "\n";
        return 2;
    }

    SeqLib::HeaderSequenceVector sq = hdr.GetHeaderSequenceVector();
    if (sq.empty()) {
        std::cerr << "panfreebayes call: BAM header has no @SQ (reference) lines\n";
        return 2;
    }
    size_t matched = 0;
    for (size_t i = 0; i < sq.size(); ++i) if (refNames.count(sq[i].Name)) ++matched;

    if (matched == 0) {
        std::cerr << "panfreebayes call: BAM and reference do not match -- none of the BAM's\n"
                  << "  reference sequences are present in " << ref << "\n"
                  << "  BAM @SQ: ";
        for (size_t i = 0; i < sq.size() && i < 4; ++i) std::cerr << (i ? ", " : "") << sq[i].Name;
        std::cerr << (sq.size() > 4 ? ", ...\n" : "\n");
        return 2;
    }
    if (matched < sq.size()) {
        std::cerr << "panfreebayes call: warning: " << (sq.size() - matched) << " of " << sq.size()
                  << " BAM reference sequences are absent from the FASTA "
                  << "(analysing the " << matched << " that match)\n";
    }
    return 0;
}

int cmd_call(int argc, char** argv) {
    panfreebayes::Options opt;
    bool sawExtra = false;

    for (int i = 0; i < argc; ++i) {
        std::string a = argv[i];
        auto need = [&](const char* what) -> std::string {
            if (i + 1 >= argc) {
                std::cerr << "panfreebayes call: " << what << " needs an argument\n";
                std::exit(2);
            }
            return std::string(argv[++i]);
        };
        if      (sawExtra)                         opt.extraArgs.push_back(a);
        else if (a == "--")                        sawExtra = true;
        else if (a == "--ref")                     opt.fasta = need("--ref");
        else if (a == "--bam")                     opt.bam = need("--bam");
        else if (a == "--pooled-continuous")       opt.pooledContinuous = true;
        else if (a == "--min-alternate-count")     opt.minAlternateCount = toInt(need("--min-alternate-count"), "--min-alternate-count");
        else if (a == "--min-alternate-fraction")  opt.minAlternateFraction = toDouble(need("--min-alternate-fraction"), "--min-alternate-fraction");
        else if (a == "--limit-coverage")          opt.limitCoverage = toInt(need("--limit-coverage"), "--limit-coverage");
        else if (a == "-h" || a == "--help")     { std::cout << USAGE; return 0; }
        else {
            std::cerr << "panfreebayes call: unrecognised argument '" << a << "'\n\n" << USAGE;
            return 2;
        }
    }

    int pf = preflight(opt.fasta, opt.bam);
    if (pf != 0) return pf;

    return panfreebayes::callVariants(opt, std::cout);
}

} // namespace

int main(int argc, char** argv) {
    signal(SIGSEGV, segfaultHandler);

    if (argc < 2) { std::cerr << USAGE; return 1; }
    std::string cmd = argv[1];

    if (cmd == "call") return cmd_call(argc - 2, argv + 2);
    if (cmd == "version" || cmd == "--version" || cmd == "-v") {
        std::cout << panfreebayes::version() << "\n";
        return 0;
    }
    if (cmd == "help" || cmd == "--help" || cmd == "-h") {
        std::cout << USAGE;
        return 0;
    }

    std::cerr << "panfreebayes: unknown command '" << cmd << "'\n\n" << USAGE;
    return 1;
}
