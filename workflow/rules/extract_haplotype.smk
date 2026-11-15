#
# Haplotype extraction -- a third way to reach a reference, alongside
# extract_region.smk (padded bubbles) and extract_path.smk (one named path).
#
# For one configured "<sample>#<hap>" prefix (config["haplotypes"], e.g.
# "HG00423#1"), extracts EVERY graph path whose name starts with
# "<sample>#<hap>#" -- i.e. all of that haplotype's contigs -- into ONE
# multi-record FASTA. That FASTA is a single reference: align_region makes
# one BAM against it and call_variants one VCF (panfreebayes streams the
# whole BAM across every contig), not one job per contig.
#
# Graph format is picked from the extension of config["graph"]:
#   .gbz  vg paths -x <graph> -F -Q '<sample>#<hap>#'   (one command; -Q
#         selects paths by name prefix). Not `-L` then `-F -p <names>`:
#         on vg v1.71.0, `-L`/`-L -H` fail with "graph does not contain any
#         paths" on a GBZ with no reference sample set, and -H can't be
#         combined with -p ("multiple selection criteria").
#         REQUIREMENT: the GBZ must have at least one reference sample
#         set, or vg paths refuses the graph entirely. Make one with:
#           vg gbwt -Z in.gbz --set-reference <SAMPLE> [--set-reference ...] \
#             --gbz-format -g out.gbz
#         Haplotype path names carry a trailing phase-block field, e.g.
#         "HG00423#1#JBHIIO010000010.1#0"; FASTA headers are kept exactly as
#         vg prints them.
#   .og   odgi paths -i <graph> -f, keeping records whose header starts
#         with the prefix (same mechanism extract_path.smk uses); the
#         match is a literal string compare (awk index() == 1), not a
#         regex, so '.' etc. in sample names can't over-match.
#
# Output goes under the same regions/<id>.fasta pattern as the other two
# extraction rules, so align_region/call_variants need no changes; the
# {hap_id} constraint (common.smk's HAP_ID_REGEX) keeps this rule's outputs
# disjoint from extract_path's and extract_region's.
#
# The contig names and lengths that went into the reference are written to
# the log (from the .fai), and the rule fails clearly if no path matches.
#
rule extract_haplotype:
    input:
        graph=config["graph"],
    output:
        fasta=f"{RESULTS}/regions/{{hap_id}}.fasta",
        fai=f"{RESULTS}/regions/{{hap_id}}.fasta.fai",
    wildcard_constraints:
        hap_id=HAP_ID_REGEX,
    params:
        prefix=lambda wc: hap_prefix_for(wc.hap_id) + "#",
    log:
        f"{RESULTS}/logs/extract_haplotype/{{hap_id}}.log",
    shell:
        r"""
        tmp={output.fasta}.tmp
        : > {log}
        case "{input.graph}" in
          *.gbz)
            vg paths -x {input.graph} -F -Q '{params.prefix}' > "$tmp" 2>> {log}
            ;;
          *.og)
            odgi paths -i {input.graph} -f 2>> {log} | awk -v p=">{params.prefix}" '
              /^>/ {{ keep = (index($0, p) == 1) }}
              keep
            ' > "$tmp"
            ;;
          *)
            echo "extract_haplotype: unsupported graph format (need .gbz or .og): {input.graph}" >&2
            exit 1
            ;;
        esac
        if ! grep -q '^>' "$tmp"; then
          rm -f "$tmp"
          echo "extract_haplotype: no graph path starts with '{params.prefix}' in {input.graph}" >&2
          exit 1
        fi
        mv "$tmp" {output.fasta}
        samtools faidx {output.fasta} 2>> {log}
        {{
          echo "haplotype {params.prefix}: $(wc -l < {output.fai} | tr -d ' ') contig(s)"
          printf 'contig\tlength\n'
          cut -f1,2 {output.fai}
        }} >> {log}
        """
