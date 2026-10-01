#
# Stage 1 -- bubble discovery.
#
# This MUST be a `checkpoint`, not a plain `rule`: Snakemake normally builds
# its whole execution DAG before running anything, which requires knowing
# every wildcard value upfront. But the number and identity of qualifying
# bubbles is only known after this awk filter actually scans the real
# deconstruct VCF -- it depends on graph/strain data, not on config.yaml.
# `checkpoint` is Snakemake's mechanism for "pause DAG construction on this
# branch, run this step, then re-evaluate using its real output" -- see
# common.smk's bubble_row()/all_bubble_ids(), which call
# checkpoints.discover_bubbles.get() to do exactly that.
#
# The awk program below is a verbatim transplant of
# scripts/panfreebayes/discover_bubbles.sh's own awk program: same metric
# (|REF_len - ALT_len|), same strict `>` threshold (default 500, confirmed
# manual-pipeline value -- not a guess), same multi-allelic handling (pick
# the ALT with the largest size difference from REF; skip symbolic ALTs like
# <INS> rather than guess their length -- itself flagged in that script as
# this implementation's own design choice, untested against real
# multi-allelic data), same descending sort by size-diff, same optional
# --limit/head -n cap.
#
# GOTCHA specific to Snakemake's `shell:` (not present in the bash version):
# `shell:` strings go through Python's .format()-style substitution for
# {input}/{output}/{params}/{wildcards}, so every literal awk `{`/`}` below
# is doubled to `{{`/`}}`. Bash's own single-quoted awk needed no such
# escaping.
#
checkpoint discover_bubbles:
    input:
        vcf=config["deconstruct_vcf"],
    output:
        tsv=f"{RESULTS}/bubbles/bubbles.tsv",
    params:
        min_size_diff=config["min_size_diff"],
        limit=lambda wc: config.get("bubble_limit") or "",
    log:
        f"{RESULTS}/logs/discover_bubbles.log",
    shell:
        r"""
        (case "{input.vcf}" in
           *.gz) gzip -dc -- {input.vcf} ;;
           *)    cat -- {input.vcf} ;;
         esac) | awk -v min_size_diff={params.min_size_diff} '
          BEGIN {{ FS = "\t"; OFS = "\t" }}
          /^#/  {{ next }}
          NF < 5 {{ next }}
          {{
            chrom = $1; pos = $2 + 0; ref = $4; alt = $5
            ref_len = length(ref)

            n = split(alt, alts, ",")
            best_diff = -1; best_alt = ""; best_alt_len = 0
            for (i = 1; i <= n; i++) {{
              a = alts[i]
              if (a == "" || substr(a, 1, 1) == "<") continue
              alt_len = length(a)
              diff = ref_len - alt_len
              if (diff < 0) diff = -diff
              if (diff > best_diff) {{ best_diff = diff; best_alt = a; best_alt_len = alt_len }}
            }}
            if (best_diff < 0) next
            if (best_diff <= min_size_diff) next

            start = pos
            end   = pos + ref_len - 1
            print chrom, start, end, ref_len, best_alt_len, best_diff, ref, best_alt
          }}
        ' | sort -t "$(printf '\t')" -k6,6nr > {output.tsv}.body 2> {log}

        {{
          printf 'path\tstart\tend\tref_len\talt_len\tsize_diff\tref\talt\n'
          if [ -n "{params.limit}" ]; then
            head -n {params.limit} {output.tsv}.body
          else
            cat {output.tsv}.body
          fi
        }} > {output.tsv}

        echo "$(wc -l < {output.tsv}.body | tr -d ' ') candidate site(s) with size-diff > {params.min_size_diff}" >> {log}
        rm -f {output.tsv}.body
        """
