/*
 * fastsimcoal2
 * Status: runnable-unverified (bioconda: fastsimcoal2, cyvcf2)
 *
 * SIMPLIFICATION: builds a single-population folded SFS directly from the
 * VCF (bin/vcf_to_sfs.py) since no popmap is provided for the two-
 * population challenges - see that script's docstring.
 *
 * CAVEAT: the bioconda `fastsimcoal2` recipe ships its binary under a
 * version-tagged name (e.g. `fsc2709`/`fsc28`) that changes across
 * releases, and we haven't independently verified the exact .tpl/.est
 * template syntax below against a real run - treat the SFS-generation
 * half as solid and the actual fastsimcoal2 invocation as a best-effort
 * starting point to check against `fsc2709 --help` / the fastsimcoal2
 * manual before trusting results.
 */
process FASTSIMCOAL2 {
    tag "${meta.challenge}"
    publishDir { "${params.outdir}/${meta.category}/${meta.challenge}/fastsimcoal2" }, mode: 'copy'
    conda { meta.conda }

    input:
    tuple val(meta), path(input_file)

    output:
    tuple val(meta), path("${meta.challenge}.fsc2_out"), emit: result

    script:
    """
    vcf_to_sfs.py ${input_file} ${meta.challenge}_MAFpop0.obs

    cat > ${meta.challenge}.tpl <<'TPL'
//Minimal 1-population model (placeholder - see module header)
1 populations
//Population effective sizes
NPOP\$
//Samples sizes
1
//Growth rates
0
//Number of migration matrices
0
//historical event
0 historical event
//Number of independent loci
1 0
//Per chromosome: data type, num loci, rec rate, mut rate
FREQ 1 0 1e-8 OUTEXP
TPL

    cat > ${meta.challenge}.est <<'EST'
// Priors and rules
[PARAMETERS]
1 NPOP\$ unif 10 100000 output
[RULES]
[COMPLEX PARAMETERS]
EST

    # Binary name is version-tagged in bioconda (fsc2709, fsc28, ...) - resolve dynamically.
    FSC_BIN=\$(compgen -c fsc | grep -E '^fsc[0-9]*\$' | head -1)
    if [ -z "\$FSC_BIN" ]; then
        echo "ERROR: no fastsimcoal2 binary (fsc*) found on PATH" >&2
        exit 1
    fi
    "\$FSC_BIN" -t ${meta.challenge}.tpl -e ${meta.challenge}.est ${meta.args} -o
    mkdir -p ${meta.challenge}.fsc2_out
    mv ${meta.challenge}/* ${meta.challenge}.fsc2_out/ 2>/dev/null || true
    """
}

