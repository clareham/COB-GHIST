/*
 * Shared helpers for the GHIST 2026 survey pipeline. Auto-compiled by
 * Nextflow from lib/ and usable from both *.config files and *.nf files,
 * which avoids relying on config-include ordering for a plain closure.
 */
class GhistUtils {

    /**
     * Resolve a GHIST data file and enforce the testing-only rule: this
     * pipeline is for pre-scoring method development, so it must never
     * read a *.final.* file (those are the held-out scoring sets).
     *
     * The data drop isn't consistently delimited - most files are
     * `<scenario>.testing.<ext>` / `<scenario>.final.<ext>`, but
     * relatedness2 is `relatedness2_testing.vcf` / `relatedness2_final.vcf`
     * (underscore). Both separators are accepted; either way, a "final"
     * token anywhere in the name is rejected outright.
     *
     * @param dataDir base directory holding the GHIST2026 data dump
     * @param name    file name, must contain a "testing" token (dot- or
     *                underscore-delimited) and no "final" token
     * @return the resolved java.io.File (may not exist - callers/channels
     *         should handle missing files, e.g. via `.filter { it.exists() }`)
     */
    static java.io.File testingFile(String dataDir, String name) {
        if (name =~ /[._]final[._]/) {
            throw new IllegalArgumentException(
                "GhistUtils.testingFile: '${name}' looks like a *.final.* (held-out scoring) file. " +
                "Only *.testing.* / *_testing.* sets may be used as pipeline input."
            )
        }
        if (!(name =~ /[._]testing[._]/)) {
            throw new IllegalArgumentException(
                "GhistUtils.testingFile: '${name}' is not a *.testing.* / *_testing.* file. " +
                "Only testing sets may be used as pipeline input - *.final.* sets are held out for scoring."
            )
        }
        return new java.io.File("${dataDir}/${name}")
    }
}
