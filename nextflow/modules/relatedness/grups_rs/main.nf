/*
 * GRUPS-rs
 * Status: stub - Rust crate, would need `cargo install` (a source build,
 * just in a different toolchain than C/GSL). Skipped per the same
 * compile-from-source policy as RAiSD. To activate: `cargo install
 * grups-rs` manually, flip `enabled: true`, and fill in the real
 * invocation (pedigree-simulation-based background estimation).
 */
process GRUPS_RS {
    tag "${meta.challenge}"

    input:
    tuple val(meta), path(input_file)

    output:
    tuple val(meta), path("${meta.challenge}.grups_rs.stub"), emit: result

    script:
    """
    echo "GRUPS-rs is not wired up yet: ${meta.reason}" >&2
    exit 1
    """
}
