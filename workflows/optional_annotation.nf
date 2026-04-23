/*
 * Phase 5: Optional Annotation
 *
 * All tools are independently skippable.
 * All database tools run locally as of v1.0.0 (offline-first).
 */

include { BLASTP_LOCAL }     from '../modules/local/blastp'
include { HHSUITE_LOCAL }    from '../modules/local/hhsuite'
include { INTERPROSCAN_LOCAL }  from '../modules/local/interproscan'
include { PROTPARAM }        from '../modules/local/protparam'

workflow OPTIONAL_ANNOTATION {
    take:
    ch_substrates  // [sample_id, substrates_filtered.tsv]
    ch_proteins    // [sample_id, proteins.faa]

    main:
    // Collect all annotation outputs into a single channel per sample
    ch_annotations = Channel.empty()

    // --- BLASTp (local only as of v1.0.0) ---
    if (!params.skip_blastp) {
        BLASTP_LOCAL(ch_substrates.join(ch_proteins))
        ch_annotations = ch_annotations.mix(BLASTP_LOCAL.out.annotations)
    }

    // --- HH-suite (local only as of v1.0.0) ---
    if (!params.skip_hhsuite) {
        HHSUITE_LOCAL(ch_substrates.join(ch_proteins))
        ch_annotations = ch_annotations.mix(HHSUITE_LOCAL.out.annotations)
    }

    // --- InterProScan (local only as of v1.0.0) ---
    if (!params.skip_interproscan) {
        INTERPROSCAN_LOCAL(ch_substrates.join(ch_proteins))
        ch_annotations = ch_annotations.mix(INTERPROSCAN_LOCAL.out.annotations)
    }

    // --- ProtParam (always local, lightweight) ---
    if (!params.skip_protparam) {
        PROTPARAM(ch_substrates.join(ch_proteins))
        ch_annotations = ch_annotations.mix(PROTPARAM.out.annotations)
    }

    emit:
    annotations = ch_annotations  // [sample_id, tool_name, annotations.csv]
}
