/*
 * HHSUITE — local HH-suite annotation via hhblits (MSA) + hhsearch (Pfam/PDB70).
 */

process HHSUITE_LOCAL {
    tag "$sample_id"
    label 'process_medium'

    input:
    tuple val(sample_id), path(substrates), path(proteins)

    output:
    tuple val(sample_id), val('hhsuite'), path("${sample_id}_hhsuite.csv"), emit: annotations

    script:
    """
    run_hhsuite.py \\
        --substrates ${substrates} \\
        --proteins ${proteins} \\
        --sample ${sample_id} \\
        --pfam-db "${params.hhsuite_pfam_db}" \\
        --pdb70-db "${params.hhsuite_pdb70_db}" \\
        --uniclust-db "${params.hhsuite_uniclust_db}" \\
        --output ${sample_id}_hhsuite.csv
    """
}
