/*
 * INTERPROSCAN — local domain annotation via interproscan.sh.
 */

process INTERPROSCAN_LOCAL {
    tag "$sample_id"
    label 'process_high'

    input:
    tuple val(sample_id), path(substrates), path(proteins)

    output:
    tuple val(sample_id), val('interproscan'), path("${sample_id}_interproscan.csv"), emit: annotations

    script:
    """
    run_interproscan.py \\
        --substrates ${substrates} \\
        --proteins ${proteins} \\
        --sample ${sample_id} \\
        --db "${params.interproscan_db}" \\
        --output ${sample_id}_interproscan.csv
    """
}
