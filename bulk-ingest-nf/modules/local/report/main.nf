process REPORT {
    label 'small'

    // The report describes a run, so it is published to the run's own place
    // rather than beside the stores. `enabled` has to be a static expression.
    publishDir path: params.report_dir ?: '.', mode: 'copy',
               enabled: params.report_dir as boolean

    input:
    path receipts

    output:
    path 'ingestion_report.json', emit: report

    script:
    """
    ingestion-report.py ${receipts} --out ingestion_report.json
    """
}
