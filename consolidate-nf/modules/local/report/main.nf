process REPORT {
    label 'small'

    // The report describes a run, so it is published to the run's own place.
    publishDir path: params.report_dir ?: '.', mode: 'copy',
               enabled: params.report_dir as boolean

    input:
    path receipts

    output:
    path 'consolidation_report.json', emit: report

    script:
    """
    consolidation-report.py ${receipts} --out consolidation_report.json
    """
}
