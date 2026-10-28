process REPORT {
    label 'small'

    // The report describes a run, so it is published to the run's own place.
    publishDir path: params.report_dir ?: '.', mode: 'copy',
               enabled: params.report_dir as boolean

    input:
    path receipts

    output:
    path 'compaction_report.json', emit: report

    script:
    """
    compaction-report.py ${receipts} --out compaction_report.json
    """
}
