process REPORT {
    label 'small'

    input:
    path receipts

    output:
    path 'ingestion_report.json', emit: report

    script:
    """
    ingestion-report.py ${receipts} --out ingestion_report.json
    """
}
