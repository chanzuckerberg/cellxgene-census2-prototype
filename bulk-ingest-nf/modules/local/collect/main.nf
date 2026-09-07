process COLLECT {
    label 'small'

    input:
    path receipts

    output:
    path 'run_summary.json', emit: summary

    script:
    """
    cxg-collect ${receipts} --out run_summary.json
    """
}
