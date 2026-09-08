#!/usr/bin/env nextflow

include { INGEST } from './modules/local/ingest'
include { REPORT } from './modules/local/report'

// The pipeline consumes a manifest that `tools/cxg-manifest.py` prepared. It never
// calls the Discover API, so no task needs public internet access.

workflow {

    if (!params.manifest) {
        error "Missing --manifest. Prepare one first:\n" +
              "  tools/cxg-manifest.py --limit 5 --sort size_asc --out manifest.jsonl"
    }

    // Executors whose workers do not share this filesystem cannot publish to a
    // local path. The profile declares that. The rule then travels with the
    // profile, and a test can check it without the platform.
    // `s3_//` is an accepted escape for `s3://`; see modules/local/ingest/main.nf.
    def outdir = params.outdir.toString().replaceFirst('^s3_//', 's3://')
    // Coerced for the same reason as the flags in modules/local/ingest: a
    // parameter that arrives as a string would otherwise always be true.
    if (params.remote_outdir_required.toString().toBoolean() && !outdir.startsWith('s3://')) {
        error "--outdir must be an s3:// URI for this execution profile.\n" +
              "  A local path cannot be published from a remote worker.\n" +
              "  Got: ${params.outdir}"
    }

    // Profiles that cannot pull a public image declare that, for the same reason
    // as remote_outdir_required: the rule travels with the profile.
    if (params.container_required.toString().toBoolean() && !params.container) {
        error "--container must be a private ECR image URI for this execution profile."
    }

    // JSONL: one dataset per line. Each record travels whole, because the
    // memory directive reads `filesize` from it.
    datasets = Channel
        .fromPath(params.manifest, checkIfExists: true)
        .splitText()
        .map { line -> line.trim() }
        .filter { it }
        .map { line -> new groovy.json.JsonSlurper().parseText(line) }
        // Without this an empty manifest gives a successful run that writes
        // nothing: no task runs, so REPORT never runs either.
        .ifEmpty { error "No dataset records in ${params.manifest}" }

    // Optional and shared by every task. An empty list stages nothing.
    ch_config  = Channel.value(params.config  ? file(params.config,  checkIfExists: true) : [])
    ch_refdata = Channel.value(params.refdata ? file(params.refdata, checkIfExists: true) : [])

    INGEST(datasets, ch_config, ch_refdata)

    REPORT(INGEST.out.receipt.collect())

    REPORT.out.report.view { "ingestion report: ${it}" }
}
