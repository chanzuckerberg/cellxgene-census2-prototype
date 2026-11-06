#!/usr/bin/env nextflow

include { CONSOLIDATE } from './modules/local/consolidate'
include { REPORT  } from './modules/local/report'

// The pipeline consumes a plan that `tools/consolidation-plan.py` prepared from a
// catalog. Each line of the plan is one group of per-dataset stores, and each group
// becomes one consolidated store. The pipeline reads only object storage.

workflow {

    if (!params.plan) {
        error "Missing --plan. Prepare one first:\n" +
              "  tools/consolidation-plan.py catalog.parquet --out plan.jsonl"
    }

    // `s3_//` is an accepted escape for `s3://`; see modules/local/consolidate/main.nf.
    def outdir = params.outdir.toString().replaceFirst('^s3_//', 's3://')
    // A parameter that arrives as a string would otherwise always be true.
    if (params.remote_outdir_required.toString().toBoolean() && !outdir.startsWith('s3://')) {
        error "--outdir must be an s3:// URI for this execution profile.\n" +
              "  A local path cannot be published from a remote worker.\n" +
              "  Got: ${params.outdir}"
    }
    if (params.container_required.toString().toBoolean() && !params.container) {
        error "--container must be a private ECR image URI for this execution profile."
    }

    // JSONL: one group per line. Each record travels whole, because the memory
    // directive reads `n_obs` from it.
    groups = Channel
        .fromPath(params.plan, checkIfExists: true)
        .splitText()
        .map { line -> line.trim() }
        .filter { it }
        .map { line -> new groovy.json.JsonSlurper().parseText(line) }
        .ifEmpty { error "No groups in ${params.plan}" }

    CONSOLIDATE(groups)

    // An ignored CONSOLIDATE task writes no receipt. Its datasets stay in their
    // per-dataset stores, which is a correct catalog. If every task is ignored,
    // fail the run instead.
    REPORT(CONSOLIDATE.out.receipt.collect().ifEmpty { error "No CONSOLIDATE task wrote a receipt" })

    REPORT.out.report.view { "consolidation report: ${it}" }
}
