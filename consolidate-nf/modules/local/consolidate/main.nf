process CONSOLIDATE {
    label 'consolidate'
    tag "${group.group}"

    // As INGEST in bulk-ingest-nf: Nextflow publishes a local destination from the
    // declared outputs, and the task writes a remote one itself. `s3_//` is the
    // escape for HealthOmics, which validates an s3:// parameter as an input.
    publishDir path: "${params.outdir.toString().replaceFirst('^s3_//', 's3://')}",
               mode: 'move', pattern: '*.zarr',
               enabled: !params.outdir.toString().matches('^s3(://|_//).*')
    publishDir path: "${params.outdir.toString().replaceFirst('^s3_//', 's3://')}",
               mode: 'copy', pattern: '*.json',
               enabled: !params.outdir.toString().matches('^s3(://|_//).*')

    // `consolidate` holds the `obs` of the group in memory and copies `X` one shard
    // at a time, so memory follows the cell count. The planner bounds the group.
    // Each retry raises the request.
    memory { def millions = (group.n_obs as long) / 1.0e6d
             def want = (Math.max(4.0d, 2.0d + 2.0d * millions) * task.attempt as int).GB
             def cap  = params.max_memory as nextflow.util.MemoryUnit
             want > cap ? cap : want }

    input:
    val group

    output:
    path "${group.group}.json", emit: receipt
    path "${group.group}.zarr", type: 'dir', optional: true, emit: store

    script:
    // Strings are quoted and escaped the POSIX way, as in bulk-ingest-nf.
    def record = groovy.json.JsonOutput.toJson(group).replace("'", "'\\''")
    def outdir = params.outdir.toString().replaceFirst('^s3_//', 's3://')
    def outarg = outdir.replace("'", "'\\''")
    def stage  = outdir.startsWith('s3://') ? '' : '--stage-local'
    def chk    = params.check.toString().toBoolean() ? '--check' : ''
    """
    consolidate-group.py \\
        --group '${record}' \\
        --outdir '${outarg}' ${stage} \\
        --threads ${task.cpus} ${chk}
    """
}
