process INGEST {
    label 'ingest'
    tag "${meta.dataset_version_id}"

    // Nextflow publishes a local destination from the task's declared outputs. A
    // remote one is written by the task straight to object storage, so publishing
    // is off and the store output is simply absent.
    // AWS HealthOmics validates any parameter whose value looks like an s3://
    // URI as an existing input object, and an output prefix is not one. So
    // `s3_//bucket/prefix` is accepted as an escape and normalized below. Both
    // spellings mean the destination is remote.
    //
    // The store moves, because copying it would duplicate gigabytes. The receipt
    // is copied, because REPORT consumes it from the work directory and 'move'
    // would take it away.
    publishDir path: "${params.outdir.toString().replaceFirst('^s3_//', 's3://')}",
               mode: 'move', pattern: '*.zarr',
               enabled: !params.outdir.toString().matches('^s3(://|_//).*')
    publishDir path: "${params.outdir.toString().replaceFirst('^s3_//', 's3://')}",
               mode: 'copy', pattern: '*.json',
               enabled: !params.outdir.toString().matches('^s3(://|_//).*')

    // Size from the manifest. The corpus is heavy-tailed: the median dataset is
    // ~0.23 GiB and the largest is ~50 GiB. A fixed request is therefore too
    // large for most tasks and too small for the largest. Each retry raises the
    // request.
    //
    // The 2.5x multiple is headroom for a harmonize() step that loads more of
    // the dataset eagerly. The current steps do not. The gene projection
    // remaps each dask row block when the store is written. Before the
    // projection, peak use was under ~10 GiB even for the largest h5ad. Measure
    // it again with the projection, and revisit the request against the task
    // metrics that Nextflow and HealthOmics record.
    memory { def gib = (meta.filesize as long) / (1024L ** 3)
             def want = (Math.max(4.0d, gib * 2.5d) * task.attempt as int).GB
             def cap  = params.max_memory as nextflow.util.MemoryUnit
             want > cap ? cap : want }

    // Room for the staged h5ad, and for the zarr when it is staged locally. The
    // request rises with the attempt, so a retry after a full disk can succeed.
    disk   { def gib = (meta.filesize as long) / (1024L ** 3)
             def want = Math.max(8.0d, gib * 2.5d) * task.attempt
             (want as int).GB }

    input:
    val meta
    // The same for every dataset. The config is optional. Reference data holds
    // the gene references and is required. It has to be staged as a workflow
    // input, because a HealthOmics task cannot fetch it.
    path config
    path refdata

    output:
    path "${meta.dataset_version_id}.json", emit: receipt
    path "${meta.dataset_version_id}.zarr", type: 'dir', optional: true, emit: store

    script:
    // Nextflow escapes a `path` value for the shell, which is why `config` and
    // `refdata` below are bare. It renders every other value verbatim, so a
    // string has to be quoted here and any quote inside it escaped the POSIX
    // way: close the quote, emit an escaped quote, reopen. A dataset title
    // carries an apostrophe often enough to matter, and so can an `--outdir`.
    // (A helper closure here does not compile: the script block cannot call
    // one, so the escape is written out at each use.)

    // The record travels as JSON so the script has the source URI and the
    // expected file size without re-reading the manifest.
    def record = groovy.json.JsonOutput.toJson(meta).replace("'", "'\\''")
    // Same rule as publishDir above, applied to where the task writes.
    def outdir = params.outdir.toString().replaceFirst('^s3_//', 's3://')
    def outarg = outdir.replace("'", "'\\''")
    def stage  = outdir.startsWith('s3://') ? '' : '--stage-local'
    // A HealthOmics parameter always arrives as a string, and a non-empty string
    // is true in Groovy, so a plain truth test would read "false" as true.
    def all    = params.all_slots.toString().toBoolean() ? '--all-slots' : ''
    def chk    = params.check.toString().toBoolean()     ? '--check'     : ''
    def cfg    = config  ? "--config ${config}"   : ''
    // `task.cpus` is an integer, so it renders verbatim and needs no quoting,
    // unlike the strings above. dask otherwise sizes its thread pool from the
    // host's CPUs, because a CPU share sets no cgroup quota to read.
    """
    cxg-ingest.py \\
        --record '${record}' \\
        --outdir '${outarg}' ${stage} \\
        --threads ${task.cpus} \\
        --refdata ${refdata} \\
        ${cfg} ${all} ${chk}
    """
}
