process INGEST {
    label 'ingest'
    tag "${meta.dataset_version_id}"

    // A local destination is published by Nextflow from the task's declared
    // outputs, which is why nothing needs to be bind-mounted. A remote one is
    // written by the task straight to object storage, so publishing is off and
    // the store output is simply absent.
    // AWS HealthOmics validates any parameter whose value looks like an s3://
    // URI as an existing input object, and an output prefix is not one. So
    // `s3_//bucket/prefix` is accepted as an escape and normalized below. Both
    // spellings mean the destination is remote.
    //
    // The store moves, because copying it would duplicate gigabytes. The receipt
    // is copied, because COLLECT consumes it from the work directory and 'move'
    // would take it away.
    publishDir path: "${params.outdir.toString().replaceFirst('^s3_//', 's3://')}/v1",
               mode: 'move', pattern: '*.zarr',
               enabled: !params.outdir.toString().matches('^s3(://|_//).*')
    publishDir path: "${params.outdir.toString().replaceFirst('^s3_//', 's3://')}/v1",
               mode: 'copy', pattern: '*.json',
               enabled: !params.outdir.toString().matches('^s3(://|_//).*')

    // Size from the manifest. The corpus is heavy-tailed: the median dataset is
    // ~0.23 GiB and the largest is ~50 GiB, so a fixed request either wastes
    // most of the fleet or dies on the tail. Each retry raises the request.
    memory { def gib = (meta.filesize as long) / (1024L ** 3)
             def want = Math.max(4.0d, gib * 2.5d) * task.attempt
             (Math.min(want, 128.0d) as int).GB }

    // Staged h5ad plus the local zarr, both on scratch before the sync.
    disk   { def gib = (meta.filesize as long) / (1024L ** 3)
             (Math.max(8.0d, gib * 2.5d) as int).GB }

    input:
    val meta
    // Optional, and the same for every dataset. Reference data has to be staged
    // as a workflow input, because a HealthOmics task cannot fetch it.
    path config
    path refdata

    output:
    path "${meta.dataset_version_id}.json", emit: receipt
    path "${meta.dataset_version_id}.zarr", type: 'dir', optional: true, emit: store

    script:
    // The record travels as JSON so the script has the source URI and the
    // expected file size without re-reading the manifest.
    def record = groovy.json.JsonOutput.toJson(meta)
    // With a local destination the task writes into its work directory and the
    // workflow publishes; with a remote one the task writes to object storage.
    def outdir = params.outdir.toString().replaceFirst('^s3_//', 's3://')
    def stage  = outdir.startsWith('s3://') ? '' : '--stage-local'
    def all    = params.all_slots ? '--all-slots' : ''
    def chk    = params.check     ? '--check'     : ''
    def cfg    = config  ? "--config ${config}"   : ''
    def ref    = refdata ? "--refdata ${refdata}" : ''
    """
    cxg-ingest \\
        --record '${record}' \\
        --outdir '${outdir}' ${stage} \\
        ${cfg} ${ref} ${all} ${chk}
    """
}
