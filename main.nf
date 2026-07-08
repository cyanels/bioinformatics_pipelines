nextflow.enable.dsl = 2

process IMPORT_RAW {

    input:
    path raw_seq_dir

    output:
    path "emp-paired-end-sequences.qza"

    script:
    """
    qiime tools import \\
        --type EMPPairedEndSequences \
        --input-path ${raw_seq_dir} \
        --output-path emp-paired-end-sequences.qza
    """
}

process DEMUX {

    input:
    path emp_paired_end_sequences
    path sample_metadata

    output:
    path "demux-full.qza", emit: seqs
    path "demux-details.qza", emit: details

    script:
    """
    qiime demux emp-paired \
        --m-barcodes-file ${sample_metadata} \
        --m-barcodes-column barcode-sequence \
        --p-rev-comp-mapping-barcodes \
        --i-seqs ${emp_paired_end_sequences} \
        --o-per-sample-sequences demux-full.qza \
        --o-error-correction-details demux-details.qza
    """
}

workflow {
    raw_seq_dir = Channel.fromPath(params.raw_seq_dir, checkIfExists: true)
    IMPORT_RAW(raw_seq_dir)
    IMPORT_RAW.out.view()
    sample_metadata = Channel.fromPath(params.sample_metadata, checkIfExists: true)
    DEMUX(IMPORT_RAW.out, sample_metadata)
    DEMUX.out.seqs.view()
    DEMUX.out.details.view()
}
