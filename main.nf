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

process DADA2_DENOISE {

    input:
    path demux_full

    output:
    path "table.qza", emit: table
    path "rep-seqs.qza", emit: rep_seqs
    path "denoising-stats.qza", emit: stats
    path "base-transition-stats.qza", emit: base_transition

    script:
    """
    qiime dada2 denoise-paired \
        --i-demultiplexed-seqs ${demux_full} \
        --p-trim-left-f ${params.trim_left_f} \
        --p-trim-left-r ${params.trim_left_r} \
        --p-trunc-len-f ${params.trunc_len_f} \
        --p-trunc-len-r ${params.trunc_len_r} \
        --p-n-threads ${task.cpus} \
        --o-table table.qza \
        --o-representative-sequences rep-seqs.qza \
        --o-denoising-stats denoising-stats.qza \
        --o-base-transition-stats base-transition-stats.qza
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
    DADA2_DENOISE(DEMUX.out.seqs)
    DADA2_DENOISE.out.table.view()
    DADA2_DENOISE.out.rep_seqs.view()
    DADA2_DENOISE.out.stats.view()
    DADA2_DENOISE.out.base_transition.view()
    
}
