nextflow.enable.dsl = 2

// Import raw multiplexed EMP-paired-end reads into a QIIME2 artifact.
// Entry point for demux_only and full workflows
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

// Demultiplex raw reads using barcode metadata.
// NOTE: --p-rev-comp-mapping-barcodes is currently hardcoded to true because
// that's what this dataset (Atacama soil, EMP-paired format) requires.
// TODO (on the horizon): expose this as params.rev_comp_barcodes so it's
// configurable per-dataset instead of assumed.
process DEMUX {
    publishDir "${params.outdir}/demux", mode: 'copy'

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

// Generate a visual summary (.qzv) of demultiplexed read quality/counts.
// This is the file to inspect (via QIIME2 View) to pick DADA2's
// trim/truncation parameters below.
process DEMUX_SUMMARIZE {
    publishDir "${params.outdir}/demux", mode: 'copy'

    input:
    path demux_seqs

    output:
    path "demux.qzv"

    script:
    """
    qiime demux summarize \
        --i-data ${demux_seqs} \
        --o-visualization demux.qzv
    """
}

// Denoise paired-end reads with DADA2: corrects sequencing errors and
// resolves reads into exact Amplicon Sequence Variants (ASVs), producing
// a feature table + representative sequences. This replaces traditional
// OTU clustering entirely — no separate clustering step is needed.
//
// trim_left_f/r and trunc_len_f/r are set in nextflow.config based on
// visual inspection of demux.qzv's quality plots (see DEMUX_SUMMARIZE above).
process DADA2_DENOISE {
    publishDir "${params.outdir}/dada2", mode: 'copy'

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

// Remove low-prevalence ASVs from the feature table before downstream
// taxonomy classification and diversity analysis.
// --p-min-samples 2 drops singletons (ASVs found in exactly 1 sample),
// which are frequently sequencing artifacts or extremely rare/noise taxa
// rather than genuine, reproducible signal.
process FILTER_TABLE {

    input:
    path asv_table

    output:
    path "asv_table_filtered.qza"

    script:
    """
    qiime feature-table filter-features \
        --i-table ${asv_table} \
        --p-min-samples 2 \
        --o-filtered-table asv_table_filtered.qza
    """
}

// Filter representative sequences to match the ASVs remaining in the
// filtered feature table (FILTER_TABLE above), keeping table and rep-seqs
// in sync.
process FILTER_SEQUENCES {

    input:
    path asv_table_filtered
    path rep_seqs

    output:
    path "asv_seqs_filtered.qza"

    script:
    """
    qiime feature-table filter-seqs \
        --i-data ${rep_seqs} \
        --i-table ${asv_table_filtered} \
        --o-filtered-data asv_seqs_filtered.qza
    """
}

// Generate a visual summary (.qzv) of filtered ASVs.
// Inspect sample and feature frequency distributions and total counts.
process FILTER_SUMMARIZE{
    publishDir "${params.outdir}/dada2", mode: 'copy', pattern: "*.qzv"

    input:
    path asv_table_filtered
    path sample_metadata

    output:
    path "feature-freqs.qza", emit: feature_freqs
    path "sample-freqs.qza", emit: sample_freqs
    path "asv_filter_summary.qzv", emit: asv_filter_summary

    script:
    """
    qiime feature-table summarize \
        --i-table ${asv_table_filtered} \
        --m-metadata-file ${sample_metadata} \
        --o-feature-frequencies feature-freqs.qza \
        --o-sample-frequencies sample-freqs.qza \
        --o-summary asv_filter_summary.qzv
    """
}
// Classify query sequences using classify-consensus-vsearch.
// Must have SILVA reference files downloaded first.
process CLASSIFY_TAXONOMY{
    publishDir "${params.outdir}/taxonomy", mode: 'copy', pattern: "taxonomy.qza"

    input:
    path asv_seqs_filtered
    path SILVA_ref_reads
    path SILVA_ref_taxon

    output:
    path "taxonomy.qza", emit: taxonomy
    path "search_results.qza", emit: search_results

    script:
    """
    qiime feature-classifier classify-consensus-vsearch \
        --i-query ${asv_seqs_filtered} \
        --i-reference-reads ${SILVA_ref_reads} \
        --i-reference-taxonomy ${SILVA_ref_taxon} \
        --p-threads ${task.cpus} \
        --p-perc-identity ${params.perc_identity} \
        --o-classification taxonomy.qza \
        --o-search-results search_results.qza
    """
}

// Generate a visual summary (.qzv) of taxonomy classifications.
// Inspect this to check assignment rates and look for mitochondria/
// chloroplast/unassigned entries before deciding on filter-table below.
process TAXONOMY_SUMMARIZE{
    publishDir "${params.outdir}/taxonomy", mode: 'copy'

    input:
    path taxonomy

    output:
    path "taxonomy_summary.qzv"

    script:
    """
    qiime metadata tabulate \
        --m-input-file ${taxonomy} \
        --o-visualization taxonomy_summary.qzv
    """
}
// Filter classified ASV table to remove mitochrondria/chloroplast/unassigned.
process FILTER_TAXONOMY{
      input:
    path asv_table_filtered
    path taxonomy

    output:
    path "taxon_table_filtered.qza"

    script:
    """
    qiime taxa filter-table \
        --i-table ${asv_table_filtered} \
        --i-taxonomy ${taxonomy} \
        --p-exclude mitochondria,chloroplast,unassigned \
        --o-filtered-table taxon_table_filtered.qza
    """
}

// Generate a visual summary (.qzv) of the taxonomy-filtered feature table.
// Compare feature/sample counts against FILTER_SUMMARIZE's output to see
// how many ASVs were removed by mitochondria/chloroplast/unassigned filtering.
process FILTER_TAXONOMY_SUMMARIZE{
    publishDir "${params.outdir}/taxonomy", mode: 'copy', pattern: "*.qzv"

    input:
    path taxon_table_filtered
    path sample_metadata

    output:
    path "taxon_filtered_feature-freqs.qza", emit: feature_freqs
    path "taxon_filtered_sample-freqs.qza", emit: sample_freqs
    path "taxon_filtered_asv_filter_summary.qzv", emit: asv_filter_summary

    script:
    """
    qiime feature-table summarize \
        --i-table ${taxon_table_filtered} \
        --m-metadata-file ${sample_metadata} \
        --o-feature-frequencies taxon_filtered_feature-freqs.qza \
        --o-sample-frequencies taxon_filtered_sample-freqs.qza \
        --o-summary taxon_filtered_asv_filter_summary.qzv
    """
}

// Filter representative sequences to match the ASVs remaining after
// taxonomy-based filtering (FILTER_TAXONOMY above)
process FILTER_TAXONOMY_SEQUENCES {

    input:
    path asv_seqs_filtered
    path taxon_table_filtered

    output:
    path "taxon_seqs_filtered.qza"

    script:
    """
    qiime feature-table filter-seqs \
        --i-data ${asv_seqs_filtered} \
        --i-table ${taxon_table_filtered} \
        --o-filtered-data taxon_seqs_filtered.qza
    """
}

// Entry point: dispatches to one of two workflows based on params.stage.
workflow {
    if (params.stage == 'demux_only') {
        demux_only()
    } else {
        full()
    }
}

// Run this workflow first, before full(), to obtain demultiplexed data
// and inspect demux.qzv's quality plots. 
// The results determine the trim/truncation parameters needed for 
// DADA2_DENOISE in the full workflow (set those in nextflow.config once identified here).
workflow demux_only{
    raw_seq_dir = Channel.fromPath(params.raw_seq_dir, checkIfExists: true)
    IMPORT_RAW(raw_seq_dir)
    sample_metadata = Channel.fromPath(params.sample_metadata, checkIfExists: true)
    DEMUX(IMPORT_RAW.out, sample_metadata)
    DEMUX_SUMMARIZE(DEMUX.out.seqs)
    }

// Full pipeline: import + demux (or reuse existing demuxed data) through
// denoising and low-prevalence filtering.
workflow full{
    // If demuxed sequences already exist from a prior demux_only run,
    // reuse them (params.demux_seqs) instead of re-running import/demux.
    // Otherwise, run import + demux from raw data.
    if (params.demux_seqs) {
        sample_metadata = Channel.fromPath(params.sample_metadata, checkIfExists: true)
        demux_seqs = Channel.fromPath(params.demux_seqs, checkIfExists: true)
    } else {
        raw_seq_dir = Channel.fromPath(params.raw_seq_dir, checkIfExists: true)
        IMPORT_RAW(raw_seq_dir)
        sample_metadata = Channel.fromPath(params.sample_metadata, checkIfExists: true)
        DEMUX(IMPORT_RAW.out, sample_metadata)
        demux_seqs = DEMUX.out.seqs
    }
    DEMUX_SUMMARIZE(demux_seqs)
    DADA2_DENOISE(demux_seqs)
    FILTER_TABLE(DADA2_DENOISE.out.table)
    FILTER_SEQUENCES(FILTER_TABLE.out, DADA2_DENOISE.out.rep_seqs)
    FILTER_SUMMARIZE(FILTER_TABLE.out, sample_metadata)
    SILVA_ref_reads = Channel.fromPath(params.SILVA_ref_reads, checkIfExists: true)
    SILVA_ref_taxon = Channel.fromPath(params.SILVA_ref_taxon, checkIfExists: true)
    CLASSIFY_TAXONOMY(FILTER_SEQUENCES.out, SILVA_ref_reads, SILVA_ref_taxon)
    TAXONOMY_SUMMARIZE(CLASSIFY_TAXONOMY.out.taxonomy)
    FILTER_TAXONOMY(FILTER_TABLE.out, CLASSIFY_TAXONOMY.out.taxonomy)
    FILTER_TAXONOMY_SUMMARIZE(FILTER_TAXONOMY.out, sample_metadata)
    FILTER_TAXONOMY_SEQUENCES(FILTER_SEQUENCES.out, FILTER_TAXONOMY.out)
    }
