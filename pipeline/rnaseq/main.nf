nextflow.enable.dsl = 2

process QC {
    tag "${sample_id}"
    cpus 4
    memory '8 GB'

    publishDir "${params.outdir}/qc/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path(reads), emit: reads_pass
    path '*_fastqc.html'
    path '*_fastqc.zip'

    script:
    """
    fastqc -t ${task.cpus} -o . ${reads[0]} ${reads[1]}
    """
}

process MAPPING {
    tag "${sample_id}"
    cpus 16
    memory '64 GB'

    input:
    tuple val(sample_id), path(reads)
    path 'star_index/*'

    output:
    tuple val(sample_id), path("${sample_id}.Aligned.sortedByCoord.out.bam"), emit: bam
    path "${sample_id}.Log.final.out", emit: log_final

    script:
    """
    STAR --runThreadN ${task.cpus} \\
         --genomeDir star_index \\
         --readFilesIn ${reads[0]} ${reads[1]} \\
         --readFilesCommand gunzip -c \\
         --outSAMtype BAM SortedByCoordinate \\
         --outFileNamePrefix ${sample_id}.
    """
}

process COUNT {
    cpus 8
    memory '32 GB'

    publishDir "${params.outdir}/counts", mode: 'copy'

    input:
    tuple path(bams), path(gtf)

    output:
    path 'gene_counts.txt'
    path 'gene_counts.txt.summary'

    script:
    def bam_args = (bams instanceof List ? bams : [bams]).collect { it.toString() }.join(' ')
    """
    featureCounts -T ${task.cpus} -a ${gtf} -o gene_counts.txt ${bam_args}
    """
}

workflow {
    read_pairs_ch = Channel.fromFilePairs("${params.fastq_dir}/*_{1,2}.fastq.gz")
    star_index_files = Channel.fromPath("${params.star_index}/*").collect()
    gtf_file = channel.fromPath(params.gtf).first()

    QC(read_pairs_ch)
    MAPPING(QC.out.reads_pass, star_index_files)

    count_in = MAPPING.out.bam
        .map { it -> it[1] }
        .collect()
        .combine(gtf_file)

    COUNT(count_in)
}
