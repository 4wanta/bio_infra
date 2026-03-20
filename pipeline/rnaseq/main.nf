nextflow.enable.dsl = 2

process TRIMMING_PE {
    tag "${sample_id}"
    cpus 4
    memory '8 GB'

    publishDir "${params.outdir}/trimmed/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_*.fastq.P.qtrim.gz"), emit: reads_pass
    path "${sample_id}_1.fastq.U.qtrim.gz", optional: true
    path "${sample_id}_2.fastq.U.qtrim.gz", optional: true

    script:
    """
    TRIMMOMATIC_JAR=\$(ls -1 /opt/conda/share/trimmomatic-*/trimmomatic*.jar 2>/dev/null | awk 'NR==1{print; exit}')
    ADAPTER_DIR=\$(ls -d /opt/conda/share/trimmomatic-*/adapters 2>/dev/null | awk 'NR==1{print; exit}')

    if [[ -z "\${TRIMMOMATIC_JAR}" || -z "\${ADAPTER_DIR}" ]]; then
        echo "Trimmomatic jar or adapters not found under /opt/conda/share" >&2
        exit 1
    fi

    java -jar "\${TRIMMOMATIC_JAR}" PE -threads ${task.cpus} -phred33 \\
        ${reads[0]} ${reads[1]} \\
        ${sample_id}_1.fastq.P.qtrim.gz ${sample_id}_1.fastq.U.qtrim.gz \\
        ${sample_id}_2.fastq.P.qtrim.gz ${sample_id}_2.fastq.U.qtrim.gz \\
        ILLUMINACLIP:\${ADAPTER_DIR}/TruSeq3-PE.fa:2:30:10 \\
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    """
}

process TRIMMING_SE {
    tag "${sample_id}"
    cpus 4
    memory '8 GB'

    publishDir "${params.outdir}/trimmed/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(read)

    output:
    tuple val(sample_id), path("${sample_id}.fastq.qtrim.gz"), emit: reads_pass

    script:
    """
    TRIMMOMATIC_JAR=\$(ls -1 /opt/conda/share/trimmomatic-*/trimmomatic*.jar 2>/dev/null | awk 'NR==1{print; exit}')
    ADAPTER_DIR=\$(ls -d /opt/conda/share/trimmomatic-*/adapters 2>/dev/null | awk 'NR==1{print; exit}')

    if [[ -z "\${TRIMMOMATIC_JAR}" || -z "\${ADAPTER_DIR}" ]]; then
        echo "Trimmomatic jar or adapters not found under /opt/conda/share" >&2
        exit 1
    fi

    java -jar "\${TRIMMOMATIC_JAR}" SE -threads ${task.cpus} -phred33 \\
        ${read} \\
        ${sample_id}.fastq.qtrim.gz \\
        ILLUMINACLIP:\${ADAPTER_DIR}/TruSeq3-SE.fa:2:30:10 \\
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    """
}

process QC_PE {
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

process QC_SE {
    tag "${sample_id}"
    cpus 4
    memory '8 GB'

    publishDir "${params.outdir}/qc/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(read)

    output:
    tuple val(sample_id), path(read), emit: reads_pass
    path '*_fastqc.html'
    path '*_fastqc.zip'

    script:
    """
    fastqc -t ${task.cpus} -o . ${read}
    """
}

process MAPPING_PE {
    tag "${sample_id}"
    cpus 16
    memory '64 GB'

    publishDir "${params.outdir}/mapping/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_sorted.bam"), emit: bam
    path "${sample_id}_sorted.bam.bai"
    path "${sample_id}_flagstat.txt"

    script:
    """
    set -euo pipefail

    # Extract HISAT2 index into the process workdir (avoid /mnt/genome caching).
    INDEX_TAR="./${params.hisat2_record}.tar.gz"
    INDEX_DIR="./${params.hisat2_record}"
    INDEX_BASE="${INDEX_DIR}/genome"

    if [ ! -e "${INDEX_BASE}.1.ht2" ]; then
        aws s3 cp "${params.hisat2_reference}/${params.hisat2_species}/${params.hisat2_record}.tar.gz" "${INDEX_TAR}"
        tar zxvf "${INDEX_TAR}"
    fi

    HISATOPT="-x ${INDEX_BASE} -p ${task.cpus}"

    # Write SAM to stdout and pipe directly to BAM+sort (avoid writing *_hits.sam).
    hisat2 \${HISATOPT} \\
        -1 ${sample_id}_1.fastq.P.qtrim.gz \\
        -2 ${sample_id}_2.fastq.P.qtrim.gz \\
        -S - \\
    | samtools view -@ ${task.cpus} -bS - \\
    | samtools sort -@ ${task.cpus} -o ${sample_id}_sorted.bam -
    samtools index ${sample_id}_sorted.bam
    samtools flagstat ${sample_id}_sorted.bam > ${sample_id}_flagstat.txt
    """
}

process MAPPING_SE {
    tag "${sample_id}"
    cpus 16
    memory '64 GB'

    publishDir "${params.outdir}/mapping/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(read)

    output:
    tuple val(sample_id), path("${sample_id}_sorted.bam"), emit: bam
    path "${sample_id}_sorted.bam.bai"
    path "${sample_id}_flagstat.txt"

    script:
    """
    set -euo pipefail

    # Extract HISAT2 index into the process workdir (avoid /mnt/genome caching).
    INDEX_TAR="./${params.hisat2_record}.tar.gz"
    INDEX_DIR="./${params.hisat2_record}"
    INDEX_BASE="${INDEX_DIR}/genome"

    if [ ! -e "${INDEX_BASE}.1.ht2" ]; then
        aws s3 cp "${params.hisat2_reference}/${params.hisat2_species}/${params.hisat2_record}.tar.gz" "${INDEX_TAR}"
        tar zxvf "${INDEX_TAR}"
    fi

    HISATOPT="-x ${INDEX_BASE} -p ${task.cpus}"

    # Write SAM to stdout and pipe directly to BAM+sort (avoid writing *_hits.sam).
    hisat2 \${HISATOPT} \\
        -U ${sample_id}.fastq.qtrim.gz \\
        -S - \\
    | samtools view -@ ${task.cpus} -bS - \\
    | samtools sort -@ ${task.cpus} -o ${sample_id}_sorted.bam -
    samtools index ${sample_id}_sorted.bam
    samtools flagstat ${sample_id}_sorted.bam > ${sample_id}_flagstat.txt
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
    gtf_file = channel.fromPath(params.gtf).first()

    def rt = (params.read_type ?: 'PE').toUpperCase()

    if (rt == 'SE') {
        reads_se_ch = Channel.fromPath("${params.fastq_dir}/*.fastq.gz")
            .filter { f ->
                def name = f.getName()
                !(name.endsWith('_1.fastq.gz') || name.endsWith('_2.fastq.gz'))
            }
            .map { f ->
                def sample_id = f.getName().replaceFirst(/\.fastq\.gz$/, '')
                tuple(sample_id, f)
            }

        TRIMMING_SE(reads_se_ch)
        QC_SE(TRIMMING_SE.out.reads_pass)
        MAPPING_SE(QC_SE.out.reads_pass)

        count_in = MAPPING_SE.out.bam
            .map { it -> it[1] }
            .collect()
            .combine(gtf_file)

        COUNT(count_in)
    } else {
        reads_pe_ch = Channel.fromFilePairs("${params.fastq_dir}/*_{1,2}.fastq.gz")

        TRIMMING_PE(reads_pe_ch)
        QC_PE(TRIMMING_PE.out.reads_pass)
        MAPPING_PE(QC_PE.out.reads_pass)

        count_in = MAPPING_PE.out.bam
            .map { it -> it[1] }
            .collect()
            .combine(gtf_file)

        COUNT(count_in)
    }
}
