nextflow.enable.dsl = 2

def ticket    = params.input_path.split('/')[0]
def fastq_dir = "${params.s3_bucket}/data/${params.input_path}"
def outdir    = "${params.s3_bucket}/data/${ticket}"

def hisat2_index_url = "${params.ref_base}/${params.hisat2_species}/${params.hisat2_record}.tar.gz"
def rsem_index_url   = "${params.ref_base}/${params.rsem_species}/rsem_${params.rsem_record}.tar.gz"

process TRIMMING_PE {
    label 'qc'
    tag "${sample_id}"
    cpus 4
    memory '8 GB'

    publishDir "${outdir}/trimmed", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_*.fastq.P.qtrim.gz"), emit: reads_pass
    path "${sample_id}_1.fastq.U.qtrim.gz", optional: true
    path "${sample_id}_2.fastq.U.qtrim.gz", optional: true

    script:
    """
    set -euo pipefail
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
    label 'qc'
    tag "${sample_id}"
    cpus 4
    memory '8 GB'

    publishDir "${outdir}/trimmed", mode: 'copy'

    input:
    tuple val(sample_id), path(read)

    output:
    tuple val(sample_id), path("${sample_id}.fastq.qtrim.gz"), emit: reads_pass

    script:
    """
    set -euo pipefail
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
    label 'qc'
    tag "${sample_id}"
    cpus 4
    memory '8 GB'

    publishDir "${outdir}/qc/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path(reads), emit: reads_pass
    path '*_fastqc.html'
    path '*_fastqc.zip'

    script:
    """
    set -euo pipefail
    fastqc -t ${task.cpus} -o . ${reads[0]} ${reads[1]}
    """
}

process QC_SE {
    label 'qc'
    tag "${sample_id}"
    cpus 4
    memory '8 GB'

    publishDir "${outdir}/qc/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(read)

    output:
    tuple val(sample_id), path(read), emit: reads_pass
    path '*_fastqc.html'
    path '*_fastqc.zip'

    script:
    """
    set -euo pipefail
    fastqc -t ${task.cpus} -o . ${read}
    """
}

process MAPPING_PE {
    label 'mapping'
    tag "${sample_id}"
    cpus 8
    memory { task.attempt <= 1 ? '24 GB' : '32 GB' }

    publishDir "${outdir}/mapping/${sample_id}", mode: 'copy'

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
    INDEX_BASE="\${INDEX_DIR}/genome"

    aws s3 cp "${hisat2_index_url}" "\${INDEX_TAR}"
    tar zxf "\${INDEX_TAR}"

    HISATOPT="-x \${INDEX_BASE} -p ${task.cpus}"

    hisat2 \${HISATOPT} \\
        -1 ${reads[0]} \\
        -2 ${reads[1]} \\
    | samtools sort -@ ${task.cpus} -O bam -o ${sample_id}_sorted.bam
    samtools index ${sample_id}_sorted.bam
    samtools flagstat ${sample_id}_sorted.bam > ${sample_id}_flagstat.txt
    """
}

process MAPPING_SE {
    label 'mapping'
    tag "${sample_id}"
    cpus 8
    memory { task.attempt <= 1 ? '24 GB' : '32 GB' }

    publishDir "${outdir}/mapping/${sample_id}", mode: 'copy'

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
    INDEX_BASE="\${INDEX_DIR}/genome"

    aws s3 cp "${hisat2_index_url}" "\${INDEX_TAR}"
    tar zxf "\${INDEX_TAR}"

    HISATOPT="-x \${INDEX_BASE} -p ${task.cpus}"

    hisat2 \${HISATOPT} \\
        -U ${read} \\
    | samtools sort -@ ${task.cpus} -O bam -o ${sample_id}_sorted.bam
    samtools index ${sample_id}_sorted.bam
    samtools flagstat ${sample_id}_sorted.bam > ${sample_id}_flagstat.txt
    """
}

process COUNT_PE {
    label 'count'
    tag "${sample_id}"
    cpus 8
    memory { task.attempt <= 1 ? '24 GB' : '32 GB' }

    publishDir "${outdir}/RSEM", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    path "${sample_id}.genes.results"
    path "${sample_id}.isoforms.results"
    path "${sample_id}.transcript_sorted.bam"
    path "${sample_id}.transcript_sorted.bam.bai"

    script:
    """
    set -euo pipefail

    # RSEM reference (rsem_run.sh: s3://.../ref/\${SPECIES}/rsem_\${RECORD}.tar.gz -> ./rsem/\${RECORD}/)
    RSEM_TAR="./rsem_${params.rsem_record}.tar.gz"
    RSEM_REF="./rsem/${params.rsem_record}"

    aws s3 cp "${rsem_index_url}" "\${RSEM_TAR}"
    tar zxf "\${RSEM_TAR}"

    rsem-calculate-expression -p ${task.cpus} \\
        --paired-end ${reads[0]} ${reads[1]} \\
        --bowtie2 \\
        "\${RSEM_REF}" \\
        ${sample_id}

    samtools sort -@ ${task.cpus} ${sample_id}.transcript.bam -o ${sample_id}.transcript_sorted.bam
    samtools index ${sample_id}.transcript_sorted.bam
    """
}

process COUNT_SE {
    label 'count'
    tag "${sample_id}"
    cpus 8
    memory { task.attempt <= 1 ? '24 GB' : '32 GB' }

    publishDir "${outdir}/RSEM", mode: 'copy'

    input:
    tuple val(sample_id), path(read)

    output:
    path "${sample_id}.genes.results"
    path "${sample_id}.isoforms.results"
    path "${sample_id}.transcript_sorted.bam"
    path "${sample_id}.transcript_sorted.bam.bai"

    script:
    """
    set -euo pipefail

    RSEM_TAR="./rsem_${params.rsem_record}.tar.gz"
    RSEM_REF="./rsem/${params.rsem_record}"

    aws s3 cp "${rsem_index_url}" "\${RSEM_TAR}"
    tar zxf "\${RSEM_TAR}"

    rsem-calculate-expression -p ${task.cpus} \\
        ${read} \\
        --bowtie2 \\
        "\${RSEM_REF}" \\
        ${sample_id}

    samtools sort -@ ${task.cpus} ${sample_id}.transcript.bam -o ${sample_id}.transcript_sorted.bam
    samtools index ${sample_id}.transcript_sorted.bam
    """
}

workflow {
    def rt = (params.read_type ?: 'PE').toUpperCase()

    if (rt == 'SE') {
        reads_se_ch = Channel.fromPath("${fastq_dir}/*.fastq.gz")
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
        MAPPING_SE(TRIMMING_SE.out.reads_pass)
        COUNT_SE(TRIMMING_SE.out.reads_pass)
    } else {
        reads_pe_ch = Channel.fromFilePairs("${fastq_dir}/*_{1,2}.fastq.gz")

        TRIMMING_PE(reads_pe_ch)
        QC_PE(TRIMMING_PE.out.reads_pass)
        MAPPING_PE(TRIMMING_PE.out.reads_pass)
        COUNT_PE(TRIMMING_PE.out.reads_pass)
    }
}
