nextflow.enable.dsl=2

def genome_params = params.genomes[params.genome_id]

process TRIMMING {
    tag "${sample_id}"
    cpus 8
    memory '16 GB'
    
    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_trim_*.fastq.gz"), emit: trimmed_reads

    script:
    """
    fastp -i ${reads[0]} -I ${reads[1]} \\
          -o ${sample_id}_trim_1.fastq.gz -O ${sample_id}_trim_2.fastq.gz \\
          -3 -q 15 -n 0 -l 100 -w ${task.cpus} --detect_adapter_for_pe
    """
}

process MAPPING {
    tag "${sample_id}"
    cpus 16
    memory '64 GB'
    
    input:
    tuple val(sample_id), path(reads)
    // 修正ポイント：受け取ったファイル群を 'idx/' フォルダ内に展開する
    path "idx/*" 

    output:
    tuple val(sample_id), path("${sample_id}.rmdup.bam"), emit: bam
    path "${sample_id}.rmdup.bam.bai"

    script:
    // 修正ポイント：インデックスのパスを 'idx/genome' のように指定する
    """
    bwa mem -t ${task.cpus} -M idx/genome ${reads[0]} ${reads[1]} | \\
    samtools sort -@ ${task.cpus} -o ${sample_id}.mapped.bam -

    picard MarkDuplicates I=${sample_id}.mapped.bam O=${sample_id}.rmdup.bam \\
           M=${sample_id}.dup_metrics.txt REMOVE_DUPLICATES=true
    
    samtools index ${sample_id}.rmdup.bam
    """
}

process MAKE_TAG_DIR {
    tag "${sample_id}"
    cpus 4
    memory '16 GB'
    
    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}_TagDir"), emit: tagdir

    script:
    """
    makeTagDirectory ${sample_id}_TagDir -single ${bam}
    """
}

process DIFFERENTIAL_ANALYSIS {
    publishDir "${params.outdir}", mode: 'copy'
    cpus 8
    memory '32 GB'

    input:
    path tagdirs
    val gsize

    output:
    path "differential/*.txt"
    path "differential/*.bed"

    script:
    """
    mkdir -p differential
    
    # 実際のサンプル名 (WT_IP, WT_Inputなど) に合わせて実行 
    findPeaks WT_IP_TagDir -style histone -gsize ${gsize} -o auto -i WT_Input_TagDir
    findPeaks MFP2KO_IP_TagDir -style histone -gsize ${gsize} -o auto -i MFP2KO_Input_TagDir
    findPeaks BMAL1KO_IP_TagDir -style histone -gsize ${gsize} -o auto -i BMAL1KO_Input_TagDir

    # マージ
    mergePeaks WT_IP_TagDir/regions.txt \\
               MFP2KO_IP_TagDir/regions.txt \\
               BMAL1KO_IP_TagDir/regions.txt > differential/merged_peaks.txt

    # 比較
    getDifferentialPeaks differential/merged_peaks.txt WT_IP_TagDir MFP2KO_IP_TagDir -F 4.0 -P 0.0001 > differential/diff_WT_vs_MFP2KO.txt
    getDifferentialPeaks differential/merged_peaks.txt WT_IP_TagDir BMAL1KO_IP_TagDir -F 4.0 -P 0.0001 > differential/diff_WT_vs_BMAL1KO.txt

    # BED作成
    grep -v ^# differential/diff_WT_vs_MFP2KO.txt | perl -ne '\$"="\\t"; @a=split; \$a[6]--; print"@a[6-8]\\n"' > differential/diff_WT_vs_MFP2KO.bed
    grep -v ^# differential/diff_WT_vs_BMAL1KO.txt | perl -ne '\$"="\\t"; @a=split; \$a[6]--; print"@a[6-8]\\n"' > differential/diff_WT_vs_BMAL1KO.bed
    """
}

workflow {
    read_pairs_ch = Channel.fromFilePairs("${params.fastq_dir}/*_{1,2}.fastq.gz")
    
    // 修正ポイント：末尾に * を付けて全インデックスファイルを取得し、collect() で一つのセットにする
    bwa_idx_files = Channel.fromPath("s3://ci-genome-tokyo/references/mm39/bwa_index/*").collect()

    TRIMMING(read_pairs_ch)
    // 修正ポイント：インデックスファイルのセットを渡す
    MAPPING(TRIMMING.out.trimmed_reads, bwa_idx_files)

    MAKE_TAG_DIR(MAPPING.out.bam)

    all_tagdirs = MAKE_TAG_DIR.out.tagdir.map{ it[1] }.collect()
    DIFFERENTIAL_ANALYSIS(all_tagdirs, genome_params.gsize)
}

