params.fasta = "s3://ci-genome-tokyo/references/mm39/bwa_index/GRCm39.fa"
params.outdir = "s3://ci-genome-tokyo/references/mm39/bwa_index/"

process prepare_index {
    container 'quay.io/biocontainers/bwa:0.7.17--h5bf99c6_8'
    publishDir "${params.outdir}", mode: 'copy'
    cpus 4
    memory '32 GB'

    input:
    path fasta

    output:
    path "genome.*"

    script:
    """
    echo "Input file: ${fasta}"
    ls -lL ${fasta}
    head -n 2 ${fasta}
    echo "Starting bwa index..."
    bwa index -a bwtsw -p genome ${fasta}
    """
}

workflow {
    fasta_ch = Channel.fromPath(params.fasta)
    prepare_index(fasta_ch)
}
