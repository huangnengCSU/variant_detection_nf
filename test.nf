#!/usr/bin/env nextflow

if( nextflow.version.matches(">= 20.07.1") ){
    nextflow.enable.dsl=2
} else {
    // Support lower version of nextflow
    nextflow.preview.dsl=2
}


process mapping{
    publishDir "${out_dir}"

    input:
        file raw_fastq
        file ref
        val nthreads
    
    output:
        path "align.sam", emit: ch_sam

    script:
    """
    minimap2 --version
    minimap2 -ax map-pb $ref $raw_fastq -t $nthreads > align.sam
    """
}

process samtools{
    publishDir "${out_dir}"

    input:
        path ch_sam
    
    output:
        path "align.sort.bam", emit: ch_bam
        path "align.sort.bam.bai", emit: ch_bai

    script:
    """
    samtools view -bS $ch_sam -o align.bam
    samtools sort align.bam -o align.sort.bam
    samtools index align.sort.bam
    """
}

raw_fastq=Channel.fromPath(params.raw_fastq)

ref=Channel.fromPath(params.ref)

out_dir=params.out_dir

nthreads=params.nthreads

workflow {
    ch_sam=mapping(raw_fastq, ref, nthreads)
    samtools(ch_sam)
}
