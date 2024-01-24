#!/usr/bin/env nextflow

if( nextflow.version.matches(">= 20.07.1") ){
    nextflow.enable.dsl=2
} else {
    // Support lower version of nextflow
    nextflow.preview.dsl=2
}


process Alignment{
    publishDir "${params.out_dir}/alignment"

    input:
        path ch_fastq
        path ch_reference 

    output:
        path "reads2ref.sorted.bam", emit: ch_bam
        path "reads2ref.sorted.bam.bai", emit: ch_bai
        path "reference.fa", emit: ch_ref
        path "reference.fa.fai", emit: ch_ref_fai

    script:
    """
    ln -s ${ch_reference} reference.fa
    samtools faidx reference.fa

    if [ ${params.tgs_type} == "ont" ]; then
        minimap2 -ax map-ont -t ${params.nthreads} reference.fa $ch_fastq > reads2ref.sam
    elif [ ${params.tgs_type} == "hifi" ]; then
        minimap2 -ax map-pb -t ${params.nthreads} reference.fa $ch_fastq > reads2ref.sam
    else
        echo "Error: tgs_type must be either 'ont' or 'hifi'"
        exit 1
    fi
    samtools view -bS -@ ${params.nthreads} reads2ref.sam -o reads2ref.bam
    samtools sort -@ ${params.nthreads} reads2ref.bam -o reads2ref.sorted.bam
    samtools index -@ ${params.nthreads} reads2ref.sorted.bam
    """
}

process Clair3Calling{
    publishDir "${params.out_dir}/variant_calling"

    input:
    path ch_bam
    path ch_bai
    path ch_ref
    path ch_ref_fai

    output:
    path "final.vcf", emit: ch_vcf

    script:
    """
    if [[ ${params.tgs_type} == "ont" ]] || [[ ${params.tgs_type} == "hifi" ]]; then
        run_clair3.sh --bam_fn=${ch_bam} --ref_fn=${ch_ref} --threads=${params.nthreads} --platform=${params.tgs_type} --model_path="/opt/models/${params.tgs_type}" --output="variant_output"
        gzip -fdc variant_output/merge_output.vcf.gz > final.vcf

    else
        echo "Error: tgs_type must be either 'ont' or 'hifi' for clair3"
        exit 1
    fi
    """
}

process PepperDeepVariantCalling{
    publishDir "${params.out_dir}/variant_calling"

    input:
    path ch_bam
    path ch_bai
    path ch_ref
    path ch_ref_fai

    output:
    path "final.vcf", emit: ch_vcf

    script:
    """
    if [ ${params.tgs_type} == "ont" ]; then
        run_pepper_margin_deepvariant call_variant -b ${ch_bam} -f ${ch_ref} -o variant_output -t ${params.nthreads} --ont_r9_guppy5_sup
        gzip -fdc variant_output/PEPPER_MARGIN_DEEPVARIANT_FINAL_OUTPUT.vcf.gz > final.vcf
    elif [ ${params.tgs_type} == "hifi" ]; then
        run_pepper_margin_deepvariant call_variant -b ${ch_bam} -f ${ch_ref} -o variant_output -t ${params.nthreads} --hifi
        gzip -fdc variant_output/PEPPER_MARGIN_DEEPVARIANT_FINAL_OUTPUT.vcf.gz > final.vcf
    else
        echo "Error: tgs_type must be either 'ont' or 'hifi' for pepper_margin_deepvariant"
        exit 1
    fi
    """
}

process NanoSNPCalling{
    publishDir "${params.out_dir}/variant_calling"

    input:
    path ch_bam
    path ch_bai
    path ch_ref
    path ch_ref_fai

    output:
    path "final.vcf", emit: ch_vcf

    script:
    """
    if [ ${params.tgs_type} == "ont" ] && [[ ! ${params.usecontig} ]]; then
        echo "run1"
        run_caller.sh -b ${ch_bam} -f ${ch_ref} -t ${params.nthreads} -c ${params.coverage} -o variant_output
        mv variant_output/merge.vcf final.vcf
    elif [ ${params.tgs_type} == "ont" ] && [ ${params.usecontig} ]; then
        echo "run2"
        run_caller.sh -b ${ch_bam} -f ${ch_ref} -t ${params.nthreads} -c ${params.coverage} -g -o variant_output
        mv variant_output/merge.vcf final.vcf
    else
        echo "Error: tgs_type must be 'ont' for nanosnp"
        exit 1
    fi
    """
}


// process VariantCalling{
//     publishDir "${params.out_dir}/variant_calling"

//     input:
//     path ch_bam
//     path ch_reference

//     output:
//     path "final.vcf", emit: ch_vcf

//     script:
//     """
//     if [ ${params.snp_caller} == "clair3" ]; then
//         if [[ ${params.tgs_type} == "ont" ]] || [[ ${params.tgs_type} == "hifi" ]]; then
//             run_clair3.sh --bam_fn=${ch_bam} --ref_fn=${ch_reference} --threads=${params.nthreads} --platform=${params.tgs_type} --model_path="/opt/models/${params.tgs_type}" --output="variant_output"
//             gzip -fdc variant_output/merge_output.vcf.gz > final.vcf

//         else
//             echo "Error: tgs_type must be either 'ont' or 'hifi' for clair3"
//             exit 1
//         fi
//     elif [ ${params.snp_caller} == "pepper" ]; then
//         if [ ${params.tgs_type} == "ont" ]; then
//             run_pepper_margin_deepvariant call_variant -b ${ch_bam} -f ${ch_reference} -o variant_output -t ${params.nthreads} --ont_r9_guppy5_sup
//             mv variant_output/PEPPER_MARGIN_DEEPVARIANT_OUTPUT.vcf final.vcf
//         elif [ ${params.tgs_type} == "hifi" ]; then
//             run_pepper_margin_deepvariant call_variant -b ${ch_bam} -f ${ch_reference} -o variant_output -t ${params.nthreads} --hifi
//             mv variant_output/PEPPER_MARGIN_DEEPVARIANT_OUTPUT.vcf final.vcf
//         else
//             echo "Error: tgs_type must be either 'ont' or 'hifi' for pepper_margin_deepvariant"
//             exit 1
//         fi
//     elif [ ${params.snp_caller} == "nanosnp" ]; then
//         run_caller.sh -b ${ch_bam} -f ${ch_reference} -t ${params.nthreads} -c ${params.coverage} -o variant_output
//         mv variant_output/merge.vcf final.vcf
//     else
//         echo "Error: snp_caller must be either 'clair3', 'pepper', or 'nanosnp'"
//         exit 1
//     fi
//     """
// }


process Evaluation{
    publishDir "${params.out_dir}/evaluation"

    input:
    path ch_vcf
    path ch_ref
    path ch_ref_fai

    output:
    path "hap_output/happy.output.summary.csv", emit: ch_happy_summary

    script:
    """
    /opt/hap.py/bin/hap.py ${params.true_vcf} ${ch_vcf} -f ${params.true_bed} -r ${ch_ref} -o "hap_output/happy.output" --pass-only --engine=vcfeval --threads=${params.nthreads}
    """
}

ch_fastq=Channel.fromPath(params.fastq)
ch_reference=Channel.fromPath(params.reference)


workflow {

    (ch_bam,ch_bai,ch_ref,ch_ref_fai)=Alignment(ch_fastq, ch_reference)
    if ( params.snp_caller == "clair3" ){
        ch_vcf=Clair3Calling(ch_bam, ch_bai, ch_ref,ch_ref_fai)
    }
    else if ( params.snp_caller == "pepper" ){
        ch_vcf=PepperDeepVariantCalling(ch_bam, ch_bai, ch_ref,ch_ref_fai)
    }
    else if ( params.snp_caller == "nanosnp" ){
        ch_vcf=NanoSNPCalling(ch_bam, ch_bai, ch_ref,ch_ref_fai)
    }
    else{
        echo "Error: snp_caller must be either 'clair3', 'pepper', or 'nanosnp'"
        exit 1
    }
}
