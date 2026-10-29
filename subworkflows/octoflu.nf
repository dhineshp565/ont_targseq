#!/usr/bin/env nextflow

process octoflu_nf {
   label "high"
    publishDir "${params.out_dir}/octoflu/" ,mode: 'copy'
    input:
    tuple(val(SampleName), path (consensus))

    output:
    path "*_output"
    path "${SampleName}_octoflu.txt",emit: clade_result

    script:
    """
    if grep -q 'InfA' '${consensus}'; then
        octoFLU.sh '${consensus}'
        sed -E 's/_[HN][0-9]+_consensus//g' *_Final_Output.txt > '${SampleName}_octoflu.txt'
    else
        touch '${SampleName}_octoflu.txt'
        mkdir -p ${SampleName}_output
    fi
   

    """
}


workflow OCTOFLU {
    take: 
        consensus

    main:
        infa_fasta_ch=consensus.map{ sample,consensus_file -> tuple(sample,consensus_file) }
        //.filter {_sample,consensus_file -> consensus_file.splitFasta(record: [header: true])
        //.any { record -> record.header.contains('InfA')}}    
        //.set {infa_fasta_ch}                   
                        
         octoflu_nf(infa_fasta_ch)
    
    emit:
    clade_result=octoflu_nf.out.clade_result
}

