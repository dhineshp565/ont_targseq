#!/usr/bin/env nextflow

process octoflu_nf {
    label 'high'

    publishDir "${params.out_dir}/octoflu/", mode: 'copy'

    input:
    tuple val(SampleName), path(consensus)
    path(octoflu_reference)

    output:
    path('*_output')
    path("${SampleName}_octoflu.txt"), emit: clade_result

    script:
 

    """
    if grep -q 'InfA' '${consensus}'; then
        octoFLU.sh '${consensus}' ${octoflu_reference}

        if ls *_Final_Output.txt >/dev/null 2>&1; then
            sed -E 's/_[HN][0-9]+_consensus//g' *_Final_Output.txt > '${SampleName}_octoflu.txt'
        else
            touch '${SampleName}_octoflu.txt'
        fi
    else
        touch '${SampleName}_octoflu.txt'
        mkdir -p '${SampleName}_output'
    fi
    """
}


workflow OCTOFLU {
    take:
    consensus
    octoflu_reference

    main:
    // Handle optional reference - map consensus with reference or empty
   
    octoflu_nf(consensus,octoflu_reference)
    

    emit:
    clade_result = octoflu_nf.out.clade_result
}