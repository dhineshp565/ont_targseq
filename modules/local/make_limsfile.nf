#!/usr/bin/env nextflow

process make_limsfile {
    publishDir "${params.out_dir}/LIMS/",mode:"copy"
    label "low"
    
    input:
    path (withseq)
    path (octoflu_out)
    path (software_version_file)
    
    output:
    path("LIMSfile_*.tsv")
    path("octoflu_*.tsv")
    
    script:
    """
    date=\$(date '+%Y-%m-%d_%H-%M-%S')
    
    awk 'FNR==1 && NR!=1 { while (/^#F/) getline; } 1 {print}' *_withseq.tsv > LIMSfile.tsv

    sed -i 's/_consensus//g' LIMSfile.tsv

    cat ${software_version_file} LIMSfile.tsv > LIMSfile_\${date}.tsv

    cat ${octoflu_out} > octoflu_\${date}.tsv
     sed -i 's/_consensus//g' octoflu_\${date}.tsv
    sed -i "1i SEQUENCE\tTYPE\tUS CLADE\tGLOBAL CLADE" octoflu_\${date}.tsv
    """
}
