#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process merge_filter_fastq {
	publishDir "${params.out_dir}/merged"
	label "low"
	input:
	tuple val(SampleName),path(SamplePath)
	val (qscore)
	output:
	tuple val(SampleName),path("${SampleName}_filtered.{fastq,fastq.gz}"),emit:reads

	script:
	"""
	count=\$(ls -1 ${SamplePath}/*.gz 2>/dev/null | wc -l)
	
	
		if [[ "\${count}" != "0" ]]
		then
			cat ${SamplePath}/*.fastq.gz > ${SampleName}.fastq.gz
			nanoq -i ${SampleName}.fastq.gz -q ${qscore} -o ${SampleName}_filtered.fastq.gz
		
		else
			count=\$(ls -1 ${SamplePath}/*.fastq 2>/dev/null | wc -l)
			if [[ "\${count}" != "0" ]]
			then
				cat ${SamplePath}/*.fastq| gzip > ${SampleName}.fastq.gz
				nanoq -i ${SampleName}.fastq.gz -q ${qscore} -o ${SampleName}_filtered.fastq.gz
			fi
		fi
	"""
}
