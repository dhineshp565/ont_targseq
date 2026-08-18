#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { make_csv } from './modules/local/make_csv'
include { merge_fastq } from './modules/local/merge_fastq'
include { porechop } from './modules/local/porechop'
include { minimap2 } from './modules/local/minimap2'
include { splitbam } from './modules/local/splitbam'
include { medaka } from './modules/local/medaka'
include { multiqc } from './modules/local/multiqc'
include { kraken2 } from './modules/local/kraken2'
include { krona_kraken } from './modules/local/krona_kraken'
include { make_report } from './modules/local/make_report'
include { htmltopdf } from './modules/local/htmltopdf'
include { blast_cons } from './modules/local/blast_cons'
include { orfipy } from './modules/local/orfipy'
include { typing } from './modules/local/typing'
include { make_limsfile } from './modules/local/make_limsfile'
include { mafft } from './modules/local/mafft'
include { iqtree } from './modules/local/iqtree'
include { ggtree } from './modules/local/ggtree'
include { extract_mapped_ref } from './modules/local/extract_mapped_ref'
include { mapped_ref_bed } from './modules/local/mapped_ref_bed'
include { bedtools } from './modules/local/bedtools'
include { igvreports } from './modules/local/igvreports'
include { seq_length } from './modules/local/seq_length'


// Include dehost subworkflow
include { DEHOST } from './subworkflows/dehost.nf'
include { AMPLICONS } from './subworkflows/amplicons.nf'
include { METAGENOMICS } from './subworkflows/metagenomics.nf'
include { QCREADS } from './subworkflows/qcreads.nf'
include { IGVREPORT } from './subworkflows/igvreport.nf'

workflow {
	// Merge fastq files and qc filtering
    QCREADS(params.input, params.qscore, params.trim_barcodes)

    // Dehost logic
    if (params.dehost) {
        host_ref = file(params.host_db)
        DEHOST(QCREADS.out.reads, host_ref)
        reads_for_alignment = DEHOST.out.dehosted_reads
    } else {
        reads_for_alignment = QCREADS.out.reads
    }
	
	reference=file(params.reference)
	primerbed=file("${baseDir}/primer.bed")
	software_version_file=file("${baseDir}/software_version.tsv")
	AMPLICONS(reads_for_alignment,reference,primerbed,params.read_count_threshold,params.consensus_mode,params.qscore)
	
	if (params.dehost) {
		stats=DEHOST.out.stats
	} else {
		stats=AMPLICONS.out.stats
	} 
	kraken2(reads_for_alignment, params.kraken_db,params.kraken_confidence)
	//bracken(kraken2.out, params.kraken_db)
	krona_kraken(kraken2.out.kraken2_raw.collect())
	
	// qc report using split bam out put
	idxstats=AMPLICONS.out.idxstats
	nanoqc=QCREADS.out.read_stats	
	multiqc(nanoqc.mix(idxstats,stats).collect())
	dbdir=file("${baseDir}/targseq")
	
	typing(AMPLICONS.out.consensus,dbdir)
	make_limsfile(typing.out.withseq.collect(),software_version_file)
	
	blast_cons(AMPLICONS.out.consensus,params.blastdb_path,params.blastdb_name)

	refdir="${baseDir}/reference_sequences"
	mafft(QCREADS.out.csv,typing.out.filtered.collect(),refdir)
	iqtree(mafft.out.collect())
	ggtree(iqtree.out.collect())
	orfipy(AMPLICONS.out.consensus)

	// Pair consensus and ORF files by sample name
	paired_consensus_orf = AMPLICONS.out.consensus.join(orfipy.out.orf)
    .map { sample, consensus, orf -> tuple(sample, consensus, orf) }
	seq_length(paired_consensus_orf)
	
	
	//generate report


	IGVREPORT(reference,AMPLICONS.out.bam)


	rmd_file=file("${baseDir}/targseq_rmdfile.Rmd")
	rmdfile_case=file("${baseDir}/targseq_rmdfile_case.Rmd")

	make_report(QCREADS.out.csv,
				krona_kraken.out.raw,
				AMPLICONS.out.mapped.collect(),
				typing.out.filtered.collect(),
				typing.out.withseq.collect(),
				blast_cons.out.blast_formatted.collect(),
				ggtree.out.png,
				rmd_file,
				IGVREPORT.out,
				orfipy.out.orf_only.collect(),
				seq_length.out.collect(),
				rmdfile_case)
	// htmltopdf(make_report.out.pdf)

}

