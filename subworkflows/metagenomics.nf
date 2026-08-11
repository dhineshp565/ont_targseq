#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { kraken2 } from '../modules/local/kraken2.nf'
include { krona_kraken } from '../modules/local/krona_kraken.nf'

process bracken {
publishDir "${params.out_dir}/bracken/",mode:"copy"
label "low"

input:
tuple val(SampleName), path(kraken_output),path(kraken_report)
path(kraken_db)

output:
tuple val(SampleName), path("${SampleName}_bracken.tsv")

script:
"""


# 1. Run Bracken
bracken -d ${kraken_db} -i ${kraken_report} -o bracken_raw.tsv -l G

if [[ ! -s bracken_raw.tsv ]]; then
    echo -e "name\\ttaxonomy_id\\ttaxonomy_lvl\\tkraken_assigned_reads\\tadded_reads\\tnew_est_reads\\tfraction_total_reads\\nNo microbial reads\\t0\\t\$BRACKEN_LEVEL\\t0\\t0\\t0\\t0.00" > bracken_raw.tsv
fi

# 2. Sort by relative abundance (column 7) and keep header
head -n 1 bracken_raw.tsv > ${SampleName}_bracken.tsv
tail -n +2 bracken_raw.tsv | sort -t\$'\t' -k7,7nr >> ${SampleName}_bracken.tsv
"""
}


process krakenuniq {


    publishDir "${params.out_dir}/krakenuniq/",mode:"copy"
	label "high"
	input:
	tuple val(SampleName),path (SamplePath)
	path(db_path)
    
	
	output:
	tuple val(SampleName), path("${SampleName}_krakenuniq.csv"), emit: krakenuniq_output
    tuple val(SampleName), path("${SampleName}_krakenuniq_report.csv"), emit: krakenuniq_report
	
	
	script:
	"""
	krakenuniq --db $db_path --threads 50 --preload-size 200G --report-file ${SampleName}_krakenuniq_report.tsv --output ${SampleName}_krakenuniq.tsv ${SamplePath}
	
	"""
}


process extract_reads {
    publishDir "${params.out_dir}/binned_reads", mode: 'copy'
    label "high"

    input:
    tuple val(SampleName), path(fastq), path(bracken), path(kraken_output), path(kraken_report)

    output:
    tuple val(SampleName), path("*.fastq"), path("${SampleName}_taxids.txt"), path("${SampleName}_taxinfo.txt")

    script:
    """
    # Get TaxIDs and names from Bracken output (skip header, get columns 2 and 1)
    tail -n +2 ${bracken} | awk -F'\t' '{print \$2"\t"\$1}' > ${SampleName}_taxinfo.txt
    
    # Extract just taxIDs for the loop
    tail -n +2 ${bracken} | cut -f 2 > ${SampleName}_taxids.txt
    
    while read -r taxid; do
        extract_kraken_reads.py \\
            -k ${kraken_output} \\
            -s ${fastq} \\
            -o ${SampleName}_\${taxid}.fastq \\
            -t \${taxid} \\
            --include-children \\
            --report ${kraken_report}
    done < ${SampleName}_taxids.txt
    """
}

process megahit {
    publishDir "${params.out_dir}/megahit_metagenomes/", mode: "copy"
    label "high"

    input:
    tuple val(SampleName), path(fastqs), path(taxids), path(taxinfo)

    output:
    tuple val(SampleName), path("${SampleName}_metagenomes.fasta")

    script:
    """
    touch ${SampleName}_metagenomes.fasta
    
    # Loop through taxids using cat - create combined assembly file
    for taxid in \$(cat ${taxids}); do
      
        megahit -r "${SampleName}_\${taxid}.fastq" -o ${SampleName}_\${taxid}_assembly --k-list 41,61,81,99 --prune-level 1 --min-contig-len 200

        # Process output - append to combined assembly file
        if [ -s "${SampleName}_\${taxid}_assembly/final.contigs.fa" ]; then
            # Rename headers to include taxID and append
             sed -E "s/>k[0-9]+_/>${SampleName}_\${taxid}_contig_/g" "${SampleName}_\${taxid}_assembly/final.contigs.fa" >> ${SampleName}_metagenomes.fasta
        fi
    done
    
    # If no contigs were assembled, add placeholder
    if [ ! -s ${SampleName}_metagenomes.fasta ]; then
        echo ">${SampleName}_no_assembly" > ${SampleName}_metagenomes.fasta
        echo "NNNN" >> ${SampleName}_metagenomes.fasta
    fi
    """
}

process refseq_masher {
    publishDir "${params.out_dir}/refseq_masher/", mode: "copy"
    label "high"

    input:
    tuple val(SampleName), path(assembly)

    output:
    path("${SampleName}_refseqmasher.tsv")

    script:
    """
    # Check if assembly file has content and valid FASTA records
    if [ -s ${assembly} ] && grep -q "^>" ${assembly}; then
        refseq_masher -vv matches --top-n-results 50 ${assembly} > ${SampleName}_refseqmasher.tsv
    else
        # Create empty result file with header if no assembly
        echo -e "sample\ttop_taxonomy_name\tdistance\tpvalue\tmatching\tfull_taxonomy\ttaxonomic_subspecies\ttaxonomic_species\ttaxonomic_genus\ttaxonomic_family\ttaxonomic_order\ttaxonomic_class\ttaxonomic_phylum\ttaxonomic_superkingdom\tsubspecies\tspecies\tgenus\tfamily\torder\tclass\tphylum\tsuperkingdom" > ${SampleName}_refseqmasher.tsv
        echo -e "${SampleName}\tNO_ASSEMBLY\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A" >> ${SampleName}_refseqmasher.tsv
    fi
    """
}

//performs blast of the consensus sequences
process blast_cons {
	publishDir "${params.out_dir}/blast/",mode:"copy"

	label "high"
	input:
	tuple val(SampleName),path(consensus)
	path (blastdb_path)
	val(blastdb_name)

	output:
	tuple val(SampleName), path ("${SampleName}_blast.tsv"),emit:blast
	script:
	"""
	# Write header first
	echo -e "queryid\tsubject_id\talignment length\tquery_coverage\t%identity\tevalue\tbitscore\tstaxids\tsscinames\tscomnames\tstitle" > "${SampleName}_blast.tsv"

	# Run BLAST and append results directly to report
	blastn -task megablast -db ${blastdb_path}/${blastdb_name} -query ${consensus} -outfmt "6 qseqid sseqid length qcovs pident evalue bitscore staxids ssciname scomnames stitle" -max_target_seqs 5 >> "${SampleName}_blast.tsv" || true

	# If only header was written (no BLAST hits), add a "NaN" placeholder row
	if [[ \$(wc -l < "${SampleName}_blast.tsv") -le 1 ]]; then
	    echo -e "NaN\tNaN\tNaN\tNaN\tNaN\tNaN\tNaN\tNaN\tNaN\tNaN\tNaN" >> "${SampleName}_blast.tsv"
	fi
	"""

}

process metagenomics_summary {
    publishDir "${params.out_dir}/merged_kraken_blast/", mode: "copy"
    label "high"

    input:
    tuple val(SampleName),path(bracken),path(blast_report)

    output:
    tuple val(SampleName), path("${SampleName}_blast_best_hits.tsv"), path("${SampleName}_kracken_blast_summary.tsv")

    script:
    """
    metagenomics_summary.py ${bracken} ${blast_report} ${SampleName}
    """
}


workflow METAGENOMICS {
    take:
	reads // Channel: tuple val(SampleName), path(fastq)
	kraken_db  // Path: kraken database
    blastdb_path // Path to blast database
    blastdb_name // Name of blast database
	kraken_confidence // Confidence threshold for kraken2
	main:
	// Run kraken2
	kraken2(reads, kraken_db, kraken_confidence)
	
	// Prepare bracken input by joining kraken outputs
	bracken_input = kraken2.out.kraken_output
	    .join(kraken2.out.kraken_report)
	    .map { sample, krakenraw, krakenReport -> tuple(sample, krakenraw, krakenReport) }
	
	// Run bracken
	bracken(bracken_input, kraken_db)
	
	// Prepare extract_reads input by joining all required channels
	extract_reads_input = reads
	    .join(kraken2.out.kraken_output)
	    .join(kraken2.out.kraken_report)
	    .join(bracken.out)
	    .map { sample, fastq, krakenCsv, krakenReport, brackenTsv -> 
	        tuple(sample, fastq, brackenTsv, krakenCsv, krakenReport) 
	    }
	
	extract_reads(extract_reads_input)
	megahit(extract_reads.out)
    //medaka_input = reads.join(megahit.out).map{sample,fastq,draft_fasta -> tuple (sample,fastq,draft_fasta)}
   // medaka(medaka_input)
    krona_kraken(kraken2.out.kraken_report.map{ _sample, file -> file }.collect())
    //refseq_masher(megahit.out)
    blast_cons(megahit.out,blastdb_path,blastdb_name)
    summary_input = bracken.out.join(blast_cons.out).map{ sample, brackentsv, blastReport -> tuple(sample, brackentsv, blastReport) }
    metagenomics_summary(summary_input)
	
	emit:
	metagenomes = megahit.out.map{_sample,consensus -> consensus}
	kraken_output = kraken2.out.kraken_output
	kraken_report = kraken2.out.kraken_report
	bracken_output = bracken.out
    krona_output = krona_kraken.out
    //refmash=refseq_masher.out
    blast_report = blast_cons.out.blast
    merged_report = metagenomics_summary.out

	
}