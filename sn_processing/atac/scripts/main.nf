#!/usr/bin/env nextflow

nextflow.enable.dsl=2
IONICE = 'ionice -c2 -n7'

// Generic data
AUTOSOMAL_REFERENCES = ['hg19': (1..22).collect({it -> 'chr' + it}),
			'hg38': (1..22).collect({it -> 'chr' + it}),
			'rn5': (1..20).collect({it -> 'chr' + it}),
			'rn6': (1..20).collect({it -> 'chr' + it}),
			'mm9': (1..19).collect({it -> 'chr' + it}),
			'mm10': (1..19).collect({it -> 'chr' + it})
]

ORGANISMS = ['hg19': 'human', 
		'hg38': 'human',
		'rn5': 'rat',
		'rn6': 'rat',
		'mm9': 'mouse',
		'mm10': 'mouse']

libraries = params.libraries.keySet()

def make_excluded_regions_arg(genome) {
	return params.blacklist[genome].collect({'--excluded-region-file ' + it}).join(' ')
}


def is_chimeric(library) {
	return get_genome(library).size() > 1
}


def get_bwa_index(genome) {
	return params.bwa_index[genome]
}


def get_genome(library) {
	return params.libraries[library].genome
}


def get_tss(genome) {
	return params.tss[genome]
}


def get_organism(genome) {
	return ORGANISMS[genome]
}


def get_chrom_sizes(genome) {
	return params.chrom_sizes[genome]
}


def get_gene_bed(genome) {
	return params.gene_bed[genome]
}


def library_to_readgroups(library) {
	return params.libraries[library].readgroups.keySet()
}


def library_and_readgroup_to_fastqs(library, readgroup) {
	return params.libraries[library].readgroups[readgroup]
}


trim_reads = []
make_barcode_corrections_in = []
fastqc_reads = []

rg_number = 0
for (library in libraries) {
	readgroups = library_to_readgroups(library)
	for (readgroup in readgroups) {
		fastqs = library_and_readgroup_to_fastqs(library, readgroup)
		first_insert = fastqs['1']
		second_insert = fastqs['2']
		barcode = fastqs['index']
		rgid = file(barcode).name.replaceFirst('_L000_R2_001.fastq.gz', '')
		trim_reads << [library, rgid, readgroup, file(first_insert), file(second_insert)]
		make_barcode_corrections_in << [library, file(barcode)]
	}
	if (rg_number > 3) {
		rg_number = 0
	}
	readgroup = readgroups[rg_number]
	fastqs = library_and_readgroup_to_fastqs(library, readgroup)
	first_insert = fastqs['1']
	second_insert = fastqs['2']
	barcode = fastqs['index']
	fastqc_reads << [library, readgroup, "1", file(first_insert)]
	fastqc_reads << [library, readgroup, "2", file(second_insert)]
	fastqc_reads << [library, readgroup, "barcode", file(barcode)]
	rg_number += 1
}

// To expedite the large AMP muscle data processing, do:
// fastqc on one readgroup per library -
// Correct barcode and make a corrected-barcode fastq in the start and count corrected barocdes per library
// Select barcodes with min count of reads
// Select barcodes during trimming and append in read 1 and 2 fastqs there
// Avoid the correct-barcode-in-bam step later on
// Update the barcode tag downstream

workflow {

	make_barcode_corrections_in_chan = Channel.from(make_barcode_corrections_in).groupTuple(sort: true)

	make_barcode_corrections(make_barcode_corrections_in_chan)

	// fastqc was done before, skip for now.
	// // add in randomly selected corrected barcode fastqs into fastqc in channel
	// barcode_fastqc_formatted = make_barcode_corrections.out.barcode_fastqc
	// 	.transpose( by: 1)
	// 	.map{it -> [it[0], it[1].name.replaceAll('_L000_R2_001.corrected.fastq', ''), "barcode", it[1]]}
	// 	.randomSample( 8, 1234321 )
	// fastqc_in = Channel.from(fastqc_reads).concat(barcode_fastqc_formatted)
	// fastqc(fastqc_in)

	// add in the corrected barcode fastq for each library, readgroup and selected list of barcodes for each library into into trim channel
	trim_formatted = make_barcode_corrections.out.barcode_trim
		.transpose( by: 1)
		.map{it -> [it[0], it[1].name.replaceAll('_L000_R2_001.corrected.fastq', ''), it[1], it[2]]}
	trim_in = Channel.from(trim_reads).combine(trim_formatted, by:[0, 1] )
	trim(trim_in)

	// map
	tmp = []
	for (library in libraries) {
		for (genome in get_genome(library)) {
			tmp << [library, genome]
		}
	}
	map_in = Channel.from(tmp).combine(trim.out, by: 0)
	map_bwa(map_in)

	//merge readgroups
	merge_readgroups(map_bwa.out.groupTuple(by: [0, 1], sort: true))

	// mark duplicates
	mark_duplicates(merge_readgroups.out)

	// ataqv
	ataqv(mark_duplicates.out.ataqv_in)
	extract_ataqv(ataqv.out.ataqv_out)
	
}

process make_barcode_corrections {
	
	storeDir "${params.results}/corrected-barcodes"
	tag "${library}"
	cpus 10
	time { 20.h * task.attempt } 
	memory { 55.GB * task.attempt } 
	errorStrategy 'retry'
	container params.container_snatac
	tag "${library}"
	
	input:
	tuple val(library), path(barcode_fastq)

	output:
	path("${library}.barcode_corrections.txt")
	path("${library}.barcode_counts.txt")
	tuple val(library), path("*.corrected.fastq"), emit: barcode_fastqc
	tuple val(library), path("*.corrected.fastq"), path("${library}.selected_barcodes.txt"), emit: barcode_trim

	"""
	${IONICE} correct-barcodes.py --threads 10 ${params.barcode_whitelist} ${barcode_fastq.join(' ')} \
	--barcode-corrections-file ${library}.barcode_corrections.txt --corrected-barcode-counts-file ${library}.barcode_counts.txt ;
	less ${library}.barcode_counts.txt | grep -v "None" | awk '{ if ((\$2>=${params.min_read_pairs})) print \$1}' > ${library}.selected_barcodes.txt
	"""

}


process fastqc {

	publishDir "${params.results}/fastqc", mode: "rellink"
	time '24h'
	memory { 12.GB * task.attempt }
	errorStrategy 'retry'
	container params.container_snatac
	tag "${library}-${readgroup}-${read}"
	
	input:
	tuple val(library), val(readgroup), val(read), path(fastq)

	output:
	tuple path(outfile_1), path(outfile_2)

	script:
	outfile_1 = fastq.getName().replaceAll('.fastq.gz', '_fastqc.html').replaceAll(".corrected_barcodes.fastq", "_fastqc.html")
	outfile_2 = fastq.getName().replaceAll('.fastq.gz', '_fastqc.zip').replaceAll(".corrected_barcodes.fastq", "_fastqc.zip")

	"""
	fastqc $fastq
	"""

}


process trim {

	storeDir "${params.results}/temp-trim"
	errorStrategy 'retry'
	maxRetries 1
	time { 12.h * task.attempt }
	memory { 5.GB * task.attempt }
	tag "${library}-${readgroup}"
	container params.container_snatac
	
	input:
	tuple val(library), val(rgid), val(readgroup), path(fastq_1), path(fastq_2), path(barcode), path(selected_barcodes)

	output:
	tuple val(library), val(readgroup), path("${library}-${readgroup}.1.trimmed.fastq.gz"), path("${library}-${readgroup}.2.trimmed.fastq.gz")

	"""
	${IONICE} cta --append-barcode $barcode --selected-barcodes ${selected_barcodes} \
	$fastq_1 $fastq_2 ${library}-${readgroup}.1.trimmed.fastq.gz ${library}-${readgroup}.2.trimmed.fastq.gz
	"""

}

process map_bwa {
	/* Set the corrected cellular barcode as tag CB */ 

	memory { 50.GB * task.attempt}
	cpus 12
	errorStrategy 'retry'
	maxRetries 1
	time { 24.hour * task.attempt }
	tag "${library}-${readgroup}-${genome}"
	storeDir "${params.results}/temp-bwa"
	container params.container_snatac
	
	input:
	tuple val(library), val(genome), val(readgroup), path(fastq_1), path(fastq_2)

	output:
	tuple val(library), val(genome), path("${library}-${readgroup}-${genome}.bam")

	"""
	bwa mem -I 200,200,5000 -M -t 12 ${get_bwa_index(genome)} ${fastq_1} ${fastq_2} |
	awk -F'\t' '{if ((\$1 !~ /^@/)) { split(\$1,rname,"_"); \$NF=\$NF"\tCB:Z:"rname[2] } print \$0 }' OFS='\t' |
	samtools sort -m 1g -@ 11 -O bam -T sort_tmp -o ${library}-${readgroup}-${genome}.bam ;
	"""

}


process merge_readgroups {

	time '36h'
	storeDir "${params.results}/temp-merge"
	cpus 12
	container params.container_snatac
	tag "${library}-${genome}"
	
	input:
	tuple val(library), val(genome), path(bams)

	output:
	tuple val(library), val(genome), path("${library}-${genome}.bam")

	"""
	samtools merge -@ 11 -O BAM ${library}-${genome}.bam ${bams.join(' ')}
	"""

}


process mark_duplicates {
	/*. Usually one would use CB, but also include as RG because ataqv is hard coded to use that for now */
	publishDir "${params.results}/mark_duplicates", mode: "rellink"
	errorStrategy 'retry'
	maxRetries 1
	time { 24.hour * task.attempt }
	memory { 50.GB * task.attempt }
	container params.container_snatac
	tag "${library}-${genome}"
	
	input:
	tuple val(library), val(genome), path("${library}-${genome}.bam")

	output:
	tuple val(library), val(genome), path("${library}-${genome}.md.bam"), path("${library}-${genome}.md.bam.bai"), emit: prune_in
	tuple val(library), val(genome), path("${library}-${genome}.md.bam"), path("${library}-${genome}.md.bam.bai"), emit: ataqv_in

	"""
	java -Xmx40g -Xms40g -jar \$PICARD_JAR MarkDuplicates TMP_DIR=. I=${library}-${genome}.bam O=${library}-${genome}.mdtemp.bam \
		READ_ONE_BARCODE_TAG=CB READ_TWO_BARCODE_TAG=CB ASSUME_SORTED=true MAX_RECORDS_IN_RAM=100000000 METRICS_FILE=${library}-${genome}.metrics VALIDATION_STRINGENCY=LENIENT ;
	samtools view -h ${library}-${genome}.mdtemp.bam |
	awk -F'\t' '{if ((\$1 !~ /^@/)) { split(\$1,rname,"_"); \$NF=\$NF"\tRG:Z:"rname[2] } print \$0 }' OFS="\t" | 
	samtools view -h - -O BAM > ${library}-${genome}.md.bam ;
	samtools index ${library}-${genome}.md.bam ;
	rm ${library}-${genome}.mdtemp.bam 
	"""

}


process ataqv {
	
	publishDir "${params.results}/ataqv", mode: "rellink"
	errorStrategy 'retry'
	maxRetries 1
	memory { 50.GB * task.attempt }
	time { 24.hour * task.attempt }
	container params.container_snatac
	tag "${library}-${genome}"
	
	input:
	tuple val(library), val(genome), path(md_bam), path(bam_index)
	
	output:
	tuple val(library), val(genome), path("${library}-${genome}.ataqv.json.gz"), emit: ataqv_out
	path("${library}-${genome}.ataqv.out")

	"""
	${IONICE} ataqv --name ${library}-${genome} --metrics-file ${library}-${genome}.ataqv.json.gz --tss-file ${get_tss(genome)} \
		${make_excluded_regions_arg(genome)} ${get_organism(genome)} $md_bam > ${library}-${genome}.ataqv.out
	"""	

}


process extract_ataqv {

	publishDir "${params.results}/ataqv", mode: "rellink"
	memory '60G'
	time '10h'
	container params.container_general
	tag "${library}-${genome}"
	
	input:
	tuple val(library), val(genome), path(ataqv)

	output:
	tuple val(library), path("${library}-${genome}.ataqv.txt")

	"""
	extractAtaqvMetric.py --files $ataqv --metrics tss_enrichment percent_hqaa hqaa total_reads \
	total_autosomal_reads percent_mitochondrial percent_autosomal_duplicate percent_duplicate max_fraction_reads_from_single_autosome > ${library}-${genome}.ataqv.txt
	"""
}

workflow.onComplete {
	if (workflow.success){
		subject = "ATAC execution complete"
	}
	else {
		subject = "ATAC execution error"
	}

	recipient = params.email

	['mail', '-s', subject, recipient].execute() << """

	Pipeline execution summary
	---------------------------
	Completed at: ${workflow.complete}
	Duration	: ${workflow.duration}
	Success		: ${workflow.success}
	workDir		: ${workflow.workDir}
	exit status : ${workflow.exitStatus}
	Error report: ${workflow.errorReport ?: '-'}
	"""
}
