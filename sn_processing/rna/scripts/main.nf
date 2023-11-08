#!/usr/bin/env nextflow

nextflow.enable.dsl=2
IONICE = 'ionice -c2 -n7'

libraries = params.libraries.keySet()

def get_star_index(genome) {
	return params.star_index[genome]
}

def get_gtf(genome) {
	return(get_star_index(genome) + '/annotation.gtf')
}

def get_chrom_sizes(genome) {
	return(get_star_index(genome) + '/chrNameLength.txt')
}
	
def get_genome(library) {
	return params.libraries[library].genome
}

def library_to_readgroups(library) {
	return params.libraries[library].readgroups.keySet()
}

def library_and_readgroup_to_fastqs(library, readgroup) {
	return params.libraries[library].readgroups[readgroup]
}


fastq_in = []
fastqc_in = []

for (library in libraries) {
	for (readgroup in library_to_readgroups(library)) {
		fastqs = library_and_readgroup_to_fastqs(library, readgroup)
		insert_read = fastqs['2']
		barcode_read = fastqs['1']
		fastqc_in << [library, readgroup, file(insert_read)]
		fastqc_in << [library, readgroup, file(barcode_read)]
		for (genome in get_genome(library)) {
			fastq_in << [library, genome, file(barcode_read), file(insert_read)]
		}
	}
}

star_in = Channel.from(fastq_in).groupTuple(by: [0,1 ])

workflow {
	// run starsolos
	starsolo(star_in)

	// get counts
	feature_counts(starsolo.out.feature_counts_in)

	// compile qc
	qc(feature_counts.out.qc_in)

	// prune bams
	prune(starsolo.out.prune_in)

	// extract rna qc
	rna_qc(qc.out)

	// plot qc
	plot_qc(rna_qc.out.collect())
}

process starsolo {
	/* For this, trying out running on /localscratch reduced runtime of a testjob by 50%
	 However, there is limited space on /localscratch.
	 Even when the running process deletes the intermediate files such as for bam sort,
	 because the process still is running, the filesystem shows that disk space is still occupied by
	 those running jobs. The resulting bam files are complained to be corrupted. 
	 To circumvent this, increase cores so that limited instances of this process run at once or include maxForks 
	 New in hg38 analyses: removed --sjdbGTFfile ${get_gtf(genome)} flag since we are already using the index.. */
	publishDir "${params.results}/starsolo/${library}-${genome}", mode: "rellink"
	memory { 120.GB * task.attempt }
	cpus 20
	time '48h'
	errorStrategy 'retry'
	maxRetries 1
	container params.container_snrna
	tag "${library}-${genome}"
	
	input:
	tuple val(library), val(genome), path(barcode_fastq), path(insert_fastq)

	output:
	tuple val(library), path("Aligned.sortedByCoord.out.bam"), path("Log.final.out"), path("Log.out"), path("Log.progress.out"), path("SJ.out.tab"), path("Solo.out"), emit: starsolo_out
	tuple val(library), val(genome), path("Aligned.sortedByCoord.out.bam"), emit: feature_counts_in
	tuple val(library), val(genome), path("Aligned.sortedByCoord.out.bam"), emit: prune_in

	script:
	soloUMIlen = params.chemistry == 'V2' ? 10 : 12

	"""
	${IONICE} STAR --soloBarcodeReadLength 0 --runThreadN 20 --genomeLoad NoSharedMemory --runRNGseed 789727 --readFilesCommand gunzip -c \
		--outSAMattributes NH HI nM AS CR CY CB UR UY UB sM GX GN --genomeDir ${get_star_index(genome)} --outSAMtype BAM SortedByCoordinate \
		--outSAMunmapped Within KeepPairs  --soloType Droplet --soloUMIlen $soloUMIlen --soloFeatures Transcript3p Gene GeneFull SJ Velocyto \
		--soloUMIfiltering MultiGeneUMI --soloCBmatchWLtype 1MM_multi_pseudocounts --soloCellFilter None --soloCBwhitelist ${params.barcode_whitelist} \
		--readFilesIn ${insert_fastq.join(',')} ${barcode_fastq.join(',')}
	"""

}


process feature_counts {

	cpus 20
	storeDir "${params.results}/temp-feature-counts"
	container params.container_snrna
	tag "${library}-${genome}"
	
	input:
	tuple val(library), val(genome), path("star.bam")

	output:
	tuple val(library), val(genome), path("${library}-${genome}.featureCounts.bam"), emit: qc_in
	

	"""
	featureCounts -a ${params.gene_gtf} -T 20 -t transcript -o ${library}-${genome} star.bam -R BAM -O -s 1; mv star.bam.featureCounts.bam ${library}-${genome}.featureCounts.bam
	"""

}


process qc {
	/* use cython import of cell class	*/
	publishDir "${params.results}/qc", mode: "rellink"
	memory { 130.GB * task.attempt }
	time '36h'
	maxRetries 1
	errorStrategy 'retry'
	container params.container_snrna
	tag "${library}-${genome}"
	
	input:
	tuple val(library), val(genome), path(bam)

	output:
	tuple val(library), val(genome), path("${library}-${genome}.qc.json.gz")

	"""
	qc-from-featurecounts.py --cell-tag CB --gene-tag XT --umi-tag UB --min-reads 100 $bam | gzip -c > ${library}-${genome}.qc.json.gz
	"""

}

process prune {

	publishDir "${params.results}/prune", mode: 'rellink', overwrite: true
	maxForks 10
	time '14h'
	container params.container_snrna
	tag "${library}-${genome}"
	
	input:
	tuple val(library), val(genome), path(bam)

	output:
	tuple path("${library}-${genome}.before-dedup.bam"), path("${library}-${genome}.before-dedup.bam.bai")

	"""
	${IONICE} samtools view -h -b -q 255 -F 4 -F 256 -F 2048 $bam > ${library}-${genome}.before-dedup.bam
	samtools index ${library}-${genome}.before-dedup.bam
	"""

}


process rna_qc {
	/* Convert json to dataframe to plot */
	
	publishDir "${params.results}/qc"
	memory '20 GB'
	container params.container_snrna
	tag "${library}-${genome}"
	
	input:
	tuple val(library), val(genome), path(jsn)

	output:
	path("${library}-${genome}.txt")

	"""
	extract-rnaseq-qc.py $jsn ${library}-${genome}.txt
	"""
}


process plot_qc {
	/* Plot RNA qc */

	storeDir "figures"
	container params.container_general
	
	input:
	path(qc)

	output:
	tuple path("fig.umis_joint.png"), path("fig.umis_cumulative_fraction.pdf"), path("fig.n_nuclei_umis.pdf"), path("fig.umis_point_density.png"), path("libs.tsv") 

	"""
	plot-qc.py --input ${qc.join(' ')} --data libs.tsv --facet-wrap 4;
	plot-qc.R libs.tsv fig.umis_point_density.png
	"""
}


workflow.onComplete {
	if (workflow.success){
		subject = "RNA execution complete"
	}
	else {
		subject = "RNA execution error"
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

