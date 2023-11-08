#!/usr/bin/env nextflow

nextflow.enable.dsl=2

def channel_glob(p) {
	chan = Channel.fromPath(p[0], checkIfExists: true).map{it -> [it.name.replaceAll(p[1], ""), it]}
    return chan
}

def channel_glob_split(p, delim) {
	chan = Channel.fromPath(p[0], checkIfExists: true).map{it -> it.name.replaceAll(p[1], "").tokenize(delim).plus([it])}
    return chan
}


def get_chrom(x, replace){
	return x.name.replaceAll(replace, "").split("\\.")[1]
}

def getfile(f){
	chan = Channel.fromPath(f)
	return chan
}

def abspath(f){
	File file = new File(f)
	path = file.getAbsolutePath()
	return path
}

workflow {
	
	// format QTL vars
	vcf = Channel.fromPath(params.qtl_vcf, checkIfExists: true)
	index = Channel.fromPath(params.qtl_vcf_index, checkIfExists: true)
	format_snps(vcf.combine(index))
	snps = format_snps.out

	//  MR Steiger
	select_mrs(getfile(params.coloc))

	mrs = select_mrs.output
		.flatten()
		.map{it -> it.name.replaceAll(".mrs.tsv","").tokenize('@').plus([it])} // chrom cluster cachunk eqchunk mrs
	mrin = mrs
		.combine(snps)
	// .take(1)
	// mrin.view()
	mrsteiger_eqtl_caqtl(mrin)
}

process select_mrs {
	publishDir "${params.results}/selected", mode: "rellink"
	errorStrategy 'finish'
	time '1h'

	input:
	path(coloc)
	
	output:
	path("*.mrs.tsv")
	
	"""
	mrs_pairs.py --coloc ${coloc} --eqtl-glob ${params.eqtl_glob} --caqtl-glob ${params.caqtl_glob}
	"""
}

process format_snps {
	/* format FUSION QTL SNPs as required by smr */
	storeDir "${params.results}/qtlvars"
	errorStrategy 'retry'
	maxRetries 2
	cpus 10
	time '2h'
	memory '10G'
	tag "qtlvars"
	
	input:
	tuple path(vcf), path(index)
	
	output:
	path("qtl.esi") 

	"""
	vcftools --gzvcf $vcf --freq --out freq
	format-freq.py freq.frq $vcf qtl.esi
	rm freq*
	"""
}

process mrsteiger_eqtl_caqtl {
	/* GWAS v caqtl smr */
	publishDir "${params.results}/eqtl-caqtl-mrsteiger", mode: "rellink"
	time '36h'
	cpus 1
	memory "50G"
	tag "${chrom}.${cluster}.${cachunk}.${eqchunk}"
	container params.container_mrs
	errorStrategy "ignore"
	
	input:
	tuple val(chrom), val(cluster), val(cachunk), val(eqchunk), path(mrs), path(qtlvars)

	output:
	tuple path("${prefix}.*.png"), path("${prefix}.*.tsv") 

	script:
	eqtl_samples = params.eqtl_samples[cluster]
	caqtl_samples = params.caqtl_samples[cluster]
	prefix = "${chrom}.${cluster}.${cachunk}.${eqchunk}"

	"""
	mrsteiger.R $mrs $qtlvars $eqtl_samples $caqtl_samples $prefix
    """
	
}

workflow.onComplete {
	if (workflow.success){
		subject = "SMR execution complete"
	}
	else {
		subject = "SMR execution error"
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

