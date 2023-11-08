#!/usr/bin/env nextflow

nextflow.enable.dsl=2

def channel_glob(p) {
	chan = Channel.fromPath(p[0])
		.map{it -> [it.name.replaceAll(p[1], ""), it]}
		.ifEmpty("!!!\nChannel empty, check ${p}\n!!!\n")
	return chan
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

nperm = params.nperm
workflow{
	// adjust eqtl and caqtl phenotypes for covariates same as in the qtl analyses
	eqpheno = channel_glob(params.eqtl_pheno)
	capheno = channel_glob(params.caqtl_pheno)
	covs = getfile(params.covs)
		
	select_cit(getfile(params.coloc))

	cit_in = select_cit.out
		.flatten()
		.map{it -> it.name.replaceAll(".cit.tsv", "").tokenize("--").plus([it])} //cluster chunk cit
		.combine(eqpheno, by:0)
		.combine(capheno, by:0)
		.combine(covs)
		// .filter{it[0] == "fusion.Type_1"}	
//		.take(1)
	cit(cit_in
//		  .take(3)
	)

	
	fdr_in = cit.out
		.groupTuple(by:0, sort:true)
	fdr_cit(fdr_in)
}


process select_cit {
	publishDir "${params.results}/cit-sh", mode: "rellink"
	errorStrategy 'finish'
	time '1h'

	input:
	path(coloc)
	
	output:
	path("*.cit.tsv")
	
	"""
	cit_pairs.py --coloc ${coloc} --chunk 200
	"""

}

process cit {
	publishDir "${params.results}/cit-perm${nperm}", mode: "rellink"
	errorStrategy 'ignore'
	time '1h'
	tag "${cluster}--${chunk}"
	
	input:
	tuple val(cluster), val(chunk), path(citdf), path("eqpheno"), path("capheno"), path(covs)
	
	output:
	tuple val(cluster), path("${prefix}.citres.rds")

	script:
	prefix = "${cluster}--${chunk}"
	vcf = abspath(params.vcf)
	
	"""
	cit.R --cit $citdf --vcf $vcf --output ${prefix}.citres.rds --eqpheno eqpheno --capheno capheno --covs $covs --nperm ${nperm}
	"""
}

process fdr_cit {
	publishDir "${params.results}/cit_fdr", mode: "rellink"
	errorStrategy 'ignore'
	time '1h'
	tag "${cluster}"
	
	input:
	tuple val(cluster), path(cits)
	
	output:
	path("${cluster}.*.tsv")

	script:
	citdir = abspath("${params.results}/cit-sh")
	
	"""
	cit-fdr.R --cit ${cits.join(',')}  --prefix ${cluster} --citdir $citdir
	"""
}

workflow.onComplete {
	if (workflow.success){
		subject = "Cit execution complete"
	}
	else {
		subject = "Cit execution error"
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

