#!/usr/bin/env nextflow

nextflow.enable.dsl=2

def channel_glob(p) {
	chan = Channel.fromPath(p[0]).map{it -> [it.name.replaceAll(p[1], ""), it]}
    return chan
}

def abspath(f){
	File file = new File(f)
	path = file.getAbsolutePath()
	return path
}

susie_dir = abspath("${params.results}/susie_window${params.window}_min_abs_cor${params.min_abs_cor}")

workflow{
	loci = channel_glob(params.loci)
		.map{it -> it[0].split("__")[0].tokenize("--").plus([it[1]])}
	beds = channel_glob(params.bed).map{it -> it + ["${it[1]}.tbi"]}
	covs = channel_glob(params.covs)
	
	start = loci
		.combine(beds, by:0)
		.combine(covs, by:0)
		.map{it -> it.plus([file("${params.covlists}/${it[1]}.txt", checkIfExists: true)])}
		// .take(1)
	adjust_phenotypes(start)
	vcf = Channel.fromList([[file(params.vcf), file(params.vcf_index)]])
	dosages(vcf)
	susie_in = adjust_phenotypes.out
		.transpose().map{it -> it.plus(it[2].name)}
		.combine(dosages.out)
		// .take(1)
	susie(susie_in)
}

process adjust_phenotypes {
	publishDir "${params.results}/adjust"
	errorStrategy 'finish'
	time '1h'
	tag "${cluster}--${covtype}"
	container params.container_qtl
	
	input:
	tuple val(cluster), val(covtype), path(eqtl), path(bed), path(index), path(cov), path(covlist) 
	
	output:
	tuple val(prefix), path("${prefix}.adjusted.bed"), path("${prefix}.pheno-list.*") 

	script:
	prefix = "${cluster}--${covtype}"
	
	"""
	QTLtools correct --bed $bed --cov $cov --normal --include-covariates $covlist --out ${prefix}.adjusted.bed  --include-phenotypes $eqtl;
	split -l 500 $eqtl ${prefix}.pheno-list.
	"""

}

process dosages {
	storeDir "${params.results}/genotype_dosages"
	errorStrategy 'finish'
	time '1h'
	
	input:
	tuple path(vcf), path(vcf_index)
	
	output:
	tuple path("genotype_dosages.tsv.gz"), path("genotype_dosages.tsv.gz.tbi")
	
	"""
	less $vcf | head -150 | grep "#CHROM" | sed -e 's:#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\tFORMAT:CHROM\\tPOS\\tREF\\tALT:g' > header.txt
	bcftools query -f "%CHROM\\t%POS\\t%REF\\t%ALT[\\t%DS]\\n" ${vcf} | cat header.txt - | bgzip > genotype_dosages.tsv.gz
    tabix -s1 -b2 -e2 -S1 genotype_dosages.tsv.gz
	"""

}

process susie {
	publishDir "${params.results}/susie_window${params.window}_min_abs_cor${params.min_abs_cor}/${prefix}"
	errorStrategy 'finish'
	time '1h'
	memory "10G"
	tag "${prefix}"
	container params.container_susie
	
	input:
	tuple val(prefix), path(adjusted), path(gene_names), val(prefix), path(dosages), path(tabix)

	output:
	tuple path("${prefix}.*susie.Rda"),  path("${prefix}.*.tsv"), path("*.png"), emit: main
	path("${prefix}.log"), emit: log

	
	"""
	susie.R --phenotype $adjusted --dosage $dosages --gene_name $gene_names --prefix $prefix --window ${params.window} --min_abs_cor ${params.min_abs_cor} &> ${prefix}.log
	"""

}


workflow.onComplete {
	if (workflow.success){
		subject = "Susie execution complete"
	}
	else {
		subject = "Susie execution error"
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

