#!/usr/bin/env nextflow

/* Run SuSiE on GWAS for each signal */

def channel_glob(p) {
	chan = Channel.fromPath(p[0]).map{it -> [it.name.replaceAll(p[1], ""), it]}
    return chan
}

nextflow.enable.dsl=2

workflow{

	// Select GWAS regions to run SuSiE on and prep sh files
	select_in = Channel.fromPath(params.config)
	select_susie(select_in)

	// run SuSiE prep
	susie_prep_in = select_susie.out.susie_prep
		.flatten()
		.map{it -> [it.name.replaceAll(".susieprep.sh", ""), it]} // locus, prep script
		.filter{! it[0].contains("chrX")}
		// .take(1)
	susie_prep(susie_prep_in)

	// Make LD matrix
	ld_in = susie_prep.out.forld // locus reference gwas
		// .filter{! it[1].contains("diamante_T2D-European__MHC:") }
	make_ld(ld_in)

	// Run SuSiE
	susie_sh = select_susie.out.susie_run
		.flatten()
		.map{it -> [it.name.replaceAll(".susie.sh", ""), it]} // locus, run_script
	susie_in = make_ld.out
		.combine(susie_prep.out.forsusie, by: 0) // locus ldmat gwas
		.combine(susie_sh, by: 0) // locus ldmat gwas run_script
	susie(susie_in
		  // .take(1)
	)
}

process select_susie {
	/* Starting with GWAS leads, define regions to run SuSiE on. Requires the config.yaml file */
	publishDir "${params.results}/susie-region", mode: "rellink"
	errorStrategy 'finish'
	time '1h'
	container params.container_susie
	
	input:
	path(config)
	
	output:
	path("*.susieprep.sh"), emit: susie_prep
	path("*.susie.sh"), emit: susie_run

	"""
	susie-regions.py --config $config
	"""
}

process susie_prep {
	/* Prep dosage files for each susie region selected */
	storeDir "${params.results}/susie-prepped"
	errorStrategy 'finish'
	time '4h'
	tag "${locus}"
	container params.container_susie

	input:
	tuple val(locus), path(prep)

	output:
	tuple val(locus), path("${locus}.ukbb-dosages.tsv"), path("${locus}.gwas.tsv"), emit: forld
	tuple val(locus), path("${locus}.gwas.tsv"), emit: forsusie
	
	
	"""
	bash $prep
	"""

}

process make_ld {
	/* make ld matrix */
	storeDir "${params.results}/ld-mat"
	errorStrategy 'finish'
	time '4h'
	memory "50G"
	tag "${locus}"
	container params.container_susie
	
	input:
	tuple val(locus), path(ukbb_ref), path(gwas)

	output:
	tuple val(locus), path("${locus}.ld.tsv")
	
	"""
	make-ld-mat.R --trait1 ${gwas} --prefix ${locus} --trait1_ld ${ukbb_ref}
	"""
}

process susie {
	/* run susie */
	publishDir "${params.results}/susie", mode: "rellink"
	publishDir "${params.results}/for-coloc", mode: "rellink", pattern: "*selected.Rda"
	publishDir "${params.results}/for-coloc", mode: "rellink", pattern: "*selected.png"
	
	errorStrategy 'ignore'
	time '10h'
	memory "50G"
	tag "${locus}"
	container params.container_susie
	
	input:
	tuple val(locus), path(ldmat), path(gwas), path(susie_run)

	output:
	tuple val(locus), path("${locus}.results.tsv"), emit: results
	tuple val(locus), path("${locus}*.Rda"), emit: rdata
	tuple val(locus), path("${locus}*.png"), emit: png
	path("${locus}.log")
	
	"""
	bash $susie_run
	ln -s .command.log ${locus}.log
	"""
}


workflow.onComplete {
	if (workflow.success){
		subject = "SuSiE execution complete"
	}
	else {
		subject = "SuSiE execution error"
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

