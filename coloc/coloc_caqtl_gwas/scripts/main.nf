#!/usr/bin/env nextflow

nextflow.enable.dsl=2

def channel_glob(p) {
	chan = Channel.fromPath(p[0]).map{it -> [it.name.replaceAll(p[1], ""), it]}
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

workflow{
	select_coloc()
	rsid = Channel.fromList([[file(params.ukbb_rsid), file(params.ukbb_rsid_index)]])
	// coloc_in = channel_glob(params.coloc_sh)
	coloc_in = select_coloc.out
		.flatten()
		.map{it -> [it.name.replaceAll(".coloc.sh", ""), it]}
		.combine(rsid)
		// .take(1)
	coloc(coloc_in
		  // .take(2)
	)
	plot_in = coloc.out.main
		.filter{! it.name.contains("none.coloc.tsv") }
		.collect()
		.map{it -> file("${params.results}/coloc")}
	plot_coloc(plot_in)
}

process select_coloc {
	publishDir "${params.results}/coloc-sh", mode: "rellink"
	errorStrategy 'finish'
	time '1h'
		
	output:
	path("*.coloc.sh")
	
	"""
	coloc_pairs.py --caqtl ${params.caqtl} --gwas ${params.gwas.join(' ')}
 	"""

}

process coloc {
	publishDir "${params.results}/coloc", mode: "rellink"
	errorStrategy 'ignore'
	time '1h'
	container params.container_coloc
	tag "$region"
	
	input:
	tuple val(region), path(sh), path(rs_vcf), path(index)
	
	output:
	path("*.coloc.tsv"), emit: main
	path("${region}.log")
	path("*coloc-sensitivity*")

	script:
	regionslice = region.split("\\.")[-1]
	
	"""
	tabix $rs_vcf $regionslice > ${region}.rsid
	bash $sh
	rm ${region}.rsid
	n=`ls *.coloc.tsv | wc -l`;
	if [[ \$n == 0 ]]; then touch ${region}.none.coloc.tsv; fi
	ln -s .command.log ${region}.log
	"""

}

process plot_coloc {
	publishDir "figures", mode: "rellink"
	errorStrategy 'finish'
	time '1h'

	input:
	path(coloc_dir)
	
	output:
	path("*.pdf")
	path("*.png")
	path("summary.tsv")
	
	"""
	plot_coloc.py --coloc ${coloc_dir} --pheno ${abspath(params.pheno)}
 	"""

}

workflow.onComplete {
	if (workflow.success){
		subject = "Coloc execution complete"
	}
	else {
		subject = "Coloc execution error"
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

