#!/usr/bin/env nextflow
nextflow.enable.dsl=2

def channel_glob(p) {
	chan = Channel.fromPath(p[0], type: 'any').map{it -> [it.name.replaceAll(p[1], "").replaceAll("atac-", ""), it]}
	return chan
}

def channel_glob_split(p, delim) {
	chan = Channel.fromPath(p[0], type: 'any')
		.map{it -> it.name.replaceAll(p[1], "").replaceAll("atac-", "").tokenize(delim).plus([it])}
	return chan
}

def getfile(p){
	return Channel.fromPath(p)
}


workflow {

	// sample-level covariates
	sample_ch = Channel.fromPath(params.sample_info)

	/*Take samples selected for any cluster, subset the vcf to prep for genotype PCA, etc. */
	sample_ids = Channel.fromPath(params.samplelists)
		.collect()
		.map{it -> [it]}
	vcf = Channel.fromList([[file(params.vcf), file(params.vcf_index)]])
	prep_in = sample_ids.combine(vcf)
	prep_vcf(prep_in)

	pca_genotype(prep_vcf.out.main)
	plot_pca_genotype(pca_genotype.out.main
					  .combine(sample_ch)
	)

	atac_counts = channel_glob(params.atac_mtx_dir) // cluster mtx
	
	sample_lists = channel_glob(params.selected_samples) // cluster samplelist
	
	summit_lists = channel_glob(params.selected_summits) // cluster features
	
	normalize_in = atac_counts
		.combine(sample_lists, by: 0)
		.combine(summit_lists, by: 0)

	normalize_tpm(normalize_in)

	pca_phenotype(normalize_tpm.out)

	covs_in = pca_phenotype.out.main
	.combine(pca_genotype.out.main)
	.combine(sample_ch)
	make_covariates(covs_in)
	make_cov_lists()

	qtl = normalize_tpm.out 
		.combine(make_covariates.out, by: 0)
		.combine(prep_vcf.out.main) // cluster, bed, index, cov_file, vcf, vcf_index
	covlists  = make_cov_lists.out.flatten().map{it -> [it.name.replaceAll(".txt",""), it]} // covtype covlist

	chunks = 1..params.total_chunks
	c_permute = Channel.fromList(chunks)

	permute_in  = qtl
		.combine(covlists) 	// val(cluster), path(bed), path(index), path(cov_file), path(vcf), path(vcf_index), val(covtype), path(covlist)
		.combine(c_permute) // val(cluster), path(bed), path(index), path(cov_file), path(vcf), path(vcf_index), val(covtype), path(covlist), val(chunk)
		// .filter{params.select_clusters.contains(it[0])}
	permute(permute_in)

}

process prep_vcf {
	
	storeDir "${params.results}/vcf"
	memory '2GB'
	time '5h'
	
	input:
	tuple path(samplelists), path(vcf), path(index)
	
	output:
	tuple path("filtered-vcf.maf0.05-hwe1e6.recode.vcf.gz"), path("filtered-vcf.maf0.05-hwe1e6.recode.vcf.gz.tbi"), emit: main
	path("keep-samples.txt")
	
	"""
	cat ${samplelists.join(' ')} | sort | uniq | grep -v SNG.1ST > keep-samples.txt ; 
	vcftools --gzvcf $vcf --maf 0.05 --min-alleles 2 --max-alleles 2 --keep keep-samples.txt --hwe 0.000001 --recode --recode-INFO-all --out filtered-vcf.maf0.05-hwe1e6 --not-chr chrX;
	bgzip filtered-vcf.maf0.05-hwe1e6.recode.vcf ;
	tabix filtered-vcf.maf0.05-hwe1e6.recode.vcf.gz
	"""
}

process pca_genotype {
	storeDir "${params.results}/pca-genotype"
	errorStrategy 'retry'
	memory '1GB'
	time '5h'
	container params.container_qtl
	
	input:
	tuple path(vcf), path(index)
	
	output:
	path("scaled-centered-genotype.pca"), emit: main
	path("scaled-centered-genotype.pca_stats")

	"""
	QTLtools pca --vcf $vcf --scale --center --out scaled-centered-genotype	 --maf 0.05 --distance 50000
	"""
}

process plot_pca_genotype {
	storeDir "${params.results}/pca-genotype-plot"
	errorStrategy 'retry'
	memory '1GB'
	time '5h'
	
	input:
	tuple path(pc), path(cluster_info)
	
	output:
	path("fig.genotype-pcs.png")
	
	"""
	plot-pcs.py --pc $pc --cluster-info ${cluster_info} --output fig.genotype-pcs.png
	"""
}

process normalize_tpm {
	/* Subest samples and summit features. Normalize to TPMs. bed.gz should have TSS coordinates */
	
	storeDir "${params.results}/tpm-matrices/"
	memory '20GB'
	time '5h'
	tag "${cluster}"
	
	input:
	tuple val(cluster), path(counts), path(samplelist), path(summitlist)

	output:
	tuple val(cluster), path("${cluster}.bed.gz"), path("${cluster}.bed.gz.tbi")

	"""
	normalize.R --peak_counts $counts --sample_list ${samplelist} --summit_list ${summitlist}  --output_bed ${cluster}.bed ; 
	less ${cluster}.bed | sed -e '1 s:chrom:#Chr:g'	 > temp; mv temp ${cluster}.bed; 
	bgzip ${cluster}.bed; tabix ${cluster}.bed.gz
	"""

}

process pca_phenotype {
	storeDir "${params.results}/pca-phenotype"
	memory '40GB'
	time '5h'
	tag "${cluster}"
	container params.container_qtl
	
	input:
	tuple val(cluster), path(bed), path(index)

	output:
	tuple val(cluster), path("${cluster}.scaled-centered-phenotype.pca"), emit: main
	path("${cluster}.scaled-centered-phenotype.pca_stats")
	
	"""
	QTLtools pca --bed $bed --scale --center --out ${cluster}.scaled-centered-phenotype 
	"""
	
}


process make_covariates {
	/*Organize covariates info such as PCs, sample info etc. */
	storeDir "${params.results}/sample-info"
	memory '4GB'
	time '5h'
	
	input:
	tuple val(cluster), path(pheno_pcs), path(geno_pcs), path(sample_covs)
	
	output:
	tuple val(cluster), path("${cluster}-covariates.tsv")

	"""
	assemble-pcs.py --pheno ${pheno_pcs} --geno ${geno_pcs} --sample-info ${sample_covs} --output ${cluster}-covariates.tsv  --cluster $cluster
	"""
}


process make_cov_lists {
	/* Make lists of covariates to include in each scan */
	publishDir "${params.results}/covariate-lists", mode: "rellink"
	memory '1GB'
	time '5h'
		
	output:
	path("covlist*.txt")

	"""
	make-covariate-lists.py --pheno-pcs ${params.n_pheno_pcs.join(' ')} --geno-pcs ${params.n_geno_pcs.join(' ')} --other-vars ${params.other_covs.join(' ')}
	"""
}


process permute {
	publishDir "${params.results}/permute-window${params.window}/", mode: 'rellink', overwrite: true
	memory '1GB'
	errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	time '8h'
	tag "${cluster}--${covtype}"
	container params.container_qtl
	
	input:
	tuple val(cluster), path(bed), path(index), path(cov_file), path(vcf), path(vcf_index), val(covtype), path(covlist), val(chunk)

	output:
	tuple val("${cluster}--${covtype}"), path("${prefix}.txt"), emit: main
	tuple val(cluster), val(covtype), path("${prefix}.txt"), emit: gettested

	script:
	prefix = "${cluster}--${covtype}--${chunk}"
	
	"""
	QTLtools cis --bed $bed --vcf ${vcf} --cov ${cov_file}  --out ${prefix}.txt --permute ${params.permutations} --chunk $chunk $params.total_chunks \
		--window ${params.window} --include-covariates ${covlist} --normal --std-err
	"""

}


workflow.onComplete {
	if (workflow.success){
		subject = "QTLtools execution complete"
	}
	else {
		subject = "QTLtools execution error"
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

