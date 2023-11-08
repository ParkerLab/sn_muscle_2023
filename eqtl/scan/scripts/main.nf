#!/usr/bin/env nextflow

/*
Prep to run eQTL from snRNA - get counts matrix for RNA with sample counts or gene TPMs
Plot QC metrics such as PCA etc. 
Select genes and samples and covariates to include
*/

nextflow.enable.dsl=2

def channel_glob(p) {
	chan = Channel.fromPath(p[0], type: 'any').map{it -> [it.name.replaceAll(p[1], ""), it]}
	return chan
}

def getfile(p){
	return Channel.fromPath(p)
}

def abspath(f){
	File file = new File(f)
	path = file.getAbsolutePath()
	return path
}

rna_counts_dir = abspath("${params.results}/input-features")
permute_formatted_dir = abspath("${params.results}/permute-formatted-window${params.window}")
covariate_dir = abspath("${params.results}/sample_info")
covariate_list_dir = abspath("${params.results}/covariate-lists")
sample_id_dir = abspath("${params.results}/sample_id_lists")

workflow{

	cluster_info = getfile(params.cluster_info)
	gene_info = getfile(params.gene_info)
	rna_counts = channel_glob(params.rna_mtx_dir)

	/*Take samples selected for any cluster, subset the vcf to prep for genotype PCA, etc. */
	select_samples(cluster_info)
	sample_ids = select_samples.out.id
		.map{it -> [it]}
	vcf = Channel.fromList([[file(params.vcf), file(params.vcf_index)]])

	prep_in = sample_ids.combine(vcf)
	prep_vcf(prep_in)

	sample_ch = getfile(params.sample_info)
	pca_genotype(prep_vcf.out.main)
	plot_pca_genotype(pca_genotype.out.main
					  .combine(sample_ch)
	)

	/* Format RNA counts  */
	cin = rna_counts
		.combine(cluster_info)
		.combine(gene_info)
	counts_matrix(cin)

	norm_in = counts_matrix.out
		.transpose()
		.map{it -> [it[1].name.replaceAll(".tsv",""), it[1]]}
	normalize_tpm(norm_in)
	
	pca_phenotype(normalize_tpm.out.tpm_out)
	covs_in = pca_phenotype.out.main
		.combine(pca_genotype.out.main)
		.combine(sample_ch)
	make_covariates(covs_in)

	// make lists of covariates to test
	make_cov_lists()
	covlists  = make_cov_lists.out.flatten().map{it -> [it.name.replaceAll(".txt",""), it]} // covtype covlist

	chunks = 1..params.total_chunks
	c_permute = Channel.fromList(chunks)

	qtl_in = normalize_tpm.out.cluster_tpms
		.combine(make_covariates.out, by: [0])
		// .filter{selected.contains(it[0])}
		.combine(prep_vcf.out.main)
		.combine(covlists) 
	.combine(c_permute) // val(cluster), path(bed), path(index), path(cov_file), path(vcf), path(vcf_index), val(covtype), path(covlist), val(chunk)
	
	permute(qtl_in)
	compile_chunks(permute.out.groupTuple(by: 0))
	scan_results(compile_chunks.out.main.collect())

	// plot results, get optimized PCs. link those permute outputs. titrate and select testable features
	plot_scan_results(scan_results.out)

	// run cis scan on final eQTL
	selected = plot_scan_results.out.selected_features
		.flatten()
		.map{it -> it.name.replaceAll(".tested-features", "").tokenize("--").plus([it])}

	// run full cis scan and index outputs
	cis_in = selected // cluster covtype featurelist
		.combine(normalize_tpm.out.cluster_tpms, by: 0) // cluster covtype featurelist bed index
		.combine(make_covariates.out, by: 0) // cluster covtype featurelist bed index cov_file
		.combine(prep_vcf.out.main) // cluster covtype featurelist bed index cov_file vcf vcf_index
		.combine(chunks)
	cis(cis_in
		// .take(2)
	)
	cis_index(cis.out)
}


process select_samples {
	/* Make selected sample lists - select samples with min 10 RNA nuclei for all clusters*/
	publishDir "${params.results}/sample-lists", pattern: "*.sample_list.txt"
	publishDir "${sample_id_dir}", pattern: "*.sample_id_list.txt"
	publishDir "${params.results}", pattern: "*.tsv"
	errorStrategy 'retry'
	memory '8GB'
	time '5h'
	
	input:
	path(info)

	output:
	path("fusion.*.sample_list.txt"), emit: sample
	path("fusion.*.sample_id_list.txt"), emit: id
	path("fusion.cluster_selected_sample_counts.tsv")

	"""
	select-samples.py --cohort FUSION --modality rna --cluster-info ${info} --min-nuclei 10 --prefix fusion
	"""
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


process counts_matrix {
	/* Using nuclei RNA counts matrix, get sum of gene counts per cluster per sample. Select out FUSION samples*/
	
	publishDir "${rna_counts_dir}"
	memory '30G'
	time '5h'
	tag "${dec}"
	
	input:
	tuple val(dec), path(countsdir), path(cinfo), path(ginfo)

	output:
	tuple val(cohort), path("${cohort}*.tsv")

	script:
	cohort = "fusion"
	
	"""
	rna-counts.py --counts-dir $countsdir --cluster-info ${cinfo}  --prefix ${cohort} --aggregate \
		--gene-info ${ginfo} --min-fusion-samples 10 --gene-types ${params.gene_types.join(' ')} --exclude-chr ${params.exclude_chrs.join(' ')}

	"""
}

process normalize_tpm {
	/* Normalize to TPMs. bed.gz should have TSS coordinates */

	publishDir "${params.results}/tpm-matrices/"
	memory '2GB'
	time '5h'
	tag "${cluster}"

	input:
	tuple val(cluster), path(counts)

	output:
	tuple val(cluster), path("${cluster}.bed.gz"), path("${cluster}.bed.gz.tbi"), path("${cluster}.top10kgenes.txt"), emit: tpm_out
	tuple val(cluster), path("${cluster}.bed.gz"), path("${cluster}.bed.gz.tbi"), emit: cluster_tpms
	path("${cluster}.median_tpm.txt")
	
	"""
	normalize.R $counts ${cluster} ; 
	less ${cluster}.bed | sed -e '1 s:X.Chr:#Chr:g' | sed -e '1 s:\\tX:\\t:g' > temp; mv temp ${cluster}.bed; 
	bgzip ${cluster}.bed; tabix ${cluster}.bed.gz
	"""

}

process pca_phenotype {

	publishDir "${params.results}/pca-phenotype/"
	memory '1GB'
	time '5h'
	tag "${cluster}"
	container params.container_qtl
	
	input:
	tuple val(cluster), path(bed), path(index), path(top10kphenotypes)

	output:
	tuple val(cluster), path("${cluster}.scaled-centered-phenotype.pca"), emit: main
	path("${cluster}.scaled-centered-phenotype.pca_stats")
	
	"""
	QTLtools pca --bed $bed --scale --center --out ${cluster}.scaled-centered-phenotype --include-phenotypes $top10kphenotypes
	"""
	
}

process make_covariates {
	publishDir "${covariate_dir}"
	memory '2GB'
	time '5h'
	tag "${cluster}"
	
	input:
	tuple val(cluster), path(pheno_pcs), path(geno_pcs), path(sample_covs)

	output:
	tuple val(cluster), path("${cluster}.covariates.tsv")

	"""
	assemble-pcs.py --pheno ${pheno_pcs} --geno ${geno_pcs} --sample-info ${sample_covs} --output ${cluster}.covariates.tsv	 --cluster ${cluster}
	"""
}

process make_cov_lists {
	/* Make lists of covariates to include in each scan */
	publishDir "${covariate_list_dir}", mode: "rellink"
	memory '1GB'
	time '5h'
	
	output:
	path("covlist*.txt")
	
	"""
	make-covariate-lists.py --pheno-pcs ${params.n_pheno_pcs.join(' ')} --geno-pcs ${params.n_geno_pcs.join(' ')} --other-vars ${params.other_covs.join(' ')}
	"""
}


process permute {
	publishDir "${params.results}/permute/", mode: 'rellink', overwrite: true
	memory '12GB'
	errorStrategy 'finish'
	time '8h'
	tag "${cluster}--${covtype}--${chunk}"
	container params.container_qtl

	input:
	tuple val(cluster), path(bed), path(index), path(cov_file), path(vcf), path(vcf_index), val(covtype), path(cov_list), val(chunk)

	output:
	tuple val("${cluster}--${covtype}"), path("${prefix}.txt")

	script:
	prefix = "${cluster}--${covtype}--${chunk}"
	
	"""
	QTLtools cis --bed $bed --vcf ${vcf} --cov ${cov_file}	--out ${prefix}.txt --permute ${params.permutations} --chunk $chunk $params.total_chunks \
		--window ${params.window} --include-covariates ${cov_list} --normal --std-err
	"""

}

process compile_chunks {
	publishDir "${permute_formatted_dir}", mode: "rellink"
	memory '8GB'
	time '5h'
	tag "${prefix}"
	
	input:
	tuple val(prefix), path(permute_outs)
	
	output:
	tuple val(prefix), path("${prefix}.permute.tsv"), emit: permute_compile_out
	path("${prefix}.permute.tsv"), emit: main
	
	"""
	qtl.py --qtl ${permute_outs.join(' ')} --output ${prefix}.permute.tsv ;
	qvalue.R --df  ${prefix}.permute.tsv --p p_beta --q qvalue_storey --output temp ; mv temp  ${prefix}.permute.tsv ;
	"""	
}

process scan_results {
	publishDir "${params.results}/scan_results", mode: "rellink"
	memory '8GB'
	time '5h'

	input:
	path(permute_outs)

	output:
	path("permute-scan*")
	
	"""
	qtl-stats.py --qtl ${permute_outs.join(' ')} --output permute-scan --threshold 0.05 --qvalue qvalue_storey ;
	"""
}

process plot_scan_results {
	/* identify optimal number of PCs. Then, identify n testable features that maximize eQTL after multiple testing correction.*/
	publishDir "${params.results}/scan_results", mode: "rellink"
	publishDir "${params.results}/permute-freeze-window${params.window}", mode: "rellink", pattern: "*.permute.tsv"
	publishDir "${params.results}/fetaure-lists", mode: "rellink", pattern: "*.tested"
	publishDir "${params.results}/fetaure-lists", mode: "rellink", pattern: "*.significant"
	publishDir "${params.results}/selected-fetaure-lists", mode: "rellink", pattern: "*-features"
	memory '8GB'
	time '5h'
	
	input:
	path(scan_result)
	
	output:
	path("permute-scan.max_signals.tsv"), emit: main
	path("*.permute.tsv"), emit: max_sig
	path("*.png")
	path("*.pdf")
	path("*.tsv")
	path("*.tested")
	path("*.significant")
	path("*.tested-features"), emit: selected_features
	path("*.significant-features")
	
	"""
	plot-pc-scan.py --n-pheno-pcs ${params.n_pheno_pcs.join(' ')} --n-geno-pcs ${params.n_geno_pcs.join(' ')} \
		--other-covs ${params.other_covs.join(' ')} --scan-result ${scan_result} --prefix permute-scan ;

	for i in `cat permute-scan.max_signals.tsv | grep -v name | cut -f1`; do p=`readlink -f ${permute_formatted_dir}/\${i}.permute.tsv`; cp \$p .; done
	titrate.py --max-signals permute-scan.max_signals.tsv --permute-dir ${permute_formatted_dir} --rna-counts-dir ${rna_counts_dir} --sample-id-dir ${sample_id_dir} 
	"""
}


process cis {
	publishDir "${params.results}/temp_cis/", mode: 'rellink', overwrite: true
	memory '1GB'
	errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	time '8h'
	tag "${cluster}--${covtype}"
	container params.container_qtl
	
	input:
	tuple val(cluster), val(covtype), path(featurelist), path(bed), path(index), path(cov_file), path(vcf), path(vcf_index), val(chunk)
	
	output:
	tuple val(prefix), path("${prefix}.txt")

	script:
	prefix = "${cluster}--${covtype}--${chunk}"
	
	"""
	QTLtools cis --bed $bed --vcf ${vcf} --cov ${cov_file}	--out ${prefix}.txt --nominal 1 --chunk $chunk $params.total_chunks \
		--window ${params.cis_window} --include-covariates ${covariate_list_dir}/${covtype}.txt --normal --std-err --include-phenotypes $featurelist
	"""

}

process cis_index {
	publishDir "${params.results}/cis-indexed/", mode: 'rellink', overwrite: true
	memory '1GB'
	errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	time '8h'
	tag "${prefix}"
	
	input:
	tuple val(prefix), path(res)

	output:
	tuple path("${prefix}.bed.gz"), path("${prefix}.bed.gz.tbi")
	
	"""
	less $res | awk -F' ' '{if ((\$NF!="NA")) print \$9,\$11-1,\$11,\$8,\$1,\$2,\$3,\$4,\$5,\$6,\$7,\$12,\$13,\$14,\$15,\$16}' OFS='\\t' | sort -k1,1 -k2,2n -k3,3n > temp; 
	echo -e "#snp_chrom\\tsnp_start\\tsnp_end\\tsnp\\tgene_name\\tchrom\\tstart_pheno\\tend_pheno\\tstrand\\tn_variants_tested\\tdistance_var_pheno\\tp_nominal\\tr2\\tslope\\tse\\tbest_hit" | 
		cat - temp > ${prefix}.bed
	bgzip ${prefix}.bed
	tabix --zero-based --sequence 1 --begin 2 --end 3 --comment=# ${prefix}.bed.gz
	"""

}

workflow.onComplete {
	if (workflow.success){
		subject = "QTLtools prep complete"
	}
	else {
		subject = "QTLtools prep error"
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

