Code for key analyses in the manuscript Single-nucleus chromatin and gene expression profiling across hundreds of skeletal muscle samples reveals context-specific regulation. Annotations in hg38.

# To run workflow:
## Specify all paths in config.yaml
## Modify scripts/nextflow.config accordingly. It currently provides configuration to run on the slurm scheduler. 
## To generate the sigularity container, use the `Singularity` file in the scripts directory to generate one as:
```
singularity build --remote container_name.sig scripts/Singularity
```
## Run the Nextflow workflow as:
```
nextflow -C scripts/nextflow.config run scripts/main.nf -params-file config.yaml -resume -with-trace trace.txt
```

	