snakemake -j 50 -k --use-conda \
	--cluster-config ../../config/cluster.yaml \
	--conda-frontend mamba \
	--use-singularity \
	--singularity-prefix .snakemake/singularity \
	--singularity-args "--bind /fs03:/fs03" \
	--cluster "sbatch \
	--cpus-per-task={threads} \
	--output={cluster.output} \
	--error={cluster.error} \
	--job-name={cluster.jobname} \
	--mail-user=muhammad.zarulhanifah@monash.edu \
	{cluster.etc}"
