rule genotype_single_gvcf:
	input:
		gvcf = rules.gatk_haplotypecaller.output.vcf,
		idx = rules.gatk_haplotypecaller.output.idx,
		ref = rules.download_reference_genome.output,
		exon_bed = config["exon_bed"]
	output:
		vcf = temp(os.path.join(out_path, "genotype_gvcf/{sample}.vcf.gz")),
		idx = temp(os.path.join(out_path, "genotype_gvcf/{sample}.vcf.gz.tbi"))
	conda:
		"../envs/gatk.yaml"
	log:
		os.path.join(out_path, "log/genotype_gvcf/{sample}.log")
	shell:
		"""
		gatk GenotypeGVCFs \
		        -R {input.ref} \
			-V {input.gvcf} \
		        -L {input.exon_bed} \
		        -O {output.vcf} 2> {log}
		"""
rule cnnscorevariants:
	input:
		vcf = rules.genotype_single_gvcf.output.vcf,
		idx = rules.genotype_single_gvcf.output.idx,
		bam = rules.gatk_applybqsr.output.bam,
		bai = rules.gatk_applybqsr.output.bai,
		ref = rules.download_reference_genome.output
	output:
		vcf = temp(os.path.join(out_path, "cnnscorevariants/{sample}.vcf.gz")),
		idx = temp(os.path.join(out_path, "cnnscorevariants/{sample}.vcf.gz.tbi")),
		_1 = touch(os.path.join(out_path, "cnnscorevariants/.temp_{sample}"))
	container:
		"docker://broadinstitute/gatk:4.1.3.0"
	log:
		os.path.join(out_path, "log/cnnscorevariants/{sample}.log")
	resources:
		mem_mb = 50000
	params:
		java_options = "-Xmx50G",
		tensor_type = "read_tensor"
	shell:
		"""
		gatk CNNScoreVariants \
			--java-options \"{params.java_options}\" \
		        -I {input.bam} \
		        -V {input.vcf} \
		        -R {input.ref} \
		        -O {output.vcf} \
			-tensor-type {params.tensor_type} 
		"""

# \"2> {log}\"
# 	conda:
# 		"../envs/gatk.yaml"

rule filtervarianttranches:
	input:
		vcf = rules.cnnscorevariants.output.vcf,
		idx = rules.cnnscorevariants.output.idx,
		resource1 = config["hapmap_resource"],
		resource2 = config["mills_resource"]
	output:
		vcf = temp(os.path.join(out_path, "filtervarianttranches/{sample}.vcf.gz")),
		idx = temp(os.path.join(out_path, "filtervarianttranches/{sample}.vcf.gz.tbi"))
	conda:
		"../envs/gatk.yaml"
	log:
		os.path.join(out_path, "log/filtervarianttranches/{sample}.log")
	params:
		info_key = "CNN_2D",
		snp_tranche = " 99.9 99.95",
		indel_tranche = " 99.4 99.0"
	shell:
		"""
		snp_tranche=$(for i in {params.snp_tranche} ; do echo " --snp-tranche $i" ; done | tr "\n" " ")
		indel_tranche=$(for i in {params.indel_tranche} ; do echo " --indel-tranche $i" ; done | tr "\n" " ")

		gatk FilterVariantTranches\
		        -V {input.vcf} \
			--resource {input.resource1} \
			--resource {input.resource2} \
			--info-key {params.info_key} \
			$snp_tranche \
			$indel_tranche \
		        -O {output.vcf} 2> {log}
		"""

rule funcotator:
	input:
		vcf = rules.filtervarianttranches.output.vcf,
		idx = rules.filtervarianttranches.output.idx,
		ref = rules.download_reference_genome.output,
		data_sources = config["funcotator_resource"]
	output:
		vcf = os.path.join(out_path, "funcotator/{sample}.vcf.gz"),
		idx = os.path.join(out_path, "funcotator/{sample}.vcf.gz.tbi")
	conda:
		"../envs/gatk.yaml"
	log:
		os.path.join(out_path, "log/funcotator/{sample}.log")
	params:
		ref_version = "hg38",
		out_file_format= "VCF"
	shell:
		"""
		gatk Funcotator \
			--variant {input.vcf} \
			--reference {input.ref} \
			--ref-version {params.ref_version} \
			--data-sources-path {input.data_sources} \
			--output {output.vcf} \
			--output-file-format {params.out_file_format} 2> {log}
		"""
