rule gatk_base_recalibrator:
	input:
		ref = rules.download_reference_genome.output,
		dict = rules.create_picard_index_dictionary.output, 
		bam = rules.sort_mark_duplicates.output.bam,
		bai = rules.sort_mark_duplicates.output.bai,
		dbsnp_ref = config["dbsnp_ref"],
		exon_bed = config["exon_bed"]
	output:
		temp(os.path.join(pre_path, "gatk_recalibrator/{sample}_recaldata.table"))
	conda:
		"../envs/gatk.yaml"
	log:
		os.path.join(pre_path, "log/gatk_base_recalibrator/{sample}.log")
	shell:
		"""
		gatk BaseRecalibrator \
		        -R {input.ref} \
		        -I {input.bam} \
		        --known-sites {input.dbsnp_ref} \
		        -L {input.exon_bed} \
		        -O {output} 2> {log}
		"""

rule gatk_applybqsr:
	input:
		bqsr_table = rules.gatk_base_recalibrator.output,
		ref = rules.download_reference_genome.output,
		bam = rules.sort_mark_duplicates.output.bam,
		bai = rules.sort_mark_duplicates.output.bai
	output:
		bam = temp(os.path.join(pre_path, "gatk_applybqsr/{sample}_recal.bam")),
		bai = temp(os.path.join(pre_path, "gatk_applybqsr/{sample}_recal.bai"))
	conda:
		"../envs/gatk.yaml"
	log:
		os.path.join(pre_path, "log/gatk_applybqsr/{sample}.log")
	shell:
		"""
		gatk ApplyBQSR \
		        -R {input.ref} \
		        -I {input.bam} \
		        -bqsr {input.bqsr_table} \
		        -O {output.bam} 2> {log}
		"""

# bam = lambda wildcards: expand(rules.gatk_applybqsr.output, sample = get_samples_given_group(wildcards.group)),
# bams=$(for bam in {input.bam} ; do echo "-I $bam" ; done | tr "\n" " ")

rule gatk_haplotypecaller:
	input:
		bam = rules.gatk_applybqsr.output.bam,
		bai = rules.gatk_applybqsr.output.bai,
		ref = rules.download_reference_genome.output,
		dbsnp_ref = config["dbsnp_ref"],
		exon_bed = config["exon_bed"]
	output:
		gvcf = os.path.join(pre_path, "gatk_haplotypecaller/{sample}_variants.g.vcf.gz"),
		idx = os.path.join(pre_path, "gatk_haplotypecaller/{sample}_variants.g.vcf.gz.tbi")
	conda:
		"../envs/gatk.yaml"
	log:
		os.path.join(pre_path, "log/gatk_haplotypecaller/{sample}.log")
	params:
		annotations = " QualByDepth",
		erc = "GVCF"
	shell:
		"""
		annotations=$(echo {params.annotations} | sed "s/ / -A /g" | sed "s/^/ -A /")

		gatk HaplotypeCaller \
		        -R {input.ref} \
		        -L {input.exon_bed} \
		        -I {input.bam} \
		        --dbsnp {input.dbsnp_ref} \
		        -O {output.gvcf} \
			$annotations \
			-ERC {params.erc} 2> {log}
		"""
