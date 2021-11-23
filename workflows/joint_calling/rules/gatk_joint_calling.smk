rule gatk_combine_gvcfs:
    input:
        gvcfs = lambda wildcards: expand(rules.gatk_haplotypecaller.output.gvcf, sample = get_samples_given_group(wildcards.group)),
        ref = rules.download_reference_genome.output,
        exon_bed = config["exon_bed"]
    output:
        gvcf = os.path.join(out_path, "gatk_combine/{group}.g.vcf.gz"),
        idx = os.path.join(out_path, "gatk_combine/{group}.g.vcf.gz.tbi")
    conda:
        "../envs/gatk.yaml"
    log:
        os.path.join(out_path, "log/gatk_combine_gvcfs/{group}.log")
    shell:
        """
        gvcfs=$(for i in {input.gvcfs} ; do echo "-V "$i ; done | tr "\\n" " ")

        gatk CombineGVCFs \
            -R {input.ref} \
            $gvcfs \
            -O {output.gvcf} \
            -L {input.exon_bed} 2> {log}
        """

rule gatk_genotype_combined_gvcf:
    input:
        gvcf = rules.gatk_combine_gvcfs.output.gvcf,
        ref = rules.download_reference_genome.output,
        exon_bed = config["exon_bed"]
    output:
        vcf = os.path.join(out_path, "gatk_genotype/{group}.vcf.gz"),
        idx = os.path.join(out_path, "gatk_genotype/{group}.vcf.gz.tbi")
    conda:
        "../envs/gatk.yaml"
    log:
        os.path.join(out_path, "log/gatk_genotype/{group}.log")
    shell:
        """
        gatk GenotypeGVCFs \
            -R {input.ref} \
            -V {input.gvcf} \
            -O {output.vcf} \
            -L {input.exon_bed} 2> {log}
        """
            
rule gatk_selectvariants_snp:
    input:
        vcf = rules.gatk_genotype_combined_gvcf.output.vcf,
        ref = rules.download_reference_genome.output,
        exon_bed = config["exon_bed"]
    output:
        vcf = os.path.join(out_path, "gatk_selectvariants_snp/{group}.vcf.gz"),
        idx = os.path.join(out_path, "gatk_selectvariants_snp/{group}.vcf.gz.tbi")
    conda:
        "../envs/gatk.yaml"
    log:
        os.path.join(out_path, "log/gatk_selectvariants_snp/{group}.log")
    params:
        select_type_to_include = "SNP"
    shell:
        """
        gatk SelectVariants \
            -R {input.ref} \
            -V {input.vcf} \
            -O {output.vcf} \
            -L {input.exon_bed} \
            --select-type-to-include {params.select_type_to_include} 2> {log}
        """

rule gatk_variant_recalibrator_snp:
    input:
        ref = rules.download_reference_genome.output,
        vcf = rules.gatk_selectvariants_snp.output.vcf,
        exon_bed = config["exon_bed"],
        resource1 = config["hapmap_resource"],
        resource2 = config["omni_resource"],
        resource3 = config["1000G_resource"],
        resource4 = config["dbsnp_resource"]
    output:
        recal = os.path.join(out_path, "gatk_variant_recalibrator_snp/{group}_snp.recal"),
        tranche_file = os.path.join(out_path, "gatk_variant_recalibrator_snp/{group}_snp.tranches"),
        fig = os.path.join(out_path, "gatk_variant_recalibrator_snp/{group}.plots.R")
    conda:
        "../envs/gatk.yaml"
    log:
        os.path.join(out_path, "log/gatk_variant_recalibrator_snp/{group}.log")
    resources:
        mem_mb = 50000
    params:
        java_options = "-DGATK_STACKTRACE_ON_USER_EXCEPTION=true -Xmx50G",
        tranches = " 100.0 99.95 99.9 99.8 99.6 99.5 99.4 99.3 99.0 98.0 97.0 90.0",
        annots = " QD FS SOR MQRankSum ReadPosRankSum MQ",
        mode = "SNP",
        max_gaussians = 4,
        max_attempts = 5,
        resource1 = "hapmap,known=false,training=true,truth=true,prior=15",
        resource2 = "omni,known=false,training=true,truth=false,prior=12",
        resource3 = "1000G,known=false,training=true,truth=false,prior=10",
        resource4 = "dbsnp,known=true,training=false,truth=false,prior=2"
    shell:
        """
        tranches=$(echo {params.tranches} | sed "s/ / -tranche /g" | sed "s/^/ -tranche /")
        annots=$(echo {params.annots} | sed "s/ / -an /g" | sed "s/^/ -an /")

        gatk VariantRecalibrator \
            --java-options \"{params.java_options}\" \
            -R {input.ref} \
            -V {input.vcf} \
            -O {output.recal} \
            -L {input.exon_bed} \
            --resource:{params.resource1} {input.resource1} \
            --resource:{params.resource2} {input.resource2} \
            --resource:{params.resource3} {input.resource3} \
            --resource:{params.resource4} {input.resource4} \
            $tranches \
            $annots \
            --mode {params.mode} \
            --max-gaussians {params.max_gaussians} \
            --rscript-file {output.fig} \
            --tranches-file {output.tranche_file} 2> {log}
        """

rule gatk_applyVQSR_snp:
    input:
        ref = rules.download_reference_genome.output,
        vcf = rules.gatk_selectvariants_snp.output.vcf,
        exon_bed = config["exon_bed"],
        recal = rules.gatk_variant_recalibrator_snp.output.recal,
        tranch = rules.gatk_variant_recalibrator_snp.output.tranche_file
    output:
        vcf = os.path.join(out_path, "gatk_VQSR_snp/{group}.vcf.gz"),
        idx = os.path.join(out_path, "gatk_VQSR_snp/{group}.vcf.gz.tbi")
    conda:
        "../envs/gatk.yaml"
    log:
        os.path.join(out_path, "log/gatk_applyVQSR_snp/{group}.log")
    params:
        truth_sensitivity_filter_level = 99.7,
        create_output_variant_index = "true",
        mode = "SNP"
    shell:
        """
        gatk ApplyVQSR \
            -R {input.ref} \
            -V {input.vcf} \
            -O {output.vcf} \
            -L {input.exon_bed} \
            --recal-file {input.recal} \
            --tranches-file {input.tranch} \
            --truth-sensitivity-filter-level {params.truth_sensitivity_filter_level} \
            --create-output-variant-index {params.create_output_variant_index} \
            -mode {params.mode} 2> {log}
       """

rule gatk_selectvariants_indel:
    input:
        vcf = rules.gatk_genotype_combined_gvcf.output.vcf,
        ref = rules.download_reference_genome.output,
        exon_bed = config["exon_bed"]
    output:
        vcf = os.path.join(out_path, "gatk_selectvariants_indel/{group}.vcf.gz"),
        idx = os.path.join(out_path, "gatk_selectvariants_indel/{group}.vcf.gz.tbi")
    conda:
        "../envs/gatk.yaml"
    log:
        os.path.join(out_path, "log/gatk_selectvariants_indel/{group}.log")
    params:
        select_type_to_include = "INDEL"
    shell:
        """
        gatk SelectVariants \
            -R {input.ref} \
            -V {input.vcf} \
            -O {output.vcf} \
            -L {input.exon_bed} \
            --select-type-to-include {params.select_type_to_include} 2> {log}
        """

rule gatk_variant_recalibrator_indel:
    input:
        ref = rules.download_reference_genome.output,
        vcf = rules.gatk_selectvariants_indel.output.vcf,
        exon_bed = config["exon_bed"],
        resource1 = config["mills_resource"],
        resource2 = config["axiompoly_resource"],
        resource3 = config["dbsnp_resource"]
    output:
        recal = os.path.join(out_path, "gatk_variant_recalibrator_indel/{group}.recal"),
        tranche_file = os.path.join(out_path, "gatk_variant_recalibrator_indel/{group}.tranches"),
        fig = os.path.join(out_path, "gatk_variant_recalibrator_indel/{group}.plots.R")
    conda:
        "../envs/gatk.yaml"
    log:
        os.path.join(out_path, "log/gatk_variant_recalibrator_indel/{group}.log")
    params:
        java_options = "-Xmx50G",
        tranches = " 100.0 99.95 99.9 99.8 99.6 99.5 99.4 99.3 99.0 98.0 97.0 90.0",
        annots = " QD FS SOR MQ MQRankSum ReadPosRankSum",
        mode = "INDEL",
        max_gaussians = 4,
        max_attempts = 5,
        resource1 = "mills,known=false,training=true,truth=true,prior=12",
        resource2 = "axiomPoly,known=false,training=true,truth=false,prior=10",
        resource3 = "dbsnp,known=true,training=false,truth=false,prior=2"
    resources:
        mem_mb = 50000
    shell:
        """
        tranches=$(for i in {params.tranches}; do echo " -tranche "$i ; done | tr "\n" " ")
        annots=$(for i in {params.annots}; do echo " -an "$i ; done | tr "\n" " " )

        gatk VariantRecalibrator \
            --java-options \"{params.java_options}\" \
            -R {input.ref} \
            -V {input.vcf} \
            -O {output.recal} \
            -L {input.exon_bed} \
            --max-attempts {params.max_attempts} \
            --resource:{params.resource1} {input.resource1} \
            --resource:{params.resource2} {input.resource2} \
            --resource:{params.resource3} {input.resource3} \
            $tranches \
            $annots \
            --mode {params.mode} \
            --max-gaussians {params.max_gaussians} \
            --rscript-file {output.fig} \
            --tranches-file {output.tranche_file} 2> {log}
        """

rule gatk_applyVQSR_indel:
    input:
        ref = rules.download_reference_genome.output,
        vcf = rules.gatk_selectvariants_indel.output.vcf,
        exon_bed = config["exon_bed"],
        recal = rules.gatk_variant_recalibrator_indel.output.recal,
        tranch = rules.gatk_variant_recalibrator_indel.output.tranche_file,
    output:
        vcf = os.path.join(out_path, "gatk_VQSR_indel/{group}.vcf.gz"),
        idx = os.path.join(out_path, "gatk_VQSR_indel/{group}.vcf.gz.tbi")
    conda:
        "../envs/gatk.yaml"
    log:
        os.path.join(out_path, "log/gatk_applyVQSR_indel/{group}.log")
    params:
        truth_sensitivity_filter_level = 99.7,
        create_output_variant_index = "true",
        mode = "INDEL"
    shell:
        """
        gatk ApplyVQSR \
            -R {input.ref} \
            -V {input.vcf} \
            -O {output.vcf} \
            -L {input.exon_bed} \
            --recal-file {input.recal} \
            --tranches-file {input.tranch} \
            --truth-sensitivity-filter-level {params.truth_sensitivity_filter_level} \
            --create-output-variant-index {params.create_output_variant_index} \
            -mode {params.mode} 2> {log}
        """

rule merge_snp_indel_vcf:
    input:
        snp = rules.gatk_applyVQSR_snp.output.vcf,
        indel = rules.gatk_applyVQSR_indel.output.vcf,
    output:
        vcf = os.path.join(out_path, "merged_vcf/{group}.vcf.gz"),
        idx = os.path.join(out_path, "merged_vcf/{group}.vcf.gz.tbi"),
    conda:
        "../envs/bwa.yaml"
    log:
        os.path.join(pre_path, "log/merge_snp_indel_vcf/{group}.log")
    shell:
        """
        picard MergeVcfs \
            I={input.snp} \
            I={input.indel} \
            O={output.vcf} 2> {log}
        """
