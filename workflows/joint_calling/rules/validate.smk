rule gatk_concordance_NA12878:
    input:
        ref = rules.download_reference_genome.output,
        intervals = config["exon_bed_nochrX"],
        truth_vcf = config["NA12878_validation_vcf"],
        truth_idx = config["NA12878_validation_vcf"] + ".tbi",
        eval_vcf = rules.merge_snp_indel_vcf.output.vcf,
        eval_idx = rules.merge_snp_indel_vcf.output.idx
    output:
        summary = os.path.join(out_path, "concordance_NA12878/gatk/{group}.tsv"),
        tpfp = os.path.join(out_path, "concordance_NA12878/gatk/{group}_tpfp.vcf.gz"),
        tpfn = os.path.join(out_path, "concordance_NA12878/gatk/{group}_tpfn.vcf.gz")
    conda:
        "../envs/gatk.yaml"
    log:
        os.path.join(out_path, "log/concordance_NA12878/gatk/{group}_concordance.log")
    shell:
        """
        gatk Concordance \
            -R {input.ref} \
            -eval {input.eval_vcf} \
            --truth {input.truth_vcf} \
            --intervals {input.intervals} \
            --summary {output.summary} \
            -tpfp {output.tpfp} \
            -tpfn {output.tpfn} 2> {log}
        """
        
rule summarize_concordance_output_to_bed:
    input:
        tpfp = rules.gatk_concordance_NA12878.output.tpfp,
        tpfn = rules.gatk_concordance_NA12878.output.tpfn,
    output:
        fp = os.path.join(out_path, "concordance_NA12878/gatk/{group}_false_positive.bed"),
        fn = os.path.join(out_path, "concordance_NA12878/gatk/{group}_false_negative.bed")
    conda:
        "../envs/bwa.yaml"
    log:
        os.path.join(out_path, "log/concordance_NA12878/gatk/{group}_bed_summary.log")
    shell:
        """
        bcftools query -f "%CHROM\t%POS\t%INFO/STATUS\n" {input.tpfp} | \
         awk '{{OFS = "\t"}}{{if ($3 != "TP") print $1,$2-1,$2,$3}}' > {output.fp}
        
        bcftools query -f "%CHROM\t%POS\t%INFO/STATUS\n" {input.tpfn} | \
         awk '{{OFS = "\t"}}{{if ($3 != "TP") print $1,$2-1,$2,$3}}' > {output.fn}
        """

rule snpsift_concordance_NA12878:
    input:
        truth_vcf = config["NA12878_validation_vcf"],
        truth_idx = config["NA12878_validation_vcf"] + ".tbi",
        eval_vcf = rules.merge_snp_indel_vcf.output.vcf,
        eval_idx = rules.merge_snp_indel_vcf.output.idx
    output:
        truth_vcf = os.path.join(out_path, "concordance_NA12878/snpsift/truth_{group}.vcf"),
        eval_vcf = os.path.join(out_path, "concordance_NA12878/snpsift/eval_{group}.vcf"),
        summary = os.path.join(out_path, "concordance_NA12878/snpsift/{group}.tsv")
    conda:
        "../envs/snpsift.yaml"
    log:
        os.path.join(out_path, "log/concordance_NA12878/snpsift/{group}.log")
    shell:
        """
        zcat {input.truth_vcf} > {output.truth_vcf}
        zcat {input.eval_vcf} > {output.eval_vcf}

        SnpSift concordance \
            -v {output.truth_vcf} \
            {output.eval_vcf} > {output.summary} 2> {log}
        """

