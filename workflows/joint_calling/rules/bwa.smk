rule build_index:
	input:
		rules.download_reference_genome.output
	output:
		bwa_index = directory(os.path.join(pre_path, "bwa_index")),
		fasta_index = os.path.join(pre_path, "GRCh38_full_analysis_set_plus_decoy_hla.fa.fai")
	conda:
		"../envs/bwa.yaml"
	log:
		os.path.join(pre_path, "log/build_index.log")
	shell:
		"""
		mkdir -p {output.bwa_index}

		bwa index -p {output.bwa_index}/bwa_index {input} 2>> {log}
		samtools faidx {input} 2>> {log}
		"""

rule create_picard_index_dictionary:
	input:
		rules.download_reference_genome.output
	output:
		os.path.join(pre_path, "GRCh38_full_analysis_set_plus_decoy_hla.dict")
	conda:
		"../envs/bwa.yaml"
	log:
		os.path.join(pre_path, "log/create_picard_index_dictionary.log")
	shell:
		"""
		picard CreateSequenceDictionary R={input} \
			O={output}
		"""

rule read_alignment:
	input:
		dir = os.path.join(config["fastq_dir_path"], "{sample}"),
		ref_genome = rules.download_reference_genome.output,
		ref_index = rules.build_index.output.bwa_index
	output:
		touch(os.path.join(pre_path, "read_alignment/{sample}/.DONE"))
	conda:
		"../envs/bwa.yaml"
	threads:
		4
	log:
		os.path.join(pre_path, "log/read_alignment/{sample}.log")
	shell:
		"""
		outdir=$(dirname {output})
		mkdir -p $outdir

		bn=$(for i in $(find {input.dir} | grep "gz") ; do basename $i | cut -d"_" -f1,2,3,4 ; done | sort -u)

		for i in $bn; do
			echo $i >> {log}
			r1={input.dir}/"$i"_1.fq.gz
			r2={input.dir}/"$i"_2.fq.gz

			tmp_fastq="$outdir"/.temp_{wildcards.sample}.fastq
			rm -rf $tmp_fastq 

			vsearch --fastx_subsample $r1 --fastqout $tmp_fastq --sample_pct 0.1 &>> {log}

			id=$(cat $tmp_fastq | awk '(NR == 1)' | cut -f 3-4 -d":" | sed "s/@//" | sed "s/:/_/g")
			lb={wildcards.sample}
			pl="ILLUMINA"
			sm={wildcards.sample}
			
			echo $id >> {log}
			echo $lb >> {log}
			echo $pl >> {log}
			echo $sm >> {log}

			bwa mem \
				-M \
				-t {threads} \
				-R "@RG\\tID:""$id""\\tLB:""$lb""\\tPL:""$pl""\\tSM:""$sm" \
				{input.ref_index}/bwa_index \
				$r1 $r2  2>> {log} | \
			samtools sort -@ {threads} -o $outdir/"$i".bam 2>> {log}
			samtools index -@ {threads} $outdir/"$i".bam
		done
		"""


rule merge_bamfiles:
	input:
		rules.read_alignment.output
	output:
		temp(os.path.join(pre_path, "merge_bamfiles/{sample}.bam"))
	conda:
		"../envs/bwa.yaml"
	log:
		os.path.join(pre_path, "log/merge_bamfiles/{sample}.log")
	threads: 16
	shell:
		"""
		indir=$(dirname {input})
		samtools merge -@ {threads} -f -o {output} $(find $indir | grep "bam$") 2> {log}
		"""
	
rule mark_duplicates:
	input:
		rules.merge_bamfiles.output
	output:
		bam = temp(os.path.join(pre_path, "marked_dup/{sample}.bam")),
		bai = temp(os.path.join(pre_path, "marked_dup/{sample}.bai")),
		metrics = temp(os.path.join(pre_path, "marked_dup/{sample}.metrics.txt"))
	conda:
		"../envs/bwa.yaml"
	log:
		os.path.join(pre_path, "log/mark_duplicates/{sample}.log")
	params:
		assume_sort_order = "coordinate",
		validation_stringency = "SILENT",
		remove_duplicates = "false",
		create_index = 	"true"
	shell:
		"""
		picard MarkDuplicates \
		 I={input} \
		 O={output.bam} \
		 ASSUME_SORT_ORDER={params.assume_sort_order} \
		 VALIDATION_STRINGENCY={params.validation_stringency} \
		 REMOVE_DUPLICATES={params.remove_duplicates} \
		 CREATE_INDEX={params.create_index} \
		 M={output.metrics} 2>> {log}
		"""

rule sort_mark_duplicates:
	input:
		rules.mark_duplicates.output.bam
	output:
		bam = temp(os.path.join(pre_path, "marked_dup/{sample}.sorted.bam")),
		bai = temp(os.path.join(pre_path, "marked_dup/{sample}.sorted.bai"))
	conda:
		"../envs/bwa.yaml"
	params:
		sort_order = "coordinate",
		create_index = 	"true"
	shell:
		"""
		picard SortSam \
		 I={input} \
    		 O={output.bam} \
    		 SORT_ORDER={params.sort_order} \
    		 CREATE_INDEX={params.create_index}
		"""	
		
