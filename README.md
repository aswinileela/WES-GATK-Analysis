# WES-GATK-Analysis

Snakemake workflow for whole-exome sequencing (WES) analysis using GATK4. Supports both single-sample and joint genotyping variant calling, with validation against the NA12878 benchmark.

## Pipeline Overview

```
Raw FASTQ → Alignment (BWA) → Preprocessing (GATK)
→ Variant Calling (HaplotypeCaller / GenomicsDBImport)
→ Annotation (Funcotator) → Validation (NA12878 concordance)
```

## Workflows

| Workflow | Description |
|---|---|
| `workflows/single_calling/` | Per-sample HaplotypeCaller variant calling |
| `workflows/joint_calling/` | Joint genotyping across multiple samples |

## Installation

```bash
# Requires Conda and Snakemake
conda install -c conda-forge -c bioconda snakemake

# Create environments (Snakemake handles this automatically)
# See workflows/*/envs/*.yaml for software details
```

## Configuration

Edit `config/config.yaml`:
- Set `fastq_dir_path` to your input FASTQ directory
- Set `metadata` to your sample TSV
- Set paths for reference genome and GATK bundle resources
- Download GRCh38 reference and GATK bundle from Broad Institute

## Usage

```bash
# Single-sample calling
cd workflows/single_calling
snakemake --use-conda --cores 16

# Submit to SLURM cluster
sbatch run_pipeline.slurm

# Joint calling
cd workflows/joint_calling
snakemake --use-conda --cores 16
```

## Input

A `metadata.tsv` file listing sample names and groups:

```
sample  group
SAMPLE1 group_A
SAMPLE2 group_A
```

## Citation

If you use this workflow, please cite:

```
Leela, A. (2023). WES-GATK-Analysis. https://github.com/aswinileela/WES-GATK-Analysis
```

## License

GPL-3.0 — see [LICENSE.md](LICENSE.md)
