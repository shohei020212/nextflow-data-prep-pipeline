# nf-data-prep-pipeline

A portable Nextflow pipeline to automate sequencing data retrieval, quality control, and genome indexing using Docker containers.

## 🚀 Overview

This pipeline streamlines the initial steps of genomics workflows:
1. **Automated Fetching**: Downloads FASTQ files from NCBI SRA (via Accession ID) or direct URLs.
2. **Reference Management**: Retrieves Reference Genome (FASTA) and Annotation (GTF) files.
3. **Quality Check**: Runs FastQC on raw reads.
4. **Genome Indexing**: Generates a STAR index for downstream alignment.

## 🛠 Prerequisites

The only tools you need to have installed locally are:

* **Nextflow** (DSL2 enabled)
* **Docker** (Ensure the Docker Engine is running)

*Note: You do NOT need to install STAR, FastQC, or sra-tools on your host machine. Nextflow will automatically pull the required Docker images.*

## 📂 Project Structure

* `main.nf`: Core workflow logic.
* `nextflow.config`: Resource management and Docker configuration.
* `params.yaml`: User-defined parameters for data sources.

## ⚙️ Configuration

### 1. Set Parameters (`params.yaml`)
Update the values to match your specific study:

```yaml
sra_id: "SRR1234567" # Supports SRA Accessions or http/ftp URLs
genome_url: "https://.../genome.fna.gz"
gtf_url: "https://.../genome.gtf.gz"
outdir: "results"
```

### 2. Resource Tuning (nextflow.config)

The configuration is pre-set for standard human genome indexing (~32GB RAM). You can adjust the cpus and memory in this file to fit your hardware.

## 🏃 How to Run

Simply run the following command in your terminal:

```bash
nextflow run main.nf -params-file params.yaml
```

## 📊 Results

Processed files are organized in the `${params.outdir}` directory:

* `fastq/`: Raw sequencing data (FASTQ)
* `fastqc/`: HTML quality reports
* `reference/`: Downloaded FASTA and GTF files
* `star_index/`: Generated STAR genome index
