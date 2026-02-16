# nf-data-prep-pipeline

A portable Nextflow pipeline to automate sequencing data retrieval, quality control, and genome indexing using Docker containers.

## 🚀 Overview

This pipeline follows a modular design:
1. **Data Import (Process 1)**: A unified process that fetches both sequencing reads (SRA/FASTQ) and genomic references (FASTA/GTF) from remote sources.
2. **Quality Control (Process 2)**: Standard FastQC assessment of the raw reads.
3. **STAR Indexing (Process 3)**: High-performance generation of the genome index using the files provided by the import step.

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
sra_id: "SRR12345678" # Supports SRA Accessions or http/ftp URLs
genome_url: "https://.../genome.fna.gz"
gtf_url: "https://.../genome.gtf.gz"
outdir: "results"
```

### 2. Resource Tuning (nextflow.config)

The configuration is pre-set for standard human genome indexing (~32GB RAM).  
You can adjust the cpus and memory in this file to fit your hardware.

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
