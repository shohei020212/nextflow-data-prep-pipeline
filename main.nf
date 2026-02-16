nextflow.enable.dsl=2

// --- Module 1: Data Import ---
process DOWNLOAD_DATA {
    publishDir "${params.outdir}/fastq", mode: 'copy', pattern: "*.fastq.gz"
    publishDir "${params.outdir}/reference", mode: 'copy', pattern: "genome.*"
    
    output:
    path "*.fastq.gz", emit: fastq
    path "genome.fna", emit: fasta
    path "genome.gtf", emit: gtf

    script:
    """
    # 1. Download Sequencing Data (FASTQ)
    if [[ "${params.sra_id}" =~ ^[SED]RR ]]; then
        prefetch ${params.sra_id}
        fasterq-dump ${params.sra_id} --split-files
        rm -rf ${params.sra_id}
    else
        wget ${params.sra_id} -O sample.fastq.gz
        gunzip sample.fastq.gz
    fi

    # 2. Download Reference Genome (FASTA)
    wget ${params.genome_url} -O genome.fna.gz
    gunzip genome.fna.gz

    # 3. Download Annotation (GTF)
    wget ${params.gtf_url} -O genome.gtf.gz
    gunzip genome.gtf.gz
    """
}

// --- Module 2: Quality Control ---
process FASTQC {
    publishDir "${params.outdir}/fastqc", mode: 'copy'

    input:
    path fastq

    output:
    path "*.{html,zip}"

    script:
    """
    fastqc ${fastq} -t ${task.cpus}
    """
}

// --- Module 3: Reference Indexing ---
process STAR_INDEX {
    publishDir "${params.outdir}/star_index", mode: 'copy'

    input:
    path fasta
    path gtf

    output:
    path "genome_index", emit: index

    script:
    """
    mkdir genome_index
    STAR --runMode genomeGenerate \\
         --genomeDir genome_index \\
         --genomeFastaFiles ${fasta} \\
         --sjdbGTFfile ${gtf} \\
         --runThreadN ${task.cpus}
    """
}

workflow {
    
    // 1. Download all required data
    DOWNLOAD_DATA()

    // 2. Perform QC on the downloaded FASTQ files
    FASTQC(DOWNLOAD_DATA.out.fastq.flatten())

    // 3. Build STAR index using the downloaded Reference files
    STAR_INDEX(DOWNLOAD_DATA.out.fasta, DOWNLOAD_DATA.out.gtf)
    
}