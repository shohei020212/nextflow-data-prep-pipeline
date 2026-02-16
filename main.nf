nextflow.enable.dsl=2

// --- Module 1: Data Import ---
process DOWNLOAD_DATA {
    publishDir "${params.outdir}/fastq", mode: 'copy', pattern: "*.fastq"
    publishDir "${params.outdir}/reference", mode: 'copy', pattern: "genome.*"
    
    output:
    path "*.fastq", emit: fastq
    path "genome.fna", emit: fasta
    path "genome.gtf", emit: gtf

    script:
    """
    # Remove fragments like # text=... from URLs just in case
    CLEAN_GENOME_URL=\$(echo "${params.genome_url}" | cut -d'#' -f1)
    CLEAN_GTF_URL=\$(echo "${params.gtf_url}" | cut -d'#' -f1)

    # 1. Download Reference Genome (FASTA)
    wget "\$CLEAN_GENOME_URL" -O genome.fna.gz
    gunzip -f genome.fna.gz

    # 2. Download Annotation (GTF)
    wget "\$CLEAN_GTF_URL" -O genome.gtf.gz
    gunzip -f genome.gtf.gz

    # 3. Download Sequencing Data (FASTQ)
    if [[ "${params.sra_id}" =~ ^[SED]RR ]]; then
        prefetch ${params.sra_id}
        fasterq-dump ${params.sra_id} --split-files
        rm -rf ${params.sra_id}
    else
        wget ${params.sra_id} -O sample.fastq.gz
        gunzip -f sample.fastq.gz
    fi
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