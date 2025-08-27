#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/****************************************************
 * Params
 ****************************************************/
params.accessions = [
    "SRR7179504","SRR7179505","SRR7179506","SRR7179507",
    "SRR7179508","SRR7179509","SRR7179510","SRR7179511",
    "SRR7179520","SRR7179521","SRR7179522","SRR7179523",
    "SRR7179524","SRR7179525","SRR7179526","SRR7179527",
    "SRR7179536","SRR7179537","SRR7179540","SRR7179541"
]

params.sampleMap = [
    "SRR7179504": "LNCAP_Normoxia_S1", "SRR7179505": "LNCAP_Normoxia_S1",
    "SRR7179506": "LNCAP_Normoxia_S1", "SRR7179507": "LNCAP_Normoxia_S1",
    "SRR7179508": "LNCAP_Normoxia_S2", "SRR7179509": "LNCAP_Normoxia_S2",
    "SRR7179510": "LNCAP_Normoxia_S2", "SRR7179511": "LNCAP_Normoxia_S2",
    "SRR7179520": "LNCAP_Hypoxia_S1", "SRR7179521": "LNCAP_Hypoxia_S1",
    "SRR7179522": "LNCAP_Hypoxia_S1", "SRR7179523": "LNCAP_Hypoxia_S1",
    "SRR7179524": "LNCAP_Hypoxia_S2", "SRR7179525": "LNCAP_Hypoxia_S2",
    "SRR7179526": "LNCAP_Hypoxia_S2", "SRR7179527": "LNCAP_Hypoxia_S2",
    "SRR7179536": "PC3_Normoxia_S1", "SRR7179537": "PC3_Normoxia_S2",
    "SRR7179540": "PC3_Hypoxia_S1", "SRR7179541": "PC3_Hypoxia_S2"
]

params.reference = "Reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
params.annotation = "Reference/Homo_sapiens.GRCh38.99.gtf"
params.genomedir = "Reference/Homo_sapiens.GRCh38.99_index"
params.outdir = "./results"

/****************************************************
 * Step 1: Download + convert SRA → FASTQ
 ****************************************************/
process sra_download {
    tag "$sra_id"
    publishDir "${params.outdir}/fastq", mode: 'copy'
    
    input:
    val sra_id
    
    output:
    path("${sra_id}*.fastq.gz")
    
    script:
    """
    set -euo pipefail
    module load sratoolkit/3.0.0
	
    echo "Downloading: ${sra_id}"
    prefetch ${sra_id}
    echo "Converting to FASTQ: ${sra_id}"
    fasterq-dump --outdir . --split-3 ${sra_id}/${sra_id}.sra
    echo "Compressing: ${sra_id}"
    gzip ${sra_id}*.fastq
    """
}

/****************************************************
 * Step 2: FastQC
 ****************************************************/
process fastqc {
    tag { fastq.baseName }
    publishDir "${params.outdir}/fastqc", mode: 'copy'
    
    input:
    path fastq
    
    output:
    path("*_fastqc.html"), emit: html
    path("*_fastqc.zip"), emit: zip
    
    script:
    """
	module load fastqc
	
    set -euo pipefail
    fastqc -t ${task.cpus} -o . ${fastq}
    """
}

/****************************************************
 * Step 3: Concatenate replicates
 ****************************************************/
process concatenate {
    tag "$sample_id"
    publishDir "${params.outdir}/concatenate", mode: 'copy'

    input:
    tuple val(sample_id), path(fastqs)

    output:
    tuple val(sample_id), path("${sample_id}.fastq.gz")

    script:
    """
    set -euo pipefail
    cat ${fastqs.join(' ')} > ${sample_id}.fastq.gz
    """
}

/****************************************************
 * Step 4: STAR index
 ****************************************************/
/*process star_index {
    tag "STAR index"
    
    input:
    path fasta
    path gtf
    
    output:
    // path(params.genomedir), emit: index_dir
	path "Homo_sapiens.GRCh38.99_index", emit: index_dir
    
	// Only responsible for publishing files to your desired folder
    publishDir "Reference", mode: 'copy'
	
	
	script:
    """
    mkdir -p ${params.genomedir}
    STAR --runThreadN 32 \\
         --runMode genomeGenerate \\
         --genomeDir ${params.genomedir} \\
         --genomeFastaFiles ${fasta} \\
         --sjdbGTFfile ${gtf} \\
         --sjdbOverhang 99
    """
}*/

/****************************************************
 * Step 5: STAR alignment
 ****************************************************/
process star_align {
    tag "$sample_id"
    publishDir "${params.outdir}/star", mode: 'copy'
    
    input:
    tuple val(sample_id), path(fastq_file)
    path index_dir
    
    output:
    tuple val(sample_id), path("${sample_id}.Aligned.sortedByCoord.out.bam")
    tuple val(sample_id), path("${sample_id}.ReadsPerGene.out.tab")
    
    script:
    """
    STAR --runThreadN 16 \\
         --genomeDir ${index_dir} \\
         --readFilesIn ${fastq_file} \\
         --readFilesCommand zcat \\
         --outFileNamePrefix ${sample_id}. \\
         --outSAMtype BAM SortedByCoordinate \\
         --quantMode GeneCounts TranscriptomeSAM \\
		 --outSAMattributes NH HI AS NM MD
    """
}

/****************************************************
 * Step 6: Qualimap
 ****************************************************/
/*
process qualimap {
    tag { bam.baseName }
    publishDir "${params.outdir}/qualimap", mode: 'copy'
    
    input:
    tuple val(sample_id), path(bam)
    
    output:
    tuple val(sample_id), path("${sample_id}_qualimap")
    
    script:
    """
    set -euo pipefail
    module load qualimap

    mkdir ${sample_id}_qualimap
    qualimap bamqc \
        -bam ${bam} \
        -outdir ${sample_id}_qualimap \
        -outformat PDF:HTML
    """
}

/****************************************************
 * Workflow
 ****************************************************/
workflow {
	// Step 1: SRA download
    accession_ch = Channel.fromList(params.accessions)
    fastq_ch = sra_download(accession_ch)
	
	// Step 2. FastQC
    fastqc(fastq_ch)
	
	// Group FASTQs into samples
samples_ch = fastq_ch
    .map { fq ->
        def id = fq.getName().replaceAll(/\.fastq(\.gz)?$/, "")
        def sampleId = params.sampleMap[id]
        if (!sampleId) throw new IllegalStateException("No mapping for accession ${id}")
        tuple(sampleId, fq)
    }
    .groupTuple()

    // Step 3: Concatenate replicates
    concatenated_ch = concatenate(samples_ch)

    // Step 4: Wrap reference + annotation in channels, build STAR index once
    reference_ch  = Channel.value(file(params.reference))
    annotation_ch = Channel.value(file(params.annotation))
	index_dir_ch  = star_index(reference_ch, annotation_ch)													   
   
    // Step 5: Align concatenated FASTQs
    star_align(concatenated_ch, index_dir_ch)
	
	// Step 6: Qualimap QC
    qualimap(aligned_ch)
	 
}
