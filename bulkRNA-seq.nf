#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/****************************************************
 * Params
 ****************************************************/
params.outdir      = "./results"
params.input       = "/mnt/d/NEXTFLOW/bulk_RNAseq/data/fastq"
params.reference   = "Reference/genome.fa"
params.annotation  = "Reference/annotation.gtf"
params.genomedir   = "./genome_index"
//params.sra_script  = "SRA_download.py"        // your existing script (produces fastq/*.fastq.gz)


// Channel: grab all FASTQ files from input folder
	fastq_ch = Channel.fromPath("${params.input}/*.fastq.gz")


//====================================================
// Step 2: FastQC QC
//====================================================
process fastqc {
    tag { fastq.baseName }          // use file name as tag
    publishDir "${params.outdir}/fastqc", mode: 'copy'

    input:
    path fastq

    output:
    path("*_fastqc.html"), emit: html
    path("*_fastqc.zip"), emit: zip
    path("*_fastqc.*"), emit: reports

    script:
    """
    set -euo pipefail
    fastqc -t ${task.cpus} -o . ${fastq}
    """
}

//====================================================
// Step 3: Map sample IDs to input FASTQ files
//====================================================
samples_ch = Channel.from([
    tuple("LNCAP_Normoxia_S1", [
        file("${params.input}/SRR7179504.fastq.gz"),
        file("${params.input}/SRR7179505.fastq.gz"),
        file("${params.input}/SRR7179506.fastq.gz"),
        file("${params.input}/SRR7179507.fastq.gz")
    ]),
    tuple("LNCAP_Normoxia_S2", [
        file("${params.input}/SRR7179508.fastq.gz"),
        file("${params.input}/SRR7179509.fastq.gz"),
        file("${params.input}/SRR7179510.fastq.gz"),
        file("${params.input}/SRR7179511.fastq.gz")
    ]),
    tuple("LNCAP_Hypoxia_S1", [
        file("${params.input}/SRR7179520.fastq.gz"),
        file("${params.input}/SRR7179521.fastq.gz"),
        file("${params.input}/SRR7179522.fastq.gz"),
        file("${params.input}/SRR7179523.fastq.gz")
    ]),
    tuple("LNCAP_Hypoxia_S2", [
        file("${params.input}/SRR7179524.fastq.gz"),
        file("${params.input}/SRR7179525.fastq.gz"),
        file("${params.input}/SRR7179526.fastq.gz"),
        file("${params.input}/SRR7179527.fastq.gz")
    ]),
    tuple("PC3_Normoxia_S1", [ file("${params.input}/SRR7179536.fastq.gz") ]),
    tuple("PC3_Normoxia_S2", [ file("${params.input}/SRR7179537.fastq.gz") ]),
    tuple("PC3_Hypoxia_S1", [ file("${params.input}/SRR7179540.fastq.gz") ]),
    tuple("PC3_Hypoxia_S2", [ file("${params.input}/SRR7179541.fastq.gz") ])
])

//====================================================
// Step 4: Concatenate files into single sample FASTQ
//====================================================
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

//====================================================
// Step 5: STAR genome index (runs once)
//====================================================
process star_index {
    tag "STAR index"
	publishDir "${params.outdir}/star_index", mode: 'copy'
    
	input:
    path fasta
    path gtf

    output:
    path(params.genomedir), emit: index_dir

    script:
    """
    mkdir -p ${params.genomedir}
    STAR --runThreadN ${task.cpus} \
         --runMode genomeGenerate \
         --genomeDir ${params.genomedir} \
         --genomeFastaFiles ${fasta} \
         --sjdbGTFfile ${gtf} \
		 --sjdbOverhang 99
    """
}

/****************************************************
 * Step 6: STAR alignment (single-end)
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
    STAR --runThreadN 32 \
         --genomeDir ${index_dir} \
         --readFilesIn ${fastq_file} \
         --readFilesCommand zcat \
         --outFileNamePrefix ${sample_id}. \
         --outSAMtype BAM SortedByCoordinate \
         --quantMode GeneCounts
    """
}


//====================================================
// Workflow
//====================================================
workflow {

	// 2. FastQC
    fastqc(fastq_ch)
	
    // 3. Concatenate single-end files
    concatenated_ch = concatenate(samples_ch)
    
	// 4. STAR index (broadcast)
    index_dir_ch = star_index(file(params.reference), file(params.annotation)).index_dir.broadcast()

    // 5. STAR indexchannel of 8 samples
    samples_ch = Channel
        .fromPath(params.reads)
        .map { file -> 
            def sample_id = file.getBaseName().replaceFirst(/\.fastq(\.gz)?$/, "")
            tuple(sample_id, file)
        }
        .take(8)  // limit to 8 samples

    // broadcast STAR index dir (already built)
    index_dir_ch = Channel.value(file(params.genomedir))

    // run alignment
    (aligned_bams_ch, gene_counts_ch) = star_align(samples_ch, index_dir_ch)


}
