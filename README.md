## RNA-sequencing pipeline and gene expression analysis </br>
The aim of this project is to create an Nextflow pipeline for RNA seq analysis. The raw sequencing data were used from Guo et al. Nature Communications 2019 and can be found in the accession GSE106305. <br>
<h3>Pipeline summary</h3>

  1) Quality control (QC)<br/>
  ** FastQC ** : A Quality Control application for FastQ files
    FastQC is a program designed to spot potential problems in high througput sequencing datasets. It runs a set of analyses on one or more raw sequence files in fastq or bam         format and produces a report which summarises the results.
  2) Align to the reference genome<br/>
  3) Post-alignment QC<br/>
  4) Visualisation<br/>
  5) Quantify transcripts<br/>
  6) Visualise tracks against the reference genome<br/>
