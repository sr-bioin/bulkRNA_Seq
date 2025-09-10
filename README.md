## RNA-sequencing pipeline and gene expression analysis </br>
The aim of this project is to create an Nextflow pipeline for RNA seq analysis. The raw sequencing data were used from Guo et al. Nature Communications 2019 and can be found in the accession GSE106305. <br>
<h3>Pipeline summary</h3>

  1) Quality control (QC)<br/>
    **FastQC**: It is a Quality Control application for FastQ files. 
    It is a program designed to spot potential problems in high througput sequencing datasets. It runs a set of analyses on one or more raw sequence files in fastq or bam format and produces a report which summarises the results. More information can be obtaiined from the website (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
  2) Align to the reference genome<br/>
    **STAR**: It aligns reads by finding the Maximal Mappable Prefix (MMP) hits between reads (or read pairs) and the genome, using a Suffix Array index. Different parts of a read can be mapped to different genomic positions, corresponding to splicing or RNA-fusions. More information can be found at https://github.com/alexdobin/STAR.
  3) Post-alignment QC<br/>
    **QualiMap**: It is a multi-threaded application built in Java and R that provides a graphical user interface to perform the quality control of alignment sequencing data. More information can be found at http://qualimap.conesalab.org/.
  
  4) Differential expression <br/>
    **DESeq2**: DESeq2 is a statistical tool that allows researchers to identify differentially expressed genes from RNA-Seq data. More information can be found at https://github.com/thelovelab/DESeq2 <br/>
     <img width="500" height="400" alt="LNCAP_volcano_plot" src="https://github.com/user-attachments/assets/6c8ce0d8-140e-4627-8e72-83f409aa2ce4" />

