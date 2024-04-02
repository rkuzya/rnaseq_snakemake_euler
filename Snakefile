import glob 
import os
import logging

# Define your input files
#SAMPLES = ["A", "B"]

TRIMFILE="/cluster/home/rkuzyakiv/project_rostyk/adapters/adapters.fa"
SAMPLES, = glob_wildcards("input/{sample}.fastq.gz")
#GENOME="/cluster/home/rkuzyakiv/project_rostyk/hisat/grch38/genome"
print(SAMPLES)



localrules: all, run_fastqc, trim, hisat_map, samtools_convert, samtools_sort,samtools_index, count_genes_htseq 

# Define the rule all to specify the final output
rule all:
    input: 
        expand("output/{sample}_fastqc.html", sample=SAMPLES) +
        expand("output/{sample}_trimmed.fastq.gz", sample=SAMPLES)+
	expand("output/{sample}.sam", sample=SAMPLES)+
	expand("output/{sample}.bam", sample=SAMPLES)+
	expand("output/{sample}_sorted.bam", sample=SAMPLES)+
	expand("output/{sample}_sorted.bam.bai", sample=SAMPLES)+
        expand("output/{sample}_htseq_counts.txt", sample=SAMPLES)
	
# Rule for running FastQC
rule run_fastqc:
    input: "input/{sample}.fastq.gz"
    output: "output/{sample}_fastqc.html"
    singularity: "https://depot.galaxyproject.org/singularity/fastqc:0.11.9--0"
    shell: "fastqc {input} --outdir=output/"

# Rule for running TRIMMOMATIC
rule trim:
    input: "input/{sample}.fastq.gz"
    output: "output/{sample}_trimmed.fastq.gz"
    singularity: "https://depot.galaxyproject.org/singularity/trimmomatic:0.39--1"
    shell: "trimmomatic SE -phred33 {input} {output} ILLUMINACLIP:"+TRIMFILE+":2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"


# Rule for running HISAT2
rule hisat_map:
    input:
        "output/{sample}_trimmed.fastq.gz"
    output:
        "output/{sample}.sam"
    singularity:
### HAD TO DOWNLOAD THE SINGULARITY ### 
#       "https://depot.galaxyproject.org/singularity/hisat2:2.1.0--py37hc9558a2_4"
       "https://depot.galaxyproject.org/singularity/hisat2:2.2.1--he1b5a44_2"
#       "/cluster/home/rkuzyakiv/project_rostyk/myRNAseq/hisat2:2.1.0--py37hc9558a2_4.1"
    shell:
        "hisat2 -x /cluster/home/rkuzyakiv/project_rostyk/myRNAseq/hisat/grch38/genome -U {input} -S {output}"

# Rule for converting SAM to BAM
rule samtools_convert:
    input:
        "output/{sample}.sam"
    output:
        "output/{sample}.bam"
    singularity:
        "https://depot.galaxyproject.org/singularity/samtools:1.9--h91753b0_8"
    shell:
        "samtools view -@ 24 -bS {input} > {output} "

# Rule for sorting BAM file 
rule samtools_sort:
    input:
        "output/{sample}.bam"
    output:
        "output/{sample}_sorted.bam"
    singularity:
        "https://depot.galaxyproject.org/singularity/samtools:1.9--h91753b0_8"
    shell:
        "samtools sort -@ 24 -T ./  -O bam {input} > {output}"

# Rule for indexing 
rule samtools_index:
    input:
        "output/{sample}_sorted.bam"
    output:
        "output/{sample}_sorted.bam.bai"
    singularity:
        "https://depot.galaxyproject.org/singularity/samtools:1.9--h91753b0_8"
    shell:
        "samtools index {input}"

# Rule for counting the reads 
#rule count_genes:
#    input:
#        "output/{sample}_sorted.bam"
#    output:
#        "output/{sample}_sorted_counts.txt"
#    singularity:
#        "https://depot.galaxyproject.org/singularity/subread:2.0.6--he4a0461_0"
#    shell:
#        "featureCounts -M -f --fraction -s 2 -T 24 -t gene -g gene_id -a hg38_annotattion/Homo_sapiens.GRCh38.111.gtf -o {output} {input}"

# Rule for counting reads using HTSeq
rule count_genes_htseq:
    input:
        bam="output/{sample}_sorted.bam"
    output:
        counts="output/{sample}_htseq_counts.txt"
    singularity:
        "https://depot.galaxyproject.org/singularity/htseq:2.0.5--py39hd5189a5_0"
    params:
        gtf="/cluster/home/rkuzyakiv/project_rostyk/myRNAseq/hg38_annotattion/Homo_sapiens.GRCh38.111.gtf",
        strandedness="reverse", # Adjust as necessary: 'yes', 'no', or 'reverse'
        id_attribute="gene_id", # or use 'gene_id' depending on your preference
        feature_type="exon", # Typically 'exon' for gene-level counts
        mode="union" # 'union', 'intersection-strict', or 'intersection-nonempty'
    shell:
        """
        htseq-count -f bam -r pos -s {params.strandedness} -i {params.id_attribute} -t {params.feature_type} --additional-attr=gene_name -m {params.mode} {input.bam} {params.gtf} > {output.counts}
        """

