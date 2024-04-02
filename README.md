**RNA-Seq Workflow using Snakemake and Singularity Containers (ETH Zurich Euler Cluster Adapted)**

This repository contains a Snakemake workflow for RNA-Seq analysis, utilizing Singularity containers for tool dependencies. The workflow includes trimming, mapping, sorting, and summarizing steps, along with generating count tables from aligned reads.

### Tools Used:
- Trimmomatic for quality trimming
- HISAT2 for read mapping
- SAMtools for SAM/BAM file manipulation
- HTSeq for counting reads mapped to genes

### Running the Workflow:
To run the workflow, execute the following command in the terminal:

### to test
```bash
snakemake -p --cores 1 --use-singularity
```
### to run 
```bash
snakemake -p --use-singularity --cluster "sbatch -N 1 -n 1 -t 60" --jobs 1
```

### Workflow Steps:
1. **Trimming**: Utilizes Trimmomatic for quality trimming of input fastq files.
2. **Mapping**: Employs HISAT2 for aligning trimmed reads to the reference genome.
3. **Sorting**: Uses SAMtools to convert and sort SAM files into BAM format.
4. **Counting Reads**: Counts reads mapped to genes using Subread's HTSeq.
5. **Generating Count Table**: Aggregates count files into a single CSV file for downstream analysis.

### Input and Output:
- Input files are expected to be in fastq.gz format located in the `input/` directory.
- Output files are organized in the 'output/ directory.

### Singularity Containers:
Singularity containers are utilized for each tool dependency, ensuring reproducibility and portability of the workflow.

For more detailed information on each rule and its functionality, please refer to the Snakefile in this repository.
