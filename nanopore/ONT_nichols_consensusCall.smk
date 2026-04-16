import os
import pandas as pd
from prepare_ref_helpers import get_genes_from_files  # Assume you moved the function to prepare_ref_helpers.py

# Constants
WORK_DIR = "/scratch/2024-08-01_TPal_ONT"
BED_DIR = f"{WORK_DIR}/bed/nichols"
REFERENCE="/Research/References/TPallidum/Pallidum_Nichols/GCA_000410535.2_ASM41053v2/GCA_000410535.2_ASM41053v2_ref/GCA_000410535.2_ASM41053v2_genomic.fa"
BED_FULL_GENOME = f"{WORK_DIR}/CP004010.2.bed"
ASSEMBLY_DIR = f"{WORK_DIR}/assembly_results/nichols"

# Load sample list
SampleList = pd.read_csv(f"{WORK_DIR}/list/nichols.sampleList.txt", header=None)
initial_samples = list(SampleList[0])

# Get unique genes from all input files
input_files = [f"{BED_DIR}/{sample}.bed" for sample in initial_samples]
initial_genes = get_genes_from_files(input_files)
#######################################################
# Define a function to filter valid samples/genes
import logging

def valid_samples_genes(samples, genes, assembly_dir):
    valid = []
    for sample in samples:
        for gene in genes:
            file_path = f"{assembly_dir}/{sample}/assembly/canu/{gene}/{sample}.{gene}.correctedReads.fasta.gz"
            if os.path.exists(file_path):
                valid.append((sample, gene))
    return valid

# List of valid sample/gene pairs
VALID_PAIRS = valid_samples_genes(initial_samples, initial_genes, ASSEMBLY_DIR)

rule all:
    input:
        expand(
            f"{ASSEMBLY_DIR}/{{sample}}/assembly/medaka/{{sample}}.{{gene}}.m.fasta",
            zip, 
            sample=[pair[0] for pair in VALID_PAIRS], 
            gene=[pair[1] for pair in VALID_PAIRS]
        )

#######################################################
##### 1. Polish with Racon ############################
#######################################################
DUMMY_FILE = f"{ASSEMBLY_DIR}/{{sample}}/assembly/canu/{{gene}}/{{sample}}.{{gene}}.finished.txt"
rule check_contigs:
    group: "rule_post_canu"
    input:
        contigs = lambda wildcards: f"{ASSEMBLY_DIR}/{wildcards.sample}/assembly/canu/{wildcards.gene}/{wildcards.sample}.{wildcards.gene}.contigs.fasta" if os.path.exists(f"{ASSEMBLY_DIR}/{wildcards.sample}/assembly/canu/{wildcards.gene}/{wildcards.sample}.{wildcards.gene}.contigs.fasta") else DUMMY_FILE,
        trimmed_reads = lambda wildcards: f"{ASSEMBLY_DIR}/{wildcards.sample}/assembly/canu/{wildcards.gene}/{wildcards.sample}.{wildcards.gene}.trimmedReads.fasta.gz" if os.path.exists(f"{ASSEMBLY_DIR}/{wildcards.sample}/assembly/canu/{wildcards.gene}/{wildcards.sample}.{wildcards.gene}.trimmedReads.fasta.gz") else DUMMY_FILE,
        corrected_reads = lambda wildcards: f"{ASSEMBLY_DIR}/{wildcards.sample}/assembly/canu/{wildcards.gene}/{wildcards.sample}.{wildcards.gene}.correctedReads.fasta.gz" if os.path.exists(f"{ASSEMBLY_DIR}/{wildcards.sample}/assembly/canu/{wildcards.gene}/{wildcards.sample}.{wildcards.gene}.correctedReads.fasta.gz") else DUMMY_FILE
    output:
        contigs = f"{ASSEMBLY_DIR}/{{sample}}/assembly/canu/{{gene}}/{{sample}}.{{gene}}.contigs.temp.fasta",
        warning = f"{ASSEMBLY_DIR}/{{sample}}/assembly/canu/{{gene}}/{{sample}}.{{gene}}.contigs.warnings.txt"
    resources:
        time_min=5, mem_mb=2000, cores=1
    run:
        contigs_exist = input.contigs != DUMMY_FILE and os.path.getsize(input.contigs) > 0
        trimmed_reads_exist = input.trimmed_reads != DUMMY_FILE and os.path.getsize(input.trimmed_reads) > 0
        corrected_reads_exist = input.corrected_reads != DUMMY_FILE and os.path.getsize(input.corrected_reads) > 0

        if contigs_exist:
            with open(output.warning, "a") as f:
                f.write(f"{input.contigs} contain contigs and that's GOOD\n")
            shell("cp {input.contigs} {output.contigs}")
        elif trimmed_reads_exist:
            with open(output.warning, "a") as f:
                f.write(f"{input.trimmed_reads} is empty or missing. Contig will be represented by the longest trimmed sequence!\n")
            shell("gunzip -c {input.trimmed_reads} | paste - - | awk -F '\\t' '{{L=length($2);if(L>M) {{M=L;R=$0;}}}} END {{print R;}}' | tr '\\t' '\\n' > {output.contigs}")
        elif corrected_reads_exist:
            with open(output.warning, "a") as f:
                f.write(f"{input.corrected_reads} is empty or missing. Contigs will be represented by the longest corrected sequence!\n")
            shell("gunzip -c {input.corrected_reads} | paste - - | awk -F '\\t' '{{L=length($2);if(L>M) {{M=L;R=$0;}}}} END {{print R;}}' | tr '\\t' '\\n' > {output.contigs}")
        else:
            with open(output.warning, "a") as f:
                f.write(f"All input files are empty or missing. Contigs will NOT BE PRESENT!\n")
            with open(f"{ASSEMBLY_DIR}/{wildcards.sample}/assembly/canu/{wildcards.gene}/{wildcards.sample}.{wildcards.gene}.warning.txt", "w") as f:
                f.write(f"contigs for {ASSEMBLY_DIR}/{wildcards.sample}/assembly/canu/{wildcards.gene}/{wildcards.sample}.{wildcards.gene} are NOT PRESENT\n")
 
rule modify_contigs_header:
    group: "rule_post_canu"
    input:
        contigs = f"{ASSEMBLY_DIR}/{{sample}}/assembly/canu/{{gene}}/{{sample}}.{{gene}}.contigs.temp.fasta",
    output:
        contigs = f"{ASSEMBLY_DIR}/{{sample}}/assembly/canu/{{gene}}/{{sample}}.{{gene}}.header_md.contigs.fasta"
    resources: time_min=5, mem_mb=2000, cores=1
    shell:
        """
        cp {input.contigs} {output.contigs}
        sed -i '1s/.*/>{wildcards.sample}.{wildcards.gene}.contigs.fasta/' {output.contigs}
        """

rule minimap2_r1_assembly:
    group: "round_1"
    input:
        contigs = f"{ASSEMBLY_DIR}/{{sample}}/assembly/canu/{{gene}}/{{sample}}.{{gene}}.header_md.contigs.fasta",
        filtered_fastq = f"{ASSEMBLY_DIR}/{{sample}}/long_read/{{sample}}.{{gene}}.length_filtered.fastq"
    output:
        sam = temp(f"{ASSEMBLY_DIR}/{{sample}}/assembly/racon/{{sample}}.{{gene}}.sorted.mini_1.sam")
    resources: time_min=60, mem_mb=16000, cores=4
    shell:
        """
        module load minimap2/2.26 samtools/1.17
        minimap2 -ax map-ont -t 2 {input.contigs} {input.filtered_fastq} | \
        samtools view -b -h - | \
        samtools sort - | \
        samtools view -h -o {output.sam}
        """

rule racon_r1:
    group: "round_1"
    input:
        filtered_fastq = f"{ASSEMBLY_DIR}/{{sample}}/long_read/{{sample}}.{{gene}}.length_filtered.fastq",
        sam = f"{ASSEMBLY_DIR}/{{sample}}/assembly/racon/{{sample}}.{{gene}}.sorted.mini_1.sam",
        contigs = f"{ASSEMBLY_DIR}/{{sample}}/assembly/canu/{{gene}}/{{sample}}.{{gene}}.header_md.contigs.fasta"
    output:
        contigs = f"{ASSEMBLY_DIR}/{{sample}}/assembly/racon/{{sample}}.{{gene}}.r1.fasta"
    resources: time_min=10, mem_mb=16000, cores=4
    shell:
        """
        module load racon/1.3.1
        racon -m 8 -x -6 -g -8 -w 500 -t 4  {input.filtered_fastq} {input.sam} {input.contigs} > {output.contigs}
        sed -i '1s/.*/>{wildcards.sample}.{wildcards.gene}.r1/' {output.contigs}
        """

rule minimap2_r2_assembly:
    group: "round_2"
    input:
        contigs = f"{ASSEMBLY_DIR}/{{sample}}/assembly/racon/{{sample}}.{{gene}}.r1.fasta",
        filtered_fastq = f"{ASSEMBLY_DIR}/{{sample}}/long_read/{{sample}}.{{gene}}.length_filtered.fastq"
    output:
        sam = temp(f"{ASSEMBLY_DIR}/{{sample}}/assembly/racon/{{sample}}.{{gene}}.sorted.mini_2.sam")
    resources: time_min=60, mem_mb=16000, cores=4
    shell:
        """
        module load minimap2/2.26 samtools/1.17
        minimap2 -ax map-ont -t 2 {input.contigs} {input.filtered_fastq} | \
        samtools view -b -h - | \
        samtools sort - | \
        samtools view -h -o {output.sam}
        """

rule racon_r2:
    group: "round_2"
    input:
        filtered_fastq = f"{ASSEMBLY_DIR}/{{sample}}/long_read/{{sample}}.{{gene}}.length_filtered.fastq",
        sam = f"{ASSEMBLY_DIR}/{{sample}}/assembly/racon/{{sample}}.{{gene}}.sorted.mini_2.sam",
        contigs = f"{ASSEMBLY_DIR}/{{sample}}/assembly/racon/{{sample}}.{{gene}}.r1.fasta"
    output:
        contigs = f"{ASSEMBLY_DIR}/{{sample}}/assembly/racon/{{sample}}.{{gene}}.r2.fasta"
    resources: time_min=60, mem_mb=16000, cores=4
    shell:
        """
        module load racon/1.3.1
        racon -m 8 -x -6 -g -8 -w 500 -t 4  {input.filtered_fastq} {input.sam} {input.contigs} > {output.contigs}
        sed -i '1s/.*/>{wildcards.sample}.{wildcards.gene}.r2/' {output.contigs}
        """

rule minimap2_r3_assembly:
    group: "round_3"
    input:
        contigs = f"{ASSEMBLY_DIR}/{{sample}}/assembly/racon/{{sample}}.{{gene}}.r2.fasta",
        filtered_fastq = f"{ASSEMBLY_DIR}/{{sample}}/long_read/{{sample}}.{{gene}}.length_filtered.fastq"
    output:
        sam = temp(f"{ASSEMBLY_DIR}/{{sample}}/assembly/racon/{{sample}}.{{gene}}.sorted.mini_3.sam")
    resources: time_min=60, mem_mb=16000, cores=4
    shell:
        """
        module load minimap2/2.26 samtools/1.17
        minimap2 -ax map-ont -t 2 {input.contigs} {input.filtered_fastq} | \
        samtools view -b -h - | \
        samtools sort - | \
        samtools view -h -o {output.sam}
        """

rule racon_r3:
    group: "round_3"
    input:
        filtered_fastq = f"{ASSEMBLY_DIR}/{{sample}}/long_read/{{sample}}.{{gene}}.length_filtered.fastq",
        sam = f"{ASSEMBLY_DIR}/{{sample}}/assembly/racon/{{sample}}.{{gene}}.sorted.mini_3.sam",
        contigs = f"{ASSEMBLY_DIR}/{{sample}}/assembly/racon/{{sample}}.{{gene}}.r2.fasta"
    output:
        contigs = f"{ASSEMBLY_DIR}/{{sample}}/assembly/racon/{{sample}}.{{gene}}.r3.fasta"
    resources: time_min=60, mem_mb=16000, cores=4
    shell:
        """
        module load racon/1.3.1
        racon -m 8 -x -6 -g -8 -w 500 -t 4  {input.filtered_fastq} {input.sam} {input.contigs} > {output.contigs}
        sed -i '1s/.*/>{wildcards.sample}.{wildcards.gene}.r3/' {output.contigs}
        """

rule minimap2_r4_assembly:
    group: "round_4"
    input:
        contigs = f"{ASSEMBLY_DIR}/{{sample}}/assembly/racon/{{sample}}.{{gene}}.r3.fasta",
        filtered_fastq = f"{ASSEMBLY_DIR}/{{sample}}/long_read/{{sample}}.{{gene}}.length_filtered.fastq"
    output:
        sam = temp(f"{ASSEMBLY_DIR}/{{sample}}/assembly/racon/{{sample}}.{{gene}}.sorted.mini_4.sam")
    resources: time_min=60, mem_mb=16000, cores=4
    shell:
        """
        module load minimap2/2.26 samtools/1.17
        minimap2 -ax map-ont -t 2 {input.contigs} {input.filtered_fastq} | \
        samtools view -b -h - | \
        samtools sort - | \
        samtools view -h -o {output.sam}
        """

rule racon_r4:
    group: "round_4"
    input:
        filtered_fastq = f"{ASSEMBLY_DIR}/{{sample}}/long_read/{{sample}}.{{gene}}.length_filtered.fastq",
        sam = f"{ASSEMBLY_DIR}/{{sample}}/assembly/racon/{{sample}}.{{gene}}.sorted.mini_4.sam",
        contigs = f"{ASSEMBLY_DIR}/{{sample}}/assembly/racon/{{sample}}.{{gene}}.r3.fasta"
    output:
        contigs = f"{ASSEMBLY_DIR}/{{sample}}/assembly/racon/{{sample}}.{{gene}}.r4.fasta"
    resources: time_min=60, mem_mb=16000, cores=4
    shell:
        """
        module load racon/1.3.1
        racon -m 8 -x -6 -g -8 -w 500 -t 4  {input.filtered_fastq} {input.sam} {input.contigs} > {output.contigs}
        sed -i '1s/.*/>{wildcards.sample}.{wildcards.gene}.r4/' {output.contigs}
        """

#######################################################
##### 2. Consnesus call with Medaka ###################
#######################################################
rule medaka_CC:
    group: "rule_medaka"
    input:
        filtered_fastq = f"{ASSEMBLY_DIR}/{{sample}}/long_read/{{sample}}.{{gene}}.length_filtered.fastq",
        contigs = f"{ASSEMBLY_DIR}/{{sample}}/assembly/racon/{{sample}}.{{gene}}.r4.fasta"
    output:
        consensus_fasta = f"{ASSEMBLY_DIR}/{{sample}}/assembly/medaka/{{sample}}.{{gene}}.m.fasta"
    params: 
        output_dir = f"{ASSEMBLY_DIR}/{{sample}}/assembly/medaka/{{sample}}.{{gene}}"
    resources: time_min=20, mem_mb=16000, cores=4
    conda: "/nas/longleaf/home/farhang/envs/medaka.yaml"
    shell:
        """
        #module load medaka/1.7.2
        medaka_consensus -i {input.filtered_fastq} -d {input.contigs} -o {params.output_dir} -t 4
        cp {params.output_dir}/consensus.fasta {output.consensus_fasta}
        sed -i '1s/.*/>{wildcards.sample}.{wildcards.gene}.m/' {output.consensus_fasta}
        """
