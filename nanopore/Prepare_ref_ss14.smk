# Import necessary modules
import os
import pandas as pd
from prepare_ref_helpers import get_genes_from_files  # Assume you moved the function to prepare_ref_helpers.py

# Constants
WORK_DIR = "/scratch/2024-08-01_TPal_ONT"
BED_DIR = f"{WORK_DIR}/bed/ss14"
ASSEMBLY_DIR = f"{WORK_DIR}/assembly_results/ss14"

# Load sample list
sample_file = f"{WORK_DIR}/list/ss14.sampleList.txt"
SampleList = pd.read_csv(sample_file, header=None)
samples = list(SampleList[0])

# Get unique genes from all input files
input_files = [f"{BED_DIR}/{sample}.bed" for sample in samples]
genes = get_genes_from_files(input_files)

# Function to generate insert size file
def process_bed(input_file, output_file):
    awk_command = f"""
    awk '
    BEGIN {{
        i = 0;
    }}
    {{
        if (i % 2 == 0) {{
            # First line of the pair
            coordinate_1 = $2;
            name = $4;
            sub(/_.*/, "", name);
        }} else {{
            # Second line of the pair
            coordinate_2 = $3;
            size_1 = coordinate_2 - coordinate_1;
            size_2 = size_1 - 250;
            size_3 = size_1 + 250;
            printf "%s\\t%d\\t%d\\t%d\\n", name, size_1, size_2, size_3 >> "{output_file}";
        }}
        i++;
    }}
    ' {input_file}
    """
    os.system(awk_command)

if __name__ == "__main__":
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    process_bed(input_file, output_file)

# Rule definitions

rule all:
    """
    Final output after processing all samples and genes.
    """
    input:
        expand(f"{ASSEMBLY_DIR}/{{sample}}/references/CP004011.1.{{sample}}.masked.fa", sample=samples),
        expand(f"{ASSEMBLY_DIR}/{{sample}}/references/tags_masked/{{gene}}.bed", sample=samples, gene=genes),
        expand(f"{ASSEMBLY_DIR}/{{sample}}/references/insertSize/{{sample}}.insert_size", sample=samples),
        expand(f"{ASSEMBLY_DIR}/{{sample}}/references/interTag/fasta/{{gene}}.fa", sample=samples, gene=genes)

#######################################################
##### Reference files #################################
#######################################################

REFERENCE = "/References/TPallidum/Pallidum_SS14/GCA_000410555.1_ASM41055v1/GCA_000410555.1_ASM41055v1_ref/GCA_000410555.1_ASM41055v1_genomic.fa"
BED_FULL_GENOME = f"{WORK_DIR}/CP004011.1.bed"

#######################################################
##### Prepare sample and gene reference files #########
#######################################################
rule mask_reference:
    group: "all_rules_group"
    input:
        bed_full_genome = BED_FULL_GENOME,
        bed = f"{WORK_DIR}/bed/ss14/{{sample}}.bed",
        reference = REFERENCE
    output:
        f"{ASSEMBLY_DIR}/{{sample}}/references/CP004011.1.{{sample}}.masked.fa"
    resources: time_min=10, mem_mb=4000, cores=1
    shell:
        """
        module load bedtools samtools bwa picard
        bedtools subtract -a {input.bed_full_genome} -b {input.bed} > {input.bed}.subtract
        bedtools maskfasta -fi {input.reference} -bed {input.bed}.subtract -fo {output}
        samtools faidx {output}
        bwa index {output}
        picard CreateSequenceDictionary REFERENCE={output}
        """

rule split_bed:
    group: "all_rules_group"
    input:
        f"{WORK_DIR}/bed/ss14/{{sample}}.bed"
    output:
        f"{ASSEMBLY_DIR}/{{sample}}/references/tags_masked/{{gene}}.bed"
    resources: time_min=5, mem_mb=2000, cores=1
    shell:
        """
        while IFS=$'\\t' read -r chrom start end tag; do
            gene=$(echo "$tag" | cut -d'_' -f1)
            echo -e "$chrom\t$start\t$end\t$tag" >> "/work/users/f/a/farhang/scratch/2024-08-01_TPal_ONT/assembly_results/ss14/{wildcards.sample}/references/tags_masked/$gene.bed"
        done < {input}
        """

rule process_interTag:
    group: "all_rules_group"
    input:
        f"{WORK_DIR}/bed/ss14/{{sample}}.bed"
    output:
        f"{ASSEMBLY_DIR}/{{sample}}/references/interTag/{{gene}}.bed"
    params:
        inter_tag_dir = f"{ASSEMBLY_DIR}/{{sample}}/references/interTag"
    resources: time_min=5, mem_mb=4000, cores=1
    shell:
        """
        read -r chrom1 start1 end1 tag1 < <(awk 'NR==1' {input})
        read -r chrom2 start2 end2 tag2 < <(awk 'NR==2' {input})

        gene=${{tag1%_*}}
        gene="${{gene}}_interTag"

        echo -e "$chrom1\\t$start1\\t$end2\\t$gene" > "{params.inter_tag_dir}/{wildcards.gene}.bed"
        """

rule get_fasta:
    group: "all_rules_group"
    input:
        bed = f"{ASSEMBLY_DIR}/{{sample}}/references/interTag/{{gene}}.bed",
        reference = REFERENCE,
    output:
        f"{ASSEMBLY_DIR}/{{sample}}/references/interTag/fasta/{{gene}}.fa"
    resources: time_min=5, mem_mb=4000, cores=1
    shell:
        """
        module load bedtools samtools fastx_toolkit
        bedtools getfasta -name -fi {input.reference} -bed {input.bed} > {output}
        sed -i 's/\r//g' {output}
        fasta_formatter -w 0 -i {output} -o {output}.reformat.fa
        mv {output}.reformat.fa {output}
        samtools faidx {output}
        """

rule process_bed:
    group: "all_rules_group"
    input:
        f"{WORK_DIR}/bed/ss14/{{sample}}.bed"
    output:
        f"{ASSEMBLY_DIR}/{{sample}}/references/insertSize/{{sample}}.insert_size"
    resources: time_min=5, mem_mb=3000, cores=1
    shell:
        """
        python prepare_ref_process_bed.py {input} {output}
        """

