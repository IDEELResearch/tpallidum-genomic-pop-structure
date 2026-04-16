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
samples = list(SampleList[0])

# Get unique genes from all input files
input_files = [f"{BED_DIR}/{sample}.bed" for sample in samples]
genes = get_genes_from_files(input_files)
#genes=["TP0620","TP0621"]
#######################################################

rule all:
    input:
        expand(f"{ASSEMBLY_DIR}/{{sample}}/assembly/canu/{{gene}}/{{sample}}.{{gene}}.finished.txt", sample=samples, gene=genes)

################################################################
##### 1. Extract first and last 200 nucleotides of ONT reads ###
################################################################

rule extract_200nc:
    group: "read_process_1"
    input:
        fastq = f"{WORK_DIR}/raw_long_read/nichols/{{sample}}.fastq"
    output:
        r1 = f"{ASSEMBLY_DIR}/{{sample}}/tmp_short_fastq/{{sample}}.R1.fastq",
        r2 = f"{ASSEMBLY_DIR}/{{sample}}/tmp_short_fastq/{{sample}}.R2.fastq"
    resources: time_min=10, mem_mb=8000, cores=2
    shell:
        """
        module load cutadapt/4.4
        cutadapt -j 2 -l 200 {input.fastq} > {output.r1}
        cutadapt -j 2 -l -200 {input.fastq} > {output.r2}
        """

################################################################
##### 2. Map reads on the masked reference #####################
################################################################

rule alignment_200nc:
    group: "read_process_2"
    input:
        r1 = f"{ASSEMBLY_DIR}/{{sample}}/tmp_short_fastq/{{sample}}.R1.fastq",
        r2 = f"{ASSEMBLY_DIR}/{{sample}}/tmp_short_fastq/{{sample}}.R2.fastq",
        masked_reference = f"{ASSEMBLY_DIR}/{{sample}}/references/CP004010.2.{{sample}}.masked.fa"
    output:
        raw_sam = temp(f"{ASSEMBLY_DIR}/{{sample}}/tmp_short_bam/{{sample}}.raw.sam"),
        unsorted_bam = temp(f"{ASSEMBLY_DIR}/{{sample}}/tmp_short_bam/{{sample}}.unsorted.bam"),
        sorted_bam = f"{ASSEMBLY_DIR}/{{sample}}/tmp_short_bam/{{sample}}.sorted.bam",
        unmapped_unsorted_bam = temp(f"{ASSEMBLY_DIR}/{{sample}}/tmp_short_bam/{{sample}}.unmapped.unsorted.bam"),
        unmapped_sorted_bam = f"{ASSEMBLY_DIR}/{{sample}}/tmp_short_bam/{{sample}}.unmapped.sorted.bam",
        bam_insert_sizes = f"{ASSEMBLY_DIR}/{{sample}}/tmp_short_bam/{{sample}}.insert-sizes.txt"
    resources: time_min=30, mem_mb=8000, cores=2
    shell:
        """
        module load minimap2/2.26 samtools/1.17
        minimap2 -ax sr -t {resources.cores} {input.masked_reference} {input.r1} {input.r2} > {output.raw_sam}
        samtools view -@ {resources.cores} -F 4 -h -b -o {output.unsorted_bam} {output.raw_sam}
        samtools sort -@ {resources.cores} -O bam -o {output.sorted_bam} {output.unsorted_bam}
        samtools index -@ {resources.cores} {output.sorted_bam}

        samtools view -@ {resources.cores} -f 4 -h -b -o {output.unmapped_unsorted_bam} {output.raw_sam}
        samtools sort -@ {resources.cores} -O bam -o {output.unmapped_sorted_bam} {output.unmapped_unsorted_bam}
        samtools index -@ {resources.cores} {output.unmapped_sorted_bam}

        samtools flagstat {output.sorted_bam}
        samtools flagstat {output.unmapped_sorted_bam}

        samtools view -f 66 {output.sorted_bam} | cut -f 9 > {output.bam_insert_sizes}
        """

################################################################
##### 3. Removing excessive soft-clips #########################
################################################################

rule soft_clip_a:
    group: "read_process_3"
    input:
        f"{ASSEMBLY_DIR}/{{sample}}/tmp_short_bam/{{sample}}.sorted.bam"
    output:
        temp(f"{ASSEMBLY_DIR}/{{sample}}/tmp_short_bam/{{sample}}.softclips.tmp.bam")
    resources: time_min=10, mem_mb=8000, cores=2
    conda: "/nas/longleaf/home/farhang/envs/ngsutils.yaml"
    shell:
        '''
        bamutils removeclipping {input} {output}
        '''

rule soft_clip_b:
    group: "read_process_3"
    input:
        f"{ASSEMBLY_DIR}/{{sample}}/tmp_short_bam/{{sample}}.softclips.tmp.bam"
    output:
        f"{ASSEMBLY_DIR}/{{sample}}/tmp_short_bam/{{sample}}.read_names_to_remove_highSoftClip.txt"
    resources: time_min=10, mem_mb=8000, cores=4
    shell:
        r"""
        module load samtools/1.17
        samtools view -@ 2 {input} | grep "ZC:f:" | \
        awk '{{ for (i = 1; i <= NF; i++) {{ if ($i ~ /ZC:f:/) {{ print $1, $i }} }} }}' | \
        sed "s/ZC:f://" | \
        awk -v MAX_SOFTCLIP="0.20" '{{ if ($2 > MAX_SOFTCLIP) {{ print $1 }} }}' | \
        cut -f 1 | sort | uniq > {output}
        """

rule soft_clip_c:
    group: "read_process_3"
    input:
        bam = f"{ASSEMBLY_DIR}/{{sample}}/tmp_short_bam/{{sample}}.sorted.bam",
        read_out = f"{ASSEMBLY_DIR}/{{sample}}/tmp_short_bam/{{sample}}.read_names_to_remove_highSoftClip.txt"
    output:
        bam_sc = f"{ASSEMBLY_DIR}/{{sample}}/tmp_short_bam/{{sample}}.highSoftClip_reads.bam",
        bam_filtered = f"{ASSEMBLY_DIR}/{{sample}}/tmp_short_bam/{{sample}}.rm_softclips.bam"
    resources: time_min=10, mem_mb=8000, cores=2
    shell:
        """
        module load picard/2.26.11
        picard FilterSamReads INPUT={input.bam} OUTPUT={output.bam_sc} READ_LIST_FILE={input.read_out} FILTER=includeReadList
        picard FilterSamReads INPUT={input.bam} OUTPUT={output.bam_filtered}  READ_LIST_FILE={input.read_out} FILTER=excludeReadList
        """

################################################################
##### 4. Removing excessive mismatches #########################
################################################################

rule remove_mismatch:
    group: "read_process_4"
    input:
        bam = f"{ASSEMBLY_DIR}/{{sample}}/tmp_short_bam/{{sample}}.rm_softclips.bam",
    output:
        bam = f"{ASSEMBLY_DIR}/{{sample}}/tmp_short_bam/{{sample}}.rm_mismatch.bam",
        failed = f"{ASSEMBLY_DIR}/{{sample}}/tmp_short_bam/{{sample}}.mismatch_reads.txt"
    resources: time_min=10, mem_mb=8000, cores=2
    conda: "/nas/longleaf/home/farhang/envs/ngsutils.yaml"
    shell:
        """
        bamutils filter {input.bam} {output.bam} -failed {output.failed} -mismatch 10
        """

################################################################
##### 5. Filter based on insert size ###########################
################################################################

rule insert_size_filter:
    group: "read_process_5"
    input:
        bam = f"{ASSEMBLY_DIR}/{{sample}}/tmp_short_bam/{{sample}}.rm_mismatch.bam"
    output:
        bam = f"{ASSEMBLY_DIR}/{{sample}}/tmp_short_bam/{{sample}}.insert_size.bam",
    params:
        depth_prefix = f"{ASSEMBLY_DIR}/{{sample}}/tmp_short_bam/{{sample}}.insert_size"
    resources: time_min=10, mem_mb=8000, cores=2
    shell:
        r"""
        module load samtools/1.17 mosdepth
        samtools view -h {input.bam} | \
        awk 'substr($0,1,1)=="@" || ($9>=1000 && $9<=6000) || ($9<=-1000 && $9>=-6000)' | \
        samtools view -b -h -o {output.bam} -
        samtools index {output.bam}
        mosdepth {params.depth_prefix} {output.bam}
        """

rule extract_reads_mapped_gene:
    group: "read_process_5"
    input:
        bam = f"{ASSEMBLY_DIR}/{{sample}}/tmp_short_bam/{{sample}}.insert_size.bam",
        bed = f"{ASSEMBLY_DIR}/{{sample}}/references/tags_masked/{{gene}}.bed"
    output:
        tmp_bam = temp(f"{ASSEMBLY_DIR}/{{sample}}/tmp_short_bam/by_gene/{{sample}}.{{gene}}.bam"),
        sorted_bam = f"{ASSEMBLY_DIR}/{{sample}}/tmp_short_bam/by_gene/{{sample}}.{{gene}}.sorted.bam"
    resources: time_min=10, mem_mb=8000, cores=2
    shell:
        """
        module load samtools/1.17
        samtools view -h -b {input.bam} -L {input.bed} > {output.tmp_bam}
        samtools sort -o {output.sorted_bam}  {output.tmp_bam}
        """

################################################################
##### 6. Extract long reads mapped to gene #####################
################################################################

rule get_reads_from_bam:
    group: "read_process_6"
    input:
        bam = f"{ASSEMBLY_DIR}/{{sample}}/tmp_short_bam/by_gene/{{sample}}.{{gene}}.sorted.bam"
    output:
        read_names = f"{ASSEMBLY_DIR}/{{sample}}/long_read/{{sample}}.{{gene}}.read_names"
    resources: time_min=10, mem_mb=4000, cores=1
    shell:
        r"""
        module load samtools/1.17
        samtools view {input.bam} | cut -f 1 | sort | uniq > {output.read_names}
        """

rule filter_long_reads:
    group: "read_process_6"
    input:
        long_reads = f"{WORK_DIR}/raw_long_read/nichols/{{sample}}.fastq",
        read_names = f"{ASSEMBLY_DIR}/{{sample}}/long_read/{{sample}}.{{gene}}.read_names"
    output:
        long_reads = f"{ASSEMBLY_DIR}/{{sample}}/long_read/{{sample}}.{{gene}}.long_reads.fastq"
    resources: time_min=10, mem_mb=4000, cores=1
    shell:
        """
        module load bbmap
        filterbyname.sh in={input.long_reads} out={output.long_reads} names={input.read_names} qin=33 overwrite=true include=t
        """

def get_min_max_from_file(file):
    min_max_values = {}
    with open(file) as f:
        for line in f:
            parts = line.strip().split('\t')
            gene = parts[0]
            min_value = parts[2]
            max_value = parts[3]
            min_max_values[gene] = (min_value, max_value)
    return min_max_values

rule cutadapt_filter:
    group: "read_process_6"
    input:
        insert_size = f"{ASSEMBLY_DIR}/{{sample}}/references/insertSize/{{sample}}.insert_size",
        fastq = f"{ASSEMBLY_DIR}/{{sample}}/long_read/{{sample}}.{{gene}}.long_reads.fastq"
    output:
        filtered_fastq = f"{ASSEMBLY_DIR}/{{sample}}/long_read/{{sample}}.{{gene}}.length_filtered.fastq",
        report = f"{ASSEMBLY_DIR}/{{sample}}/long_read/{{sample}}.{{gene}}.cutadapt_report"
    params:
        too_short = f"{ASSEMBLY_DIR}/{{sample}}/long_read/{{sample}}.{{gene}}.tooShort",
        too_long = f"{ASSEMBLY_DIR}/{{sample}}/long_read/{{sample}}.{{gene}}.tooLong"
    resources: time_min=10, mem_mb=4000, cores=1
    run:
        min_max_values = get_min_max_from_file(input.insert_size)
        min_value, max_value = min_max_values[wildcards.gene]
        shell("""
           module load cutadapt/4.4
           cutadapt -m {min_value} -M {max_value} --too-short-output {params.too_short} --too-long-output {params.too_long} -o {output.filtered_fastq} {input.fastq} > {output.report}
        """)

#######################################################
##### 7. Alignment with Canu ##########################
#######################################################
rule canu_assembly:
    group: "canu_aln"
    input:
        filtered_fastq = f"{ASSEMBLY_DIR}/{{sample}}/long_read/{{sample}}.{{gene}}.length_filtered.fastq"
    output:
        finish_txt = f"{ASSEMBLY_DIR}/{{sample}}/assembly/canu/{{gene}}/{{sample}}.{{gene}}.finished.txt"
    params:
        output_dir = lambda wildcards: f"{ASSEMBLY_DIR}/{wildcards.sample}/assembly/canu/{wildcards.gene}",
        output_prefix = lambda wildcards: f"{wildcards.sample}.{wildcards.gene}"
    resources: time_min=30, mem_mb=16000, cores=4
#    conda: "/nas/longleaf/home/farhang/envs/canu2.yaml"
    shell:
        """
        module load canu/2.0
        canu maxInputCoverage=3000 maxThreads=4 maxMemory=16g -p {params.output_prefix} -d {params.output_dir} stopOnLowCoverage=1 genomeSize=6k useGrid=false -fast -nanopore {input.filtered_fastq} || true
        
        touch {output.finish_txt}
        """

