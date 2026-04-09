import pandas as pd
MyList = pd.read_csv("list/tpa.nichols.sampleList.txt", header=None)
SAMPLES = list(MyList[0])


rule all:
	input:
		expand("VCF/nichols/{sample}_HaploCaller.nichols.pass.vcf.gz", sample=SAMPLES)


#######################################################
##### T. Pallidum subsp Nichols reference #############
#######################################################

REF="/References/TPallidum/Pallidum_Nichols/GCA_000410535.2_ASM41053v2/GCA_000410535.2_ASM41053v2_ref/GCA_000410535.2_ASM41053v2_genomic.fa"

#######################################################
##### 1. Tools PATH ###################################
#######################################################

PICARD="Tools/picard-2.2.4/picard.jar"
STRAINSEEKER="Tools/strainseeker/seeker.pl"

#######################################################
##### 1. Pre-alignment QC #############################
#######################################################

##### Adapter trimming

rule adapter_trim:
	group: "trim"
	input:
		fq_1 = "raw_fastq/nichols/{sample}_R1.fastq.gz",
		fq_2 = "raw_fastq/nichols/{sample}_R2.fastq.gz"
	output:
		fq_1 = "processed_fastq/nichols/adapter_trimmed/{sample}_R1.paired_trim.fq.gz",
		fq_2 = "processed_fastq/nichols/adapter_trimmed/{sample}_R2.paired_trim.fq.gz",
		fq_1_unpaired = "processed_fastq/nichols/adapter_trimmed/{sample}_R1.unpaired.fq.gz",
		fq_2_unpaired = "processed_fastq/nichols/adapter_trimmed/{sample}_R2.unpaired.fq.gz",
		log_file = 		"processed_fastq/nichols/adapter_trimmed/{sample}_trim_log.txt"
	resources: time_min=1440, mem_mb=32000, cores=8
	params: "ILLUMINACLIP:/apps/trimmomatic/0.36/Trimmomatic-0.36/adapters/TruSeq3-PE.fa:2:30:10:8:TRUE LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"
	shell:
		r"""
		module load trimmomatic/0.36
		trimmomatic PE -threads 8 -trimlog {output.log_file} {input.fq_1} {input.fq_2} {output.fq_1} {output.fq_1_unpaired} \
		{output.fq_2} {output.fq_2_unpaired} {params}
		"""
##### a. Host genome removal
rule fq_rm_host:
	group: "trim"
	input:
		fq_1 = "processed_fastq/nichols/adapter_trimmed/{sample}_R1.paired_trim.fq.gz",
		fq_2 = "processed_fastq/nichols/adapter_trimmed/{sample}_R2.paired_trim.fq.gz"
	output:
		fq_clean = temp("processed_fastq/nichols/host_filtered/{sample}.host_filtered.fq.gz"),
		fq_cnt = "processed_fastq/nichols/host_filtered/{sample}.contamination.fq.gz"
	resources: time_min=1440, mem_mb=64000, cores=16
	params: "threads=16 -Xmx64g minid=0.95 maxindel=3 bandwidthratio=0.16 bandwidth=12 quickmatch fast minhits=2 pigz unpigz"
	shell:
		r"""
		module load bbmap/39.01
		bbmap.sh build=1 {params} \
		path=References/bbmap_masked_ref/hg19/ \
		in={input.fq_1} in2={input.fq_2} outm={output.fq_cnt} outu={output.fq_clean}
		"""

##### b. Filtering out broken reads
rule fq_rm_broken:
	group: "trim"
	input:
		"processed_fastq/nichols/host_filtered/{sample}.host_filtered.fq.gz"
	output:
		fq_1 = "processed_fastq/nichols/host_filtered/{sample}_R1.host_filtered.fq.gz",
		fq_2 = "processed_fastq/nichols/host_filtered/{sample}_R2.host_filtered.fq.gz"
	resources: time_min=360, mem_mb=32000, cores=8
	params: "-Xmx32g"
	shell:
		r"""
		module load bbmap/39.01
		repair.sh tossbrokenreads {params} in={input} out={output.fq_1} out2={output.fq_2}
		"""

##### c. Estimating bacterial contamination
rule fq_seqtk_sub:
	group: "trim"
	input:
		fq_1 = "processed_fastq/nichols/host_filtered/{sample}_R1.host_filtered.fq.gz",
		fq_2 = "processed_fastq/nichols/host_filtered/{sample}_R2.host_filtered.fq.gz",
		fq_clean = "processed_fastq/host_filtered/{sample}.host_filtered.fastq.gz"
	output:
		temp("processed_fastq/nichols/strainseeker/{sample}.sub")
	resources: time_min=180, mem_mb=32000, cores=8
	shell:
		r"""
		module load seqtk/1.3-r106
		seqtk sample -s100 {input.fq_clean} 100000 > {output}
		"""
rule strainseeker:
	group: "trim"
	input:
		"processed_fastq/nichols/strainseeker/{sample}.sub"
	output:
		"processed_fastq/nichols/strainseeker/{sample}.seeker.txt"
	resources: time_min=180, mem_mb=16000, cores=2
	shell:
		'module load r/4.1.3 perl/5.18.2 ;'
               'perl ' + STRAINSEEKER +  ' -i {input} \
                -d /Tools/strainseeker/ss_db_w32_4324 \
                -o {output}'

#######################################################
##### 2. Alignment ####################################
#######################################################

##### a. Read alignmnet with BWA
rule bwa_align:
	group: "alignment"
	input:
		fa = REF,
		strainseeker = "processed_fastq/nichols/strainseeker/{sample}.seeker.txt",
		fq_1 = "processed_fastq/nichols/host_filtered/{sample}_R1.host_filtered.fq.gz",
		fq_2 = "processed_fastq/nichols/host_filtered/{sample}_R2.host_filtered.fq.gz"
	output:
		temp("bam/nichols/aln/{sample}.nichols.raw.bam")
	resources: time_min=120, mem_mb=12000, cores=4
	params:
		rg=r"@RG\tID:{sample}\tSM:{sample}\tLB:{sample}_L"
	shell:
		r"""
		module load bwa/0.7.17 samtools/1.21
		bwa mem -R '{params.rg}' -t 2 {input.fa} {input.fq_1} {input.fq_2} | samtools view -@ 2 -b -o {output}
		"""

##### b. Fixmate correction
rule bam_fixmate:
	group: "alignment"
	input:
		"bam/nichols/aln/{sample}.nichols.raw.bam",
	output:
		temp("bam/nichols/aln/{sample}.nichols.fixmate.bam")
	resources: time_min=60, mem_mb=8000, cores=4
	shell:
		r"""
		module load samtools/1.21
		samtools fixmate -@ {resources.cores} -O bam {input}  {output}
		"""

##### c. Sort bam
rule bam_sort:
	group: "alignment"
	input:
		"bam/nichols/aln/{sample}.nichols.fixmate.bam"
	output:
		temp("bam/nichols/aln/{sample}.nichols.sorted.bam")
	resources: time_min=120, mem_mb=12000, cores=4
	shell:
		r"""
		module load samtools/1.21
		samtools sort -T {wildcards.sample} -@ {resources.cores} -O bam -o {output} {input}
		"""

rule bam_index:
	group: "alignment"
	input:
		"bam/nichols/aln/{sample}.nichols.sorted.bam"
	output:
		temp("bam/nichols/aln/{sample}.nichols.sorted.bam.bai")
	resources: time_min=30, mem_mb=8000, cores=4
	shell:
		r"""
		module load samtools/1.21
		samtools index -@ {resources.cores} {input}
		"""

##### d. Read filtering based on read flags
rule flagstat:
	group: "alignment"
	input:
		bam = "bam/nichols/aln/{sample}.nichols.sorted.bam",
		bai = "bam/nichols/aln/{sample}.nichols.sorted.bam.bai"
	output:
		temp("bam/nichols/aln/{sample}.nichols.flagstat")
	resources: time_min=60, mem_mb=8000, cores=4
	shell:
		r"""
		module load samtools/1.21
		samtools flagstat -@ {resources.cores} {input.bam} > {output}
		"""

rule filter_flags:
	group: "alignment"
	input:
		bam = "bam/nichols/aln/{sample}.nichols.sorted.bam",
		bai = "bam/nichols/aln/{sample}.nichols.sorted.bam.bai",
		flagstat = "bam/nichols/aln/{sample}.nichols.flagstat"
	output:
		temp("bam/nichols/aln/{sample}.nichols.flag_filt.bam")
	resources: time_min=60, mem_mb=8000, cores=4
	shell:
		r"""
		module load samtools/1.21
		samtools view -@ {resources.cores} -h -F 4 -F 256 -f 2 -b -o {output} {input.bam}
		""" 

rule index_flag_filter:
	group: "alignment"
	input:
		"bam/nichols/aln/{sample}.nichols.flag_filt.bam"
	output:
		temp("bam/nichols/aln/{sample}.nichols.flag_filt.bam.bai")
	resources: time_min=60, mem_mb=8000, cores=4
	shell:
		r"""
		module load samtools/1.21
		samtools index -@ {resources.cores} {input}
		"""

##### f. Marking duplicates
rule mark_dups:
	group: "alignment"
	input:
		bam = "bam/nichols/aln/{sample}.nichols.flag_filt.bam",
		bai = "bam/nichols/aln/{sample}.nichols.flag_filt.bam.bai"
	output:
		bam = temp("bam/nichols/aln/{sample}.nichols.dedup.bam"),
		metrics = "bam/nichols/aln/{sample}.dedup.metrics",
		logfile = "bam/nichols/aln/{sample}.nic.log.txt"
	resources: time_min=120, mem_mb=8000
	threads: 4
	log: "bam/nichols/aln/{sample}.nic.log.txt"
	shell:
		r"""
		set -euo pipefail
		module load java/17.0.2
		
		#Use node-local temp if available, otherwise /tmp
		TMP_ROOT="${{TMPDIR:-/tmp}}"
        java -Xmx8g -XX:ParallelGCThreads={threads} \
			-Djava.io.tmpdir="$TMP_ROOT" \
			-jar {PICARD} MarkDuplicates \
			INPUT={input.bam} \
			OUTPUT={output.bam} \
			METRICS_FILE={output.metrics} \
			TMP_DIR="$TMP_ROOT" \
			&> {log}
		"""

rule index_dedup:
	group: "alignment"
	input:
		"bam/nichols/aln/{sample}.nichols.dedup.bam"
	output:
		temp("bam/nichols/aln/{sample}.nichols.dedup.bam.bai")
	resources: time_min=60, mem_mb=8000, cores=4
	shell:
		r"""
		module load samtools/1.21
		samtools index -@ {resources.cores} {input}
		"""

##### g. Fixmate correction
rule bam_fixmate_b:
	group: "alignment"
	input:
		bam = "bam/nichols/aln/{sample}.nichols.dedup.bam",
		bai = "bam/nichols/aln/{sample}.nichols.dedup.bam.bai"
	output:
		bam = temp("bam/nichols/aln/{sample}.nichols.matefixed.bam"),
		bai = temp("bam/nichols/aln/{sample}.nichols.matefixed.bai"),
		logfile = "bam/nichols/aln/{sample}.nichols.matefixed.log.txt"
	resources: time_min=120, mem_mb=8000
	threads: 4
	log: "bam/nichols/aln/{sample}.nichols.matefixed.log.txt"
	shell:
		r"""
		set -euo pipefail
		module load java/17.0.2
		
		TMP_ROOT="${{TMPDIR:-/tmp}}" 
		java -Xmx8g -XX:ParallelGCThreads={threads} \
			-Djava.io.tmpdir="$TMP_ROOT" \
			-jar {PICARD} FixMateInformation \
			INPUT={input.bam} \
			OUTPUT={output.bam} \
			CREATE_INDEX=true \
			TMP_DIR="$TMP_ROOT" \
			&> {log}
		"""
		
		
#######################################################
##### 3. Post-alignment QC ############################
#######################################################

##### a. Removing excessive mismatches
rule remove_mismatch:
	group: "post_alignment"
	input:
		bam = "bam/nichols/aln/{sample}.nichols.matefixed.bam"
#		bai = "bam/nichols/aln/{sample}.nichols.matefixed.bai"
	output:
		bam = temp("bam/nichols/postfilter/{sample}.nichols.rm_mismatch.bam"),
		failed = "bam/nichols/postfilter/{sample}.nichols.mismatch_reads.txt"
	resources: time_min=30, mem_mb=16000, cores=2
	conda: "/nas/longleaf/home/farhang/envs/ngsutils.yaml"
	shell:
		r"""
		bamutils filter {input.bam} {output.bam} -failed {output.failed} -mismatch 5
		"""

##### b. Removing short alignments
rule short_alignment:
	group: "post_alignment"
	input:
		"bam/nichols/postfilter/{sample}.nichols.rm_mismatch.bam"
	output:
		temp("bam/nichols/postfilter/{sample}.nichols.minlen.bam")
	resources: time_min=30, mem_mb=16000, cores=2
	conda: "/nas/longleaf/home/farhang/envs/ngsutils.yaml"
	shell:
		r"""
		bamutils filter {input} {output} -minlen 35
		"""

##### c. Removing soft-clips
rule soft_clip_a:
	group: "post_alignment"
	input:
		"bam/nichols/postfilter/{sample}.nichols.minlen.bam"
	output:
		temp("bam/nichols/postfilter/{sample}.nichols.softclips.tmp.bam")
		
	resources: time_min=30, mem_mb=32000, cores=2
	conda: "/nas/longleaf/home/farhang/envs/ngsutils.yaml"
	shell:
		r"""
		bamutils removeclipping {input} {output}
		"""

rule soft_clip_b:
	group: "post_alignment"
	input:
		"bam/nichols/postfilter/{sample}.nichols.softclips.tmp.bam"
	output:
		"bam/nichols/postfilter/{sample}.read_names_to_remove_highSoftClip.txt"
	resources: time_min=60, mem_mb=32000, cores=6
	shell:
		r"""
		module load samtools/1.21
		samtools view -@ 4 {input} | grep "ZC:f:" | \
		awk '{{ for (i = 1; i <= NF; i++) {{ if ($i ~ /ZC:f:/) {{ print $1, $i }} }} }}' | \
		sed "s/ZC:f://" | \
		awk -v MAX_SOFTCLIP="0.05" '{{ if ($2 > MAX_SOFTCLIP) {{ print $1 }} }}' | \
		cut -f 1 | sort | uniq > {output}
		"""

rule soft_clip_c:
	group: "post_alignment"
	input:
		bam = "bam/nichols/postfilter/{sample}.nichols.minlen.bam",
		read_out = "bam/nichols/postfilter/{sample}.read_names_to_remove_highSoftClip.txt"
	output:
		bam_sc = temp("bam/nichols/postfilter/{sample}.nichols.softclip_filtOut_reads.bam"),
		bam_filtered = temp("bam/nichols/postfilter/{sample}.nichols.rm_softclips.bam")
	resources: time_min=60, mem_mb=32000, cores=2
	shell:
		r"""
		module load java/17.0.2
		java -Xmx32g -jar {PICARD} FilterSamReads INPUT={input.bam} OUTPUT={output.bam_sc} READ_LIST_FILE={input.read_out} FILTER=includeReadList
		java -Xmx32g -jar {PICARD} FilterSamReads INPUT={input.bam} OUTPUT={output.bam_filtered}  READ_LIST_FILE={input.read_out} FILTER=excludeReadList
		"""

##### d. Removing chimeric alignment
rule chimeric_aln:
	group: "post_alignment"
	input:
		"bam/nichols/postfilter/{sample}.nichols.rm_softclips.bam"
	output:
		bam=temp("bam/nichols/postfilter/{sample}.nichols.rm_chim.bam"),
		bam_filtOut_reads=temp("bam/nichols/postfilter/{sample}.nichols.chim_reads.bam")
	resources: time_min=60, mem_mb=32000, cores=2
	shell:
		r"""
		module load samtools/1.21
		samtools view -@ 2 -h -f 2048 -b  -o {output.bam_filtOut_reads} {input}
		samtools view -@ 2 -h -F 2048 -b  -o {output.bam} {input}
		"""

#### e. Removing hard-clips
rule hard_clip:
	group: "post_alignment"
	input:
		"bam/nichols/postfilter/{sample}.nichols.rm_chim.bam"
	output:
		bam=temp("bam/nichols/postfilter/{sample}.nichols.rm_hardclips.bam")
	resources: time_min=60, mem_mb=32000, cores=4
	shell:
		r"""
		module load samtools/1.21 
		samtools view -@ 2 -h {input} | awk '$6 !~ /H/' | samtools view -bh -o {output.bam} -
		"""

##### f. MAPQ40 filter
rule MAPQ_filt:
	group: "post_alignment"
	input:
		"bam/nichols/postfilter/{sample}.nichols.rm_hardclips.bam"
	output:
		bam=temp("bam/nichols/postfilter/{sample}.nichols.MAPQ40.bam")
	resources: time_min=60, mem_mb=32000, cores=4
	shell:
		r"""
		module load samtools/1.21
		samtools view -@ 2 -h -bq 40 {input} > {output.bam}
		"""

##### g. Remove singletone
rule rm_singletone:
	group: "post_alignment"
	input:
		"bam/nichols/postfilter/{sample}.nichols.MAPQ40.bam"
	output:
		"bam/nichols/postfilter/{sample}.nichols.final.bam"
	resources: time_min=60, mem_mb=32000, cores=4
	shell:
		r"""
		module load samtools/1.21
		samtools view -@ 2 -h -F 8 -b -o {output} {input} 
		"""

rule index_final_bam:
	group: "post_alignment"
	input:
			"bam/nichols/postfilter/{sample}.nichols.final.bam"
	output:
			"bam/nichols/postfilter/{sample}.nichols.final.bam.bai"
	resources: time_min=120, mem_mb=32000, cores=2
	shell:
			r"""
			module load samtools/1.21
			samtools index -@ {resources.cores} {input} {output}
			"""


#######################################################
##### 4. Variant calling ##############################
#######################################################

##### a. Generating individual gVCF file
rule make_gvcf:
	group: "variantCall"
	input:
		fa = REF,
		bam = "bam/nichols/postfilter/{sample}.nichols.final.bam",
		bai = "bam/nichols/postfilter/{sample}.nichols.final.bam.bai"
	output:
		"VCF/nichols/{sample}_HaploCaller.nichols.raw.g.vcf.gz"
	resources: time_min=240, mem_mb=20000, cores=4
	conda: "/envs/gatk4.yaml"
	params:
		java_opts = "-Xmx20G",
		gatk = "-ploidy 1 -stand-call-conf 30 -ERC GVCF"
	shell:
		r"""
		gatk --java-options "{params.java_opts}"  HaplotypeCaller -R {input.fa} {params.gatk} -I {input.bam} -O {output}
		"""

##### b. GATK variant calling
rule make_raw_vcf:
	group: "variantCall"
	input:
		fa = REF,
		gvcf = "VCF/nichols/{sample}_HaploCaller.nichols.raw.g.vcf.gz"
	output:
		"VCF/nichols/{sample}_HaploCaller.nichols.raw.vcf.gz"
	resources: time_min=60, mem_mb=8000, cores=1
	conda: "/envs/gatk4.yaml"
	params:
		java_opts = "-Xmx8G"
	shell:
		r"""
		gatk --java-options "{params.java_opts}"  GenotypeGVCFs -R {input.fa} -V {input.gvcf} -O {output}
		"""

rule filter_vcf:
	group: "variantCall"
	input:
		fa = REF,
		vcf = "VCF/nichols/{sample}_HaploCaller.nichols.raw.vcf.gz"
	output:
		"VCF/nichols/{sample}_HaploCaller.nichols.qual.vcf.gz"
	resources: time_min=60, mem_mb=8000, cores=1
	conda: "/envs/gatk4.yaml"
	params:
		java_opts = "-Xmx8G"
	shell:
		r"""
		gatk --java-options '{params.java_opts}' VariantFiltration \
		-R {input.fa} \
		-V {input.vcf} \
		-O {output} \
		--filter-expression 'QD < 2.0' --filter-name 'QD2' \
		--filter-expression 'MQ < 40.0' --filter-name 'MQ40' \
		--filter-expression 'FS > 60.0' --filter-name 'FS60' \
		--filter-expression 'SOR > 4.0' --filter-name 'SOR4' \
		--filter-expression 'MQRankSum <  -12.5' --filter-name 'MQRankSum-12.5' \
		--filter-expression 'ReadPosRankSum < -8.0' --filter-name 'ReadPosRankSum-8' \
		--filter-expression 'DP < 3.0' --filter-name 'LowCoverage'
		"""
rule selectVCF_pass:
	group: "variantCall"
	input:
		fa = REF,
		vcf = "VCF/nichols/{sample}_HaploCaller.nichols.qual.vcf.gz"
	output:
		"VCF/nichols/{sample}_HaploCaller.nichols.pass.vcf.gz"
	resources: time_min=30, mem_mb=8000, cores=1
	conda: "/envs/gatk4.yaml"
	params:
		java_opts = "-Xmx8G"
	shell:
		r"""
		gatk  --java-options "{params.java_opts}" SelectVariants -R {input.fa} -V {input.vcf} -O {output} -select "vc.isNotFiltered()"
		"""
