# tpallidum-genomic-pop-structure
Workflows for variant calling of Treponema pallidum subsp. pallidum whole-genome Illumina data and Nanopore-based analysis of Tpr gene family, developed for the manuscript "Treponema pallidum subsp. pallidum genetic population structure and relevance to syphilis prevention and treatment."

## Repository overview
This repository contains pipelines for:
- Calling variants from WGS Illumina reads mapped to the Nichols (CP004010.2) and SS14 (CP004011.1) reference genomes (illumina folder).
- Analyzing Tpr gene family using Oxford Nanopore sequencing data, based on Pospíšilová et al. (2025) (nanopore folder).

### Illumina Variant Calling
The tpa_VariantCalling_Nichols.smk and tpa_VariantCalling_SS14.smk workflows perform read QC, host removal, alignment, post-alignment filtering, and GATK4-based variant calling for the Nichols and SS14 sublineages.

#### Main processing steps
* __Pre-alignment QC:__ adapter trimming (Trimmomatic), host read removal and broken read repair (bbmap), contamination estimation (seqtk + StrainSeeker).

* __Alignment:__ BWA-MEM mapping to the Nichols and SS14 reference, sorting, flagstat, flag-based filtering, duplicate marking, and mate-fix using Picard and samtools.

* __Post-alignment QC:__ removal of high-mismatch, short, soft-/hard-clipped, chimeric, low-quality (MAPQ < 40), and singleton reads using bamutils, samtools, and Picard.

* __Variant calling:__ GATK4 HaplotypeCaller variant calling, and variant filtration.

### Nanopore Tpr Gene Family Analysis
This analysis focuses on the Tpr gene family using Oxford Nanopore sequencing, adapted from Pospíšilová et al. (2025):
(https://journals.asm.org/doi/full/10.1128/msphere.00213-25)

#### Implementation in this repository
* __Reads processing:__ Basecalling and demultiplexing were conducted using Guppy v4.4.1 (high-accuracy model; Q20 threshold). Cutadapt was used for adapter trimming. Reads were mapped to masked Nichols (CP004010.2) and SS14 (CP004011.1) reference genomes with minimap2. Aligned reads were processed with SAMtools, NGSUtils, and Picard; reads with >20% soft clipping or >10% mismatch rates were excluded. Reads of 1–6 kb and spanning ROI ±500 nt were retained.

* __Assembly:__ Consensus sequences were assembled with Canu v1.8 and polished using Racon v1.4.13 and Medaka v0.7.1.

## Dependencies
The workflow assumes a Linux/HPC environment with module and conda support.
