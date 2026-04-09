# tpallidum-genomic-pop-structure
Workflows for variant calling of Treponema pallidum subsp. pallidum whole-genome Illumina data, developed for the manuscript “Treponema pallidum subsp. pallidum genetic population structure and relevance to syphilis prevention and treatment.”

## Repository overview
This repository contains pipelines for calling variants from WGS Illumina reads mapped to the Nichols (CP004010.2) and SS14 (CP004011.1) reference genomes.

The tpa_VariantCalling_Nichols.smk workflow performs read QC, host removal, alignment, post-alignment filtering, and GATK4-based variant calling for the Nichols sublineage.

## Main processing steps
* __Pre-alignment QC:__ adapter trimming (Trimmomatic), host read removal and broken read repair (bbmap), contamination estimation (seqtk + StrainSeeker).

* __Alignment:__ BWA-MEM mapping to the Nichols and SS14 reference, sorting, flagstat, flag-based filtering, duplicate marking, and mate-fix using Picard and samtools.

* __Post-alignment QC:__ removal of high-mismatch, short, soft-/hard-clipped, chimeric, low-quality (MAPQ < 40), and singleton reads using bamutils, samtools, and Picard.

* __Variant calling:__ GATK4 HaplotypeCaller variant calling, and variant filtration.

## Dependencies
The workflow assumes a Linux/HPC environment with module and conda support.
