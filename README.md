# gatk4-GenotypeGVCFs-nf
Joint calling of gVCF, following GATK4 [Best Practices](https://software.broadinstitute.org/gatk/best-practices/workflow?id=11145).

<img src="https://us.v-cdn.net/5019796/uploads/editor/mz/tzm69d8e2spl.png" width="800" />

## Description

Whole cohort variant calling (joint genotyping).

## Dependencies 

1. This pipeline is based on [nextflow](https://www.nextflow.io). As we have several nextflow pipelines, we have centralized the common information in the [IARC-nf](https://github.com/IARCbioinfo/IARC-nf) repository. Please read it carefully as it contains essential information for the installation, basic usage and configuration of nextflow and our pipelines.
2. [GATK4 executables](https://software.broadinstitute.org/gatk/download/)
3. References (genome in fasta, dbSNP vcf, 1000 Genomes vcf, Mills and 1000 Genomes Gold Standard vcf, Axiom Exome Plus genotypes all populations poly VCF file, hapmap VCF file, 1000G omni VCF file), available in [GATK Bundle](https://software.broadinstitute.org/gatk/download/bundle).

## Input

- `--input` : your cohort gVCF files (between quotes) `--input "test_*.g.vcf"`)
- `--output_dir` : the folder that will contain your aligned, recalibrated, analysis-ready BAM file(s).
- `--cohort` : name of your cohort (default = cohort). 
- `--ref_fasta` : your reference in FASTA. 
- `--dbsnp` : dbSNP VCF file. 
- `--onekg` : 1000 Genomes High Confidence SNV VCF file. 
- `--mills` : Mills and 1000 Genomes Gold Standard SID VCF file. 
- `--axiom` : Axiom Exome Plus genotypes all populations poly VCF file. 
- `--hapmap` : hapmap VCF file. 
- `--omni` :1000G omni VCF file. 
- `--gatk_exec` : the full path to your GATK4 binary file.

A nextflow.config is also included, please modify it for suitability outside our pre-configured clusters ([see Nexflow configuration](https://www.nextflow.io/docs/latest/config.html#configuration-file)).

## Usage for Cobalt cluster
```
nextflow run iarcbioinfo/gatk4-GenotypeGVCFs.nf -profile cobalt --input "/data/test_*.g.vcf" --output_dir /data/myJointCall --cohort test --ref_fasta /ref/Homo_sapiens_assembly38.fasta --gatk_exec /bin/gatk-4.0.5.0/gatk --dbsnp /ref/dbsnp_146.hg38.vcf.gz --onekg /ref/1000G_phase1.snps.high_confidence.hg38.vcf.gz --mills Mills_and_1000G_gold_standard.indels.hg38.vcf.gz --axiom /ref/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz --hapmap /ref/hapmap_3.3.hg38.vcf.gz --omni 1000G_omni2.5.hg38.vcf.gz
```


