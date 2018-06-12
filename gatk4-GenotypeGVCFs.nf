#!/usr/bin/env nextflow

// Copyright (C) 2018 IARC/WHO

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

params.help = null

log.info ""
log.info "-------------------------------------------------------------------------"
log.info "  gatk4-GenotypeGVCFs v1: Exact Joint Genotyping GATK4 Best Practices         "
log.info "-------------------------------------------------------------------------"
log.info "Copyright (C) IARC/WHO"
log.info "This program comes with ABSOLUTELY NO WARRANTY; for details see LICENSE"
log.info "This is free software, and you are welcome to redistribute it"
log.info "under certain conditions; see LICENSE for details."
log.info "-------------------------------------------------------------------------"
log.info ""

if (params.help)
{
    log.info "---------------------------------------------------------------------"
    log.info "  USAGE                                                 "
    log.info "---------------------------------------------------------------------"
    log.info ""
    log.info "nextflow run iarcbioinfo/gatk4-GenotypeGVCFs-nf [OPTIONS]"
    log.info ""
    log.info "Mandatory arguments:"
    log.info "--input                         VCF FILES                 All cohort gVCF files (between quotes)"
    log.info "--output_dir                    OUTPUT FOLDER             Output for VCF file"
    log.info "--cohort                        STRING                    Cohort name"
    log.info "--ref_fasta                     FASTA FILE                Reference FASTA file"
    log.info "--gatk_exec                     BIN PATH                  Full path to GATK4 executable"
    log.info "--dbsnp                         VCF FILE                  dbSNP VCF file"
    log.info "--mills                         VCF FILE                  Mills and 1000G gold standard indels VCF file"
    log.info "--axiom                         VCF FILE                  Axiom Exome Plus genotypes all populations poly VCF file"
    log.info "--hapmap                        VCF FILE                  hapmap VCF file"
    log.info "--omni                          VCF FILE                  1000G omni VCF file"
    log.info "--onekg                         VCF FILE                  1000G phase1 snps high confidence VCF file"
    log.info ""
    log.info "Optional arguments:"
    log.info "--interval_list                 INTERVAL_LIST FILE        Interval.list file For target"
    exit 1
}

//
// Parameters Init
//
params.input         = null
params.output_dir    = "."
params.cohort        = "cohort"
params.ref_fasta     = null
params.gatk_exec     = null
params.interval_list = null
params.dbsnp         = null
params.mills         = null
params.axiom         = null
params.hapmap        = null
params.omni          = null
params.onekg         = null

//
// Optional Argument Treatment
//
if (params.interval_list)
{
	interval_list_arg = "-L ${params.interval_list} "
}
else
{
	interval_list_arg = " "
}

//
// Parse Input Parameters
//
gvcf_ch = Channel
			.fromPath(params.input)

gvcf_idx_ch = Channel
			.fromPath(params.input)
			.map { file -> file+".idx" }
			
GATK                              = params.gatk_exec
ref                               = file(params.ref_fasta)
dbsnp_resource_vcf                = file(params.dbsnp)
mills_resource_vcf                = file(params.mills)
axiomPoly_resource_vcf            = file(params.axiom)
hapmap_resource_vcf               = file(params.hapmap)
omni_resource_vcf                 = file(params.omni)
one_thousand_genomes_resource_vcf = file(params.onekg)

// ExcessHet is a phred-scaled p-value. We want a cutoff of anything more extreme
// than a z-score of -4.5 which is a p-value of 3.4e-06, which phred-scaled is 54.69
excess_het_threshold = 54.69

//
// Process launching GenomicsDBImport to gather all VCFs
//
process GenomicsDBImport {

	cpus 4 
	memory '24 GB'
	time '12h'
	
	tag "${params.cohort}"

    input:
	file (gvcf) from gvcf_ch.toList()
	file (gvcf_idx) from gvcf_idx_ch.toList()

	output:
    file "*.tar" into gendb_ch
	
    script:
	"""
	${GATK} GenomicsDBImport --java-options "-Xmx4g -Xms4g" \
	${gvcf.collect { "-V $it " }.join()} \
	$interval_list_arg \
	--genomicsdb-workspace-path ${params.cohort}
	
	tar -cf ${params.cohort}.tar ${params.cohort}
	"""
}	



//
// Process launching GenotypeGVCFs on the previously created genDB
//
process GenotypeGVCFs {

	cpus 4 
	memory '24 GB'
	time '12h'
	
	tag "${params.cohort}"

	publishDir params.output_dir, mode: 'copy'

    input:
	file (workspace_tar) from gendb_ch

	output:
    set file("${params.cohort}.vcf"), file("${params.cohort}.vcf.idx") into vcf_ch

    script:
	"""
	tar -xf ${workspace_tar}
    WORKSPACE=\$( basename ${workspace_tar} .tar)

    ${GATK} --java-options "-Xmx5g -Xms5g" \
     GenotypeGVCFs \
     -R ${ref} \
     -O ${params.cohort}.vcf \
     -D ${dbsnp_resource_vcf} \
     -G StandardAnnotation \
     --only-output-calls-starting-in-intervals \
     --use-new-qual-calculator \
     -V gendb://\$WORKSPACE \
     $interval_list_arg

	"""
}	



//
// Process Hard Filtering on ExcessHet
//
process HardFilter {

	cpus 1
	memory '24 GB'
	time '12h'
	
	tag "${params.cohort}"

    input:
	set file (vcf), file (vcfidx) from vcf_ch

	output:
    set file("${params.cohort}.filtered.vcf"), file("${params.cohort}.filtered.vcf.idx") into (vcf_snv_ch, vcf_sid_ch, vcf_recal_ch)

    script:
	"""
	${GATK} --java-options "-Xmx3g -Xms3g" \
      VariantFiltration \
      --filter-expression "ExcessHet > ${excess_het_threshold}" \
      --filter-name ExcessHet \
      --exclude-filtered \
      -V ${vcf} \
      -O ${params.cohort}.filtered.vcf

	"""
}	


//
// Process SID recalibration
//
process SID_VariantRecalibrator {

	cpus 1
	memory '24 GB'
	time '12h'
	
	tag "${params.cohort}"

    input:
	set file (vcf), file (vcfidx) from vcf_sid_ch

	output:
    set file("${params.cohort}.sid.recal"),file("${params.cohort}.sid.recal.idx"),file("${params.cohort}.sid.tranches") into sid_recal_ch

    script:
	"""
    ${GATK} --java-options "-Xmx24g -Xms24g" \
      VariantRecalibrator \
      -R ${ref} \
      -V ${vcf} \
      --output ${params.cohort}.sid.recal \
      --tranches-file ${params.cohort}.sid.tranches \
      --trust-all-polymorphic \
      -an QD -an DP -an FS -an SOR -an ReadPosRankSum -an MQRankSum -an InbreedingCoeff \
      -mode INDEL \
      --max-gaussians 4 \
      -resource mills,known=false,training=true,truth=true,prior=12:${mills_resource_vcf} \
      -resource axiomPoly,known=false,training=true,truth=false,prior=10:${axiomPoly_resource_vcf} \
      -resource dbsnp,known=true,training=false,truth=false,prior=2:${dbsnp_resource_vcf}
	
	"""
}	



//
// Process SNV recalibration
//
process SNV_VariantRecalibrator {

	cpus 1
	memory '90 GB'
	time '12h'
	
	tag "${params.cohort}"

    input:
	set file (vcf), file (vcfidx) from vcf_snv_ch

	output:
    set file("${params.cohort}.snv.recal"),file("${params.cohort}.snv.recal.idx"),file("${params.cohort}.snv.tranches") into snv_recal_ch

    script:
	"""
    ${GATK} --java-options "-Xmx90g -Xms90g" \
      VariantRecalibrator \
      -R ${ref} \
      -V ${vcf} \
      --output ${params.cohort}.snv.recal \
      --tranches-file ${params.cohort}.snv.tranches \
      --trust-all-polymorphic \
      -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP -an InbreedingCoeff \
      -mode SNP \
      --max-gaussians 6 \
      -resource hapmap,known=false,training=true,truth=true,prior=15:${hapmap_resource_vcf} \
      -resource omni,known=false,training=true,truth=true,prior=12:${omni_resource_vcf} \
      -resource 1000G,known=false,training=true,truth=false,prior=10:${one_thousand_genomes_resource_vcf} \
      -resource dbsnp,known=true,training=false,truth=false,prior=7:${dbsnp_resource_vcf}
	
	"""
}	



//
// Process Apply SNV and SID recalibrations
//
process ApplyRecalibration {

	cpus 1 
	memory '7 GB'
	time '12h'
	
	tag "${params.cohort}"

	publishDir params.output_dir, mode: 'copy'

    input:
	set file (input_vcf), file (input_vcf_idx) from vcf_recal_ch
	set file (indels_recalibration), file (indels_recalibration_idx), file (indels_tranches) from sid_recal_ch
	set file (snps_recalibration), file (snps_recalibration_idx), file (snps_tranches) from snv_recal_ch

	output:
    set file("${params.cohort}.recalibrated.vcf"),file("${params.cohort}.recalibrated.vcf.idx") into vcf_final_ch

    script:
	"""
    ${GATK} --java-options "-Xmx5g -Xms5g" \
      ApplyVQSR \
      -O tmp.indel.recalibrated.vcf \
      -V ${input_vcf} \
      --recal-file ${indels_recalibration} \
      --tranches-file ${indels_tranches} \
      --truth-sensitivity-filter-level 99.0 \
      --exclude-filtered \
      --create-output-variant-index true \
      -mode INDEL

    ${GATK} --java-options "-Xmx5g -Xms5g" \
      ApplyVQSR \
      -O ${params.cohort}.recalibrated.vcf \
      -V tmp.indel.recalibrated.vcf \
      --recal-file ${snps_recalibration} \
      --tranches-file ${snps_tranches} \
      --truth-sensitivity-filter-level 99.5 \
      --exclude-filtered \
      --create-output-variant-index true \
      -mode SNP
		
	"""
}	






