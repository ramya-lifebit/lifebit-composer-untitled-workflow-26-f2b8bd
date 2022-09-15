nextflow.enable.dsl=2

include { gwas_vcf_regenie_1 } from './modules/gwas_vcf_regenie_1/module.nf'
include { traits_gcta_greml_1 } from './modules/traits_gcta_greml_1/module.nf'

workflow {
input1 = Channel.fromPath(params.ch_user_input_vcf)
input2 = Channel.fromPath(params.ch_king_reference_data)
input3 = Channel.fromPath(params.ch_input_pheno_transform)
input4 = Channel.fromPath(params.ch_high_ld_regions)
input5 = Channel.fromPath(params.ch_gwas_cat)
input6 = Channel.fromPath(params.ch_ld_scores)
input7 = Channel.fromPath(params.ch_pheno)
gwas_vcf_regenie_1(input1, input2, input3, input4, input5, input6, input7)
traits_gcta_greml_1(gwas_vcf_regenie_1.out.output1, gwas_vcf_regenie_1.out.output2)
}
