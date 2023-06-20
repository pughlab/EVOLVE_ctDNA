## Contents
# perl
run_agena_basecounts.pl # run basecounts to extract genotype of Agena SNPs
run_bam_readcount.pl # run bam-readcount for known SNP positions (identified from prior whole-exome sequencing of matched patient tissues)
run_bedtools_coverage.pl # run bedtools coverage for target intervals (CHARM + EVOLVE panel)
run_get_insert_size.pl # run to extract fragment sizes from BAM files (deprecated; use Picard CollectInsertSize instead)

# R
agena_fingerprinting.R # format and contrast agena SNP vafs
assess_contamination.R # compares sample pairs with high contamination
compare_tumour_content_estimates.R # compare all tumour content estimations
compare_tumour_content_max_vaf.R # use maximum somatic VAF to estimate tumour content
contrast_fragment_sizes_CCNE1.R # compare fragmentation profiles between CCNE1-amp and non-amp
custom_format_ensemble_mutations.R # generate ensemble mutation calls (no matched normals)
extract_fragment_sizes.R # extract global fragment sizes using R biovizBase
extract_fragment_sizes__sites.R # extract fragment sizes at specific sites using R biovizBase
extract_triallelic_sites_mutect.R # use MuTect (v1) stats file to identify triallelic sites
filter_mavis_to_regions.R # filter MAVIS SV calls to target regions
find_reversions.R # compare variant calls to identify potential reversion sites
find_triallelic_snps.R # use bam-readcounts for germline SNP to find triallelic variants
format_and_plot_tumour_content_bamreadcount.R # create plots for tumour content and clearance
format_mutations_for_cbio.R # format mutations calls for cbioportal (combine WXS and ctDNA)
format_readcount_output.R
plot_coverage.R
plot_reversions.R
plot_wxs_ctdna_comparison.R
run_panelCN_mops.R # custom version from pipeline-suite; use to evaluate alternate parameters
validate_mutation_calls.R # contrast WXS and ctDNA variants

# bash
run_spliceAI.sh # run spliceAI on somatic mutation sites


## assess contamination
perl run_agena_basecounts.pl -d /path/to/project/configs/evolve_ctDNA_raw_bams.yaml -t /path/to/project/configs/dna_pipeline_config.yaml -p /path/to/project/extdata/agena_snps_2col.txt -o /path/to/project/fingerprinting --no-wait -c slurm

Rscript agena_fingerprinting.R

## assess tumour content
perl run_bam_readcount.pl -d /path/to/project/configs/evolve_ctDNA_combined_bams.yaml -t /path/to/project/dna_pipeline_config.yaml -o /path/to/project/tumour_content -p /path/to/project/extdata/EVOLVE_panel_covered_mutations.bed --no-wait -c slurm

Rscript format_and_plot_tumour_content_bamreadcount.R

# compare above tumour content estimates with maximum VAF (run after ensemble calls)
# this will fill in most NAs
# also creates nice plots
Rscript compare_tumour_content_max_vaf.R

## examine fragmentation profiles as a way to evaluate tumour content
# check global (all reads) and reads covering high-confidence mutations
perl run_get_insert_size.pl -t /path/to/project/configs/dna_pipeline_config.yaml -d /path/to/project/configs/evolve_ctDNA_combined_bams.yaml -o /path/to/project/fragment_size -c slurm --dry-run

# compare above tumour content estimates (TP53, max VAF, short fragments)
Rscript compare_tumour_content_estimates.R

## generate ensemble mutation calls
module load R

Rscript custom_format_ensemble_mutations.R -p EVOLVE_ctDNA -o /path/to/project/Ensemble_calls --mutect mutect.tsv --mutect2 mutect2.tsv --varscan varscan.tsv --vardict vardict.tsv --strelka strelka.tsv --pindel pindel.tsv

# get sensitivity metrics for SNV/INDEL tools + ENSEMBLE (vs exome); this will also produce a VennDiagram
Rscript validate_mutation_calls.R
 
# format mutation calls for cbioportal (merge with exome data [validated] and novel pathogenic [discovery]
Rscript format_mutations_for_cbio.R

## find reversions
# get readcounts for germline positions (from WXS)
perl run_bam_readcount.pl -d /path/to/project/configs/evolve_ctDNA_combined_bams.yaml -t /path/to/project/configs/dna_pipeline_config.yaml -o /path/to/project/germline_snps -p /path/to/project/extdata/germline_snp_positions.tsv --no-wait -c slurm

# find triallelic variants
Rscript find_triallelic_snps.R -p EVOLVE_ctDNA -d /path/to/project/germline_snps -t /path/to/project/HaplotypeCaller/CPSR/2022-06-15_EVOLVE_ctDNA_mutations_for_cbioportal.tsv

# find potential splice mutations
run_spliceAI.sh

# run find_reversions.R
# please install findReversions package (git@github.com:pughlab/findReversions.git) prior to running
module load R/3.6.1

Rscript find_reversions.R -d /path/to/project/Reversions -p EVOLVE_ctDNA -q /path/to/project/extdata/EVOLVE_exome_mutations.txt -m /path/to/project/HaplotypeCaller/cohort/VCF2MAF/2022-07-15_EVOLVE_ctDNA_mutations_for_cbioportal.tsv -v /path/to/project/Mavis/filtered_mavis_output.txt -t /path/to/project/extdata/CHARM-MMR_plus_EVOLVE_hg38.bed -x /path/to/project/extdata/known_exome_reversions.maf -s /path/to/project/SpliceAI/2022-07-14_EVOLVE_ctDNA_ensemble_mutations_SPLICEAI.txt -a /path/to/project/germline_snps2022-07-14_EVOLVE_ctDNA_triallelic_SNPs.tsv
