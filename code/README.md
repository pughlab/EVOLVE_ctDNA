# Contents
## perl
- run_agena_basecounts.pl : run basecounts to extract genotype of Agena SNPs
- run_bam_readcount.pl : run bam-readcount for known SNP positions (identified from prior whole-exome sequencing of matched patient tissues)
- run_bedtools_coverage.pl : run bedtools coverage for target intervals (CHARM + EVOLVE panel)
- run_get_insert_size.pl : run to extract fragment sizes from BAM files (deprecated; use Picard CollectInsertSize instead)

## R
- agena_fingerprinting.R : format and contrast agena SNP vafs
- assess_contamination.R : compares sample pairs with high contamination
- compare_tumour_content_estimates.R : compare all tumour content estimations
- compare_tumour_content_max_vaf.R : use maximum somatic VAF to estimate tumour content
- collect_global_fragment_sizes.R : collect and format fragment sizes (global)
- collect_mutation_fragment_sizes.R : collect and format fragment sizes (mutation-specific)
- contrast_fragment_sizes_CCNE1.R : compare fragmentation profiles between CCNE1-amp and non-amp
- custom_format_ensemble_mutations.R : generate ensemble mutation calls (no matched normals)
- extract_fragment_sizes.R : extract global fragment sizes using R biovizBase
- extract_fragment_sizes__sites.R : extract fragment sizes at specific sites using R biovizBase
- extract_triallelic_sites_mutect.R : use MuTect (v1) stats file to identify triallelic sites
- filter_mavis_to_regions.R : filter MAVIS SV calls to target regions
- find_reversions.R : compare variant calls to identify potential reversion sites
- find_triallelic_snps.R : use bam-readcounts for germline SNP to find triallelic variants
- format_and_plot_tumour_content_bamreadcount.R : create plots for tumour content and clearance
- format_mutations_for_cbio.R : format mutations calls for cbioportal (combine WXS and ctDNA)
- format_readcount_output.R
- plot_coverage.R
- plot_reversions.R
- plot_wxs_ctdna_comparison.R
- run_panelCN_mops.R : custom version from pipeline-suite; use to evaluate alternate parameters
- plot_CN_mops_data.R : plot CCNE1 copy-number findings (per-exon)
- validate_mutation_calls.R : contrast WXS and ctDNA variants
- make_validation_venn.R : make venn diagram to visualize overlap between variant calls from tissue exome and ctDNA

## bash
- run_spliceAI.sh : run spliceAI on somatic mutation sites
