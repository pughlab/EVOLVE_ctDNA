## Manuscript Figures
# organize and plot clinical information and sample timelines (Figure 1)
Rscript plot_sample_timeline.R


# examine estimated tumour fraction of circulating DNA over time (Figure 2)
Rscript manuscript_figure2.R


# examine fragmentation profiles as a way to evaluate tumour content (Figure 3 | Supplementary Figure 5)
Rscript examine_fragment_sizes_global.R
Rscript examine_fragment_sizes_mutations.R


# create SNV (and CNV) heatmap for both wxs and ctDNA (Figure 4)
Rscript plot_mutation_summary.R


# examine overlap in variant calls between WES of matched tumour tissue and ctDNA (Supplementary Figure 6B)
make_validation_venn.R


# plot CCNE1 copy-number results (Supplementary Figure 7)
plot_CN_mops_data.R


