# Manuscript Figures
### Main figures
1) Figure 1: organize and plot clinical information and sample timelines
<code><pre>Rscript plot_sample_timeline.R</code></pre>


2) Figure 2: examine estimated tumour fraction of circulating DNA over time
<code><pre>Rscript plot_tumour_content.R</code></pre>


3) Figure 3: examine fragmentation profiles as a way to evaluate tumour content (also Supplementary Figure 5)
<code><pre>Rscript examine_fragment_sizes_global.R
Rscript examine_fragment_sizes_mutations.R</code></pre>


4) Figure 4: create SNV (and CNV) heatmap for both wxs and ctDNA
<code><pre>Rscript plot_mutation_summary.R</code></pre>

### Supplemental Figures
Supplementary Figure 6B: examine overlap in variant calls between WES of matched tumour tissue and ctDNA
<code><pre>make_validation_venn.R</code></pre>


Supplementary Figure 7: plot CCNE1 copy-number results
<code><pre>plot_CN_mops_data.R</code></pre>


