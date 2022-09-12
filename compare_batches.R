### compare_batches.R ##############################################################################
# Compare data between batches using replicate sample.

### PREPARE SESSION ################################################################################
# import libraries
library(BoutrosLab.plotting.general);

source('/cluster/home/sprokope/git/analysis/helper_functions/session.functions.R');

working.dir <- '/cluster/projects/ovgroup/projects/OV_Superset/EVOLVE/ctDNA/pipeline_suite/Ensemble_calls';

setwd(working.dir);

### READ DATA ######################################################################################
# read in mutation data
maf <- read.delim('ensemble_mutation_data.tsv');

### FORMAT DATA ####################################################################################
# trim to replicate sample
maf <- maf[grepl('EVO-400-007_ctDNA_Screening', maf$Tumor_Sample_Barcode),];

# calculate VAF
maf$t_vaf <- maf$t_alt_count / maf$t_depth;

plot.data <- reshape(
	maf[,c('Hugo_Symbol','Chromosome','Start_Position','End_Position','Variant_Type','Tumor_Sample_Barcode','t_vaf')],
	direction = 'wide',
	timevar = 'Tumor_Sample_Barcode',
	idvar = c('Hugo_Symbol','Chromosome','Start_Position','End_Position','Variant_Type')
	);

colnames(plot.data) <- gsub('_ctDNA','',gsub('t_vaf\\.','',colnames(plot.data)));

# get correlation results
cor.data <- cor(plot.data[,grepl('Screening',colnames(plot.data))], use = 'pairwise');

### PLOT DATA ######################################################################################
batches <- rep(1, ncol(cor.data));
batches[grepl('Screening2',colnames(cor.data))] <- 2;

types <- rep(1, ncol(cor.data));
types[grepl('DCS_SSCS', colnames(cor.data))] <- 2;
types[grepl('DCS$', colnames(cor.data))] <- 3;
types[grepl('allUNIQUE', colnames(cor.data))] <- 4;

# make some covariates
smp.covariates <- list(
	rect = list(
		col = 'black',
		fill = c('grey60','grey20')[match(batches, c(1,2))],
		lwd = 1
		),
	rect = list(
		col = 'black',
		fill = default.colours(4, 'pastel')[match(types, c(1,2,3,4))],
		lwd = 1
		)
	);

# make the plot legend (mutation type/consequence)
functional.legend <- legend.grob(
	legends = list(
		legend = list(
			colours = c('grey60','grey20'),
			labels = c('2021','2022'),
			title = 'Batch'
			),
		legend = list(
			colours = default.colours(4, 'pastel'),
			labels = c('SSCS','DCS_SSCS','DCS','allUNIQUE'),
			title = 'Evidence Level'
			)
		),
	title.just = 'left',
	title.fontface = 'plain',
	label.cex = 0.8,
	title.cex = 0.8,
	layout = c(1,2),
	size = 2
	);

# create heatmap
create.heatmap(
	cor.data,
	covariates = smp.covariates,
	covariates.top = smp.covariates,
	xaxis.lab = NULL,
	yaxis.lab = NULL,
	xaxis.cex = 0.5,
	yaxis.cex = 0.5,
	xaxis.fontface = 'plain',
	yaxis.fontface = 'plain',
	xaxis.tck = 0,
	yaxis.tck = 0,
	colourkey.cex = 1,
	inside.legend = list(fun = functional.legend, x = 1.2, y = 0.4),
	right.padding = 10,
	axes.lwd = 1,
	filename = generate.filename('EVOLVE_ctDNA','batch_comparison','png'),
	height = 5,
	width = 6,
	resolution = 800
	);

### SAVE SESSION INFO ##############################################################################
save.session.profile(generate.filename('CompareBatch','SessionProfile','txt'));
