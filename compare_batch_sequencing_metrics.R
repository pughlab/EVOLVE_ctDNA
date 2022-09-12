### compare_batches.R ##############################################################################
# Compare data between batches using replicate sample.

### PREPARE SESSION ################################################################################
# import libraries
library(BoutrosLab.plotting.general);
library(ggfortify);             # pca autoplot

source('/cluster/home/sprokope/git/analysis/helper_functions/session.functions.R');

base.dir <- '/cluster/projects/ovgroup/projects/OV_Superset/EVOLVE/ctDNA/pipeline_suite/';
output.dir <- '/cluster/projects/ovgroup/projects/OV_Superset/EVOLVE/ctDNA/pipeline_suite/batch_effects';

setwd(base.dir);

### READ DATA ######################################################################################
# read in clinical information
sample.info <- read.delim('configs/2022-07-27_sample_info_with_batch.txt', stringsAsFactors = F);

# read in metric data
cov.data <- read.delim('BAMQC/Coverage/2022-06-13_EVOLVE_ctDNA_Coverage_summary.tsv', row.names = 1);

### FORMAT DATA ####################################################################################
# format clinical information
sample.info[nrow(sample.info)+1,] <- c('EVO-400-007','EVO-400-007_ctDNA_Screening2','Screening',NA,'baseline',NA);
sample.info[79,4] <- 0;
sample.info[79,6] <- 2;

# format coverage metrics
cov.data <- cov.data[grepl('allUNIQUE',rownames(cov.data)),];
rownames(cov.data) <- gsub('_allUNIQUE','',rownames(cov.data));

# combine data
master.matrix <- merge(
	sample.info, 
	cov.data,
	by.x = 'Sample',
	by.y = 'row.names'
	);

### PCA ############################################################################################
# indicate some factors of interest
features.list <- c('Batch','Group','Patient');

# plot PCA
pca.raw <- prcomp(master.matrix$total);

for (feature in features.list) {

	png(file = paste0('raw__pca_', feature),'png'), res = 200,
		units = 'in', height = 6, width = 6);
	print(
		autoplot(pca.raw,
			data = covariate.data,
			shape = FALSE,
			colour = feature,
			main = paste0("Raw TPM:", feature))
		);
	dev.off();
	}



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
