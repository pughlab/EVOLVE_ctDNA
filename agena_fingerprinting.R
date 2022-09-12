### agena_fingerprinting.R #########################################################################
# Use Agena SNPs (included in the CHARM/EVOLVE panel) to evaluate contamination.

### PREAMBLE #######################################################################################
library(BoutrosLab.plotting.general);
library(plyr);

source('/cluster/home/sprokope/git/analysis/helper_functions/session.functions.R');

# indicate directory
working.dir <- '/cluster/projects/ovgroup/projects/OV_Superset/EVOLVE/ctDNA/pipeline_suite/fingerprinting';

# indicate files to get
agena.snps <- '/cluster/projects/ovgroup/projects/OV_Superset/EVOLVE/ctDNA/pipeline_suite/fingerprinting/agena_snps.txt';
clinical <- '/cluster/projects/ovgroup/projects/OV_Superset/EVOLVE/ctDNA/pipeline_suite/configs/2022-07-27_sample_info_with_batch.txt';

setwd(working.dir);

### READ DATA ######################################################################################
# read in locations of agena snps
agena <- read.delim(agena.snps);

# read in results data (VAF for each Agena SNP)
results <- join_all(lapply(
	list.files(pattern = 'frequencies.tsv', recursive = TRUE),
	function(i) {
		smp <- sub('_allUNIQUE_agenaSNP_frequencies.tsv','', basename(i));
		tmp <- read.delim(i, stringsAsFactors = FALSE);
		colnames(tmp)[5] <- smp;
		return(tmp);
		}
	));

# read in sample info
sample.info <- read.delim(clinical, stringsAsFactors = FALSE);
sample.info$Sample <- gsub('_ctDNA','',sample.info$Sample);

### FORMAT DATA ####################################################################################
# format data
results.long <- reshape(
	results,
	direction = 'long',
	varying = list(5:ncol(results)),
	v.names = 'VAF',
	timevar = 'Sample',
	times = colnames(results)[5:ncol(results)]
	);
rownames(results.long) <- 1:nrow(results.long);
results.long <- results.long[,c('chr','pos','Sample','VAF')];

# extract VAF for each base
results.long$A <- sapply(results.long$VAF, function(i) { unlist(strsplit(i,'\\|'))[1] } );
results.long$T <- sapply(results.long$VAF, function(i) { unlist(strsplit(i,'\\|'))[2] } );
results.long$G <- sapply(results.long$VAF, function(i) { unlist(strsplit(i,'\\|'))[3] } );
results.long$C <- sapply(results.long$VAF, function(i) { unlist(strsplit(i,'\\|'))[4] } );
results.long$INS <- sapply(results.long$VAF, function(i) { unlist(strsplit(i,'\\|'))[5] } );
results.long$DEL <- sapply(results.long$VAF, function(i) { unlist(strsplit(i,'\\|'))[6] } );

# convert to numeric values
for (i in c('A','T','G','C','INS','DEL')) {
	results.long[,i] <- as.numeric(results.long[,i]);
	}

results.long[which(results.long == 'NaN', arr.ind = T)] <- NA;

# reshape again for correlations
tmp <- reshape(
	results.long[,-4],
	direction = 'long',
	varying = list(4:ncol(results.long[,-4])),
	v.names = 'VAF',
	times = colnames(results.long)[5:ncol(results.long)],
	timevar = 'Base'
	);

results.wide <- reshape(
#	tmp[!grepl('Screening2',tmp$Sample),1:5],
	tmp[,1:5],
	direction = 'wide',
	idvar = c('chr','pos','Base'),
	timevar = 'Sample'
	);

results.wide <- results.wide[order(results.wide$chr, results.wide$pos, results.wide$Base),];
colnames(results.wide) <- gsub('VAF\\.','', colnames(results.wide));
rownames(results.wide) <- 1:nrow(results.wide);

cor.results <- cor(results.wide[,4:ncol(results.wide)], use = 'pairwise');
rownames(cor.results) <- gsub('_ctDNA','', rownames(cor.results));
colnames(cor.results) <- gsub('_ctDNA','', colnames(cor.results));

# add clinical for missing sample(s)
missing.smps <- setdiff(colnames(cor.results), sample.info$Sample);
for (smp in missing.smps) {
	sample.info[nrow(sample.info)+1,] <- c(substr(smp,0,11), smp, unlist(strsplit(smp,'_'))[2], NA, NA, NA);
	}

#sample.info[which(sample.info$Sample == 'EVO-400-007_Screening2'),4:6] <- c(0,'baseline',2);

sample.info$Sample <- factor(sample.info$Sample, levels = colnames(cor.results));
sample.info <- sample.info[order(sample.info$Sample),];

# organize covariates
batch.groups <- rep(0,nrow(cor.results));
batch.groups[which(rownames(cor.results) %in% sample.info[which(sample.info$Batch == 2),]$Sample)] <- 1;

time.groups <- rep(0,nrow(cor.results));
time.groups[which(rownames(cor.results) %in% sample.info[which(sample.info$Group == 'on.trial'),]$Sample)] <- 1;
time.groups[which(rownames(cor.results) %in% sample.info[which(sample.info$Group == 'EOT'),]$Sample)] <- 2;

smp.groups <- rep(0,nrow(cor.results));
#smp.groups[which(rownames(cor.results) == 'EVO-009-009_C13D1')] <- 1;
smp.groups[grepl('EVO-009-009', rownames(cor.results))] <- 1;
smp.groups[grepl('EVO-009-020', rownames(cor.results))] <- 2;
#smp.groups[which(rownames(cor.results) == 'EVO-009-010_C10D1')] <- 3;
smp.groups[grepl('EVO-009-010', rownames(cor.results))] <- 3;
smp.groups[grepl('EVO-009-016', rownames(cor.results))] <- 4;
smp.groups[grepl('EVO-400-007_Screening', rownames(cor.results))] <- 5;

covariate.colours <- list(
	rect = list(
		col = 'transparent',
		fill = c('grey60','grey20')[match(batch.groups, c(0,1))],
		lwd = 1
		),
	rect = list(
		col = 'transparent',
		fill = c('skyblue','seagreen','pink')[match(time.groups, c(0,1,2))],
		lwd = 1
		),
	rect = list(
		col = 'transparent',
		fill = c('grey80','darkgreen','green','darkorchid4','darkorchid1','dodgerblue2')[match(smp.groups, c(0,1:5))],
		lwd = 1
		)
	);

# make the plot legend (mutation type/consequence)
functional.legend <- legend.grob(
	legends = list(
		legend = list(
			colours = c('darkgreen','green','darkorchid4','darkorchid1','dodgerblue2','grey80'),
			labels = c('EVO-009-009','EVO-009-020','EVO-009-010','EVO-009-016','EVO-400-007_Screening','other'),
			title = 'Sample sets'
			),
		legend = list(
			colours = c('skyblue','seagreen','pink'),
			labels = c('baseline','on-trial','time of progression'),
			title = 'Timepoint'
			),
		legend = list(
			colours = c('grey60','grey20'),
			labels = c('2021','2022'),
			title = 'Batch'
			)
		),
	title.just = 'left',
	title.fontface = 'plain',
	label.cex = 0.8,
	title.cex = 0.8,
	layout = c(1,3),
	x = 0.91, y = 0.75,
	size = 1.5
	);

my.heatmap <- create.heatmap(
	cor.results,
	yaxis.lab = NA,
	xaxis.lab = NA,
	yaxis.cex = 0.5,
	xaxis.cex = 0.5,
	yaxis.fontface = 'plain',
	xaxis.fontface = 'plain',
	axes.lwd = 1,
	axis.xlab.padding = 1,
	bottom.padding = 2,
	right.padding = 28,
	colour.scheme = c('blue','white','red'),
	at = seq(0,1,0.001),
	colour.centering.value = 0.5,
	colourkey.cex = 1,
	xaxis.covariates = covariate.colours,
	xaxis.covariates.y = -0.005
	);

write.plot(
	my.heatmap,
	additional.trellis.objects = list(grid.draw(functional.legend)),
	additional.trellis.locations = list(xleft = 1, xright = 1, ytop = 1, ybottom = 0.5),
	filename = generate.filename('EVOLVE_ctDNA','agenaSNP_correlation','png'),
	height = 10,
	width = 12,
	resolution = 800
	);

# create 'minimal' version
cor.results <- cor.results[!grepl("Screening2", rownames(cor.results)),!grepl("Screening2", colnames(cor.results))];

smp.groups <- rep(0,nrow(cor.results));
smp.groups[grepl('EVO-009-009', rownames(cor.results))] <- 1;
smp.groups[grepl('EVO-009-020', rownames(cor.results))] <- 2;
smp.groups[grepl('EVO-009-010', rownames(cor.results))] <- 3;
smp.groups[grepl('EVO-009-016', rownames(cor.results))] <- 4;

covariate.colours <- list(
	rect = list(
		col = 'transparent',
		fill = c('grey80','darkgreen','green','darkorchid4','darkorchid1')[match(smp.groups, c(0,1:4))],
		lwd = 1
		)
	);

create.heatmap(
	cor.results,
	yaxis.lab = NA,
	xaxis.lab = NA,
	yaxis.cex = 0.6,
	xaxis.cex = 0.6,
	yaxis.fontface = 'plain',
	xaxis.fontface = 'plain',
	axes.lwd = 1,
	axis.xlab.padding = 1,
	bottom.padding = 2,
	left.padding = 2,
	colour.scheme = c('blue','white','red'),
	at = seq(0,1,0.001),
	colour.centering.value = 0.5,
	colourkey.cex = 1,
	xaxis.covariates = covariate.colours,
	xaxis.covariates.y = 0.002,
	yaxis.covariates = covariate.colours,
	yaxis.covariates.x = -0.035,
	filename = generate.filename('EVOLVE_ctDNA','agenaSNP_correlation_minimal','png'),
	height = 12,
	width = 12,
	resolution = 800
	);

### SAVE SESSION INFO ##############################################################################
save.session.profile(generate.filename('CompareAgena','SessionProfile','txt'));
