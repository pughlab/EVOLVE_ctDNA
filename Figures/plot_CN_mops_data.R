### plot_cnmops_profile_ccne1.R ####################################################################
# Use CN calls from panelCNmops to plot CN profile for CCNE1 to identify amplifications

### PREPARE SESSION ################################################################################
# import libraries
library(BoutrosLab.plotting.general);

source('/cluster/home/sprokope/git/analysis/helper_functions/session.functions.R');

starting.dir <- '/cluster/projects/ovgroup/projects/OV_Superset/EVOLVE/ctDNA/pipeline_suite';
#working.dir <- 'panelCNmops/sp_v1';
working.dir <- '/cluster/projects/ovgroup/projects/OV_Superset/EVOLVE/ctDNA/pipeline_suite/version2/panelCNmops';

setwd(starting.dir);

### READ DATA ######################################################################################
# get clinical data
sample.info <- read.delim('configs/2022-09-06_sample_info_with_batch.txt');

# get cn.MOPS output
setwd(working.dir);

cn.data.long <- do.call(rbind, lapply(
	list.files(pattern = 'cohort_panelcn.mops_results.tsv', recursive = TRUE), read.delim)
	);

# get cn.MOPS interval annotations
anno.data <- read.delim('batch1/PanelOfNormals/formatted_countWindows.bed');
colnames(anno.data) <- c('Chromosome','Start','End','Name','Gene','Exon');

### FORMAT DATA ####################################################################################
# format CN data
cn.data.long$CN <- as.numeric(gsub('CN','',as.character(cn.data.long$CN)));

# extract CN calls for CCNE1 only
cn.data.long <- cn.data.long[which(cn.data.long$Gene == 'CCNE1'),];

# remove duplicate/overlapped region
cn.data.long <- cn.data.long[which(cn.data.long$Exon != 'CCNE1'),];

# convert lowQual calls to NA
cn.data.long[which(cn.data.long$lowQual == 'lowQual'),]$CN <- NA;

# summarize CN mops using weighted mean
weighted.cn <- aggregate(
	CN ~ Sample,
	cn.data.long,
	function(i) {
		j <- table(i);
		sum(as.numeric(names(j))*j)/sum(j)
		}
	);

weighted.cn$Call <- 0;
weighted.cn[which(weighted.cn$CN > 2.05),]$Call <- 1;

to.write <- merge(sample.info, weighted.cn);

write.table(
	to.write,
	file = generate.filename('EVOLVE_ctDNA', '_summarized_CCNE1_calls','tsv'),
	row.names = FALSE,
	col.names = TRUE,
	sep = '\t'
	);	

# format data for plotting
cn.data <- reshape(
	cn.data.long[,c('Sample','Gene','Chr','Start','End','CN')],
	direction = 'wide',
	idvar = c('Chr','Start','End','Gene'),
	timevar = 'Sample'
	);

cn.data <- merge(
	anno.data[,-4],
	cn.data,
	by.x = c('Chromosome','Start','End','Gene'),
	by.y = c('Chr','Start','End','Gene')
	);

cn.data <- cn.data[order(cn.data$Chromosome, cn.data$Start, cn.data$End),];

# organize for plotting
colnames(cn.data) <- gsub('CN-','',gsub('\\.','-',gsub('_ctDNA','',colnames(cn.data))));
rownames(cn.data) <- paste0(cn.data$Exon, '_', cn.data$Start);

# format clinical for sample ordering
sample.info$Patient <- as.character(sample.info$Patient);
sample.info$Sample <- gsub('_ctDNA','',sample.info$Sample);
sample.info <- sample.info[order(sample.info$Patient, sample.info$Timepoint),];

keep.samples <- intersect(sample.info$Sample, colnames(cn.data));

sample.info <- sample.info[which(sample.info$Sample %in% keep.samples),];
sample.order <- sample.info$Sample;

plot.cn.data <- data.frame(t(cn.data[,sample.order]));

# indicate where to put lines
patient.breaks <- get.line.breaks(rev(sample.info$Patient));

# indicate amplified samples
ccne1.amps <- data.frame(
	WES = rep(NA, nrow(plot.cn.data)),
	MOPS = rep(0, nrow(plot.cn.data))
	);
rownames(ccne1.amps) <- rownames(plot.cn.data);

wxs.with.ccne1.amps <- c('EVO-009-004','EVO-009-006','EVO-009-007','EVO-009-009','EVO-009-011','EVO-400-007','EVO-400-008');
ccne1.amps[which(rownames(ccne1.amps) %in% sample.info[which(sample.info$Group == 'baseline'),]$Sample),]$WES <- 0;
ccne1.amps[which(rownames(ccne1.amps) %in% sample.info[which(sample.info$Patient %in% wxs.with.ccne1.amps & sample.info$Group == 'baseline'),]$Sample),]$WES <- 1;

mops.amplified <- gsub('_ctDNA','',weighted.cn[which(weighted.cn$CN > 2.05),]$Sample);
ccne1.amps[mops.amplified,]$MOPS <- 1;

# create heatmap
create.heatmap(
	plot.cn.data-2,
	cluster.dimensions = 'none',
	same.as.matrix = TRUE,
	covariates = list(
		rect = list(
			col = 'black',
			fill = rev(c('white','red')[match(ccne1.amps$MOPS, c(0,1))]),
			lwd = 0
			)
		),
	covariates.grid.border = list(col = 'black', lwd = 1),
	covariates.grid.row = list(col = 'black', lwd = 1),
	covariates.grid.col = list(col = 'black', lwd = 1),
	covariates.row.lines = patient.breaks-0.5,
	yaxis.lab = NA,
	xaxis.lab = NA,
	yaxis.cex = 0.6,
	xaxis.cex = 0.6,
	xaxis.tck = 0,
	yaxis.tck = 0,
	xlab.label = expression('CCNE1 Exon Copy-Number'),
	xlab.cex = 1,
	axis.xlab.padding = 1,
	xaxis.fontface = 'plain',
	yaxis.fontface = 'plain',
	axes.lwd = 1,
	row.colour = 'black',
	grid.row = TRUE,
	force.grid.row = TRUE,
	row.lines = patient.breaks,
	print.colour.key = TRUE,
	fill.colour = 'grey80',
	at = seq(-2.5,2.5,1), # } else { seq(-2,2,0.1) },
	colour.scheme = c('blue','white','red'),
	colourkey.labels.at = seq(-2,2,1), # } else { seq(-2,2,1) },
	colourkey.labels = seq(-2,2,1), # } else { seq(-2,2,1) },
	colourkey.cex = 1,
	height = 11, #if (length(all.samples) > 12) { 8 } else { 5 },
	width = 8,
	resolution = 800,
	filename = generate.filename('EVOLVE_ctDNA', 'panel_cna_landscape','png')
	);

# add WES calls for comparison
per.exon <- create.heatmap(
	plot.cn.data-2,
	cluster.dimensions = 'none',
	same.as.matrix = TRUE,
	yaxis.lab = NA,
	xaxis.lab = NA,
	yaxis.cex = 0.6,
	xaxis.cex = 0.6,
	xaxis.tck = 0,
	yaxis.tck = 0,
	xlab.label = expression('CCNE1 Exon Copy-Number'),
	xlab.cex = 1,
	axis.xlab.padding = 1,
	xaxis.fontface = 'plain',
	yaxis.fontface = 'plain',
	axes.lwd = 1,
	row.colour = 'black',
	grid.row = TRUE,
	force.grid.row = TRUE,
	row.lines = patient.breaks,
	print.colour.key = TRUE,
	fill.colour = 'grey80',
	at = seq(-2.5,2.5,1), # } else { seq(-2,2,0.1) },
	colour.scheme = c('blue','white','red'),
	colourkey.labels.at = seq(-2,2,1), # } else { seq(-2,2,1) },
	colourkey.labels = seq(-2,2,1), # } else { seq(-2,2,1) },
	colourkey.cex = 1
	);

summarized.plot <- create.heatmap(
	ccne1.amps,
	cluster.dimensions = 'none',
	same.as.matrix = TRUE,
	yaxis.lab = NULL,
	xaxis.lab = c('WES','ctDNA panel'),
	xaxis.cex = 0.6,
	xaxis.tck = 0,
	yaxis.tck = 0,
	xlab.label = '',
	xaxis.fontface = 'plain',
	axes.lwd = 1,
	row.colour = 'black',
	grid.row = TRUE,
	grid.col = TRUE,
	force.grid.row = TRUE,
	force.grid.col = TRUE,
	row.lines = patient.breaks,
	print.colour.key = FALSE,
	colour.scheme = c('white','red')
	);

create.multipanelplot(
	plot.objects = list(per.exon, summarized.plot),
	plot.objects.heights = 1,
	plot.objects.widths = c(8,1),
	layout.height = 1,
	layout.width = 2,
	x.spacing = 2,
	left.legend.padding = 1,
	right.legend.padding = 1,
	top.legend.padding = 0,
	bottom.legend.padding = 0,
	xlab.axis.padding = 3,
	height = 11,
	width = 8,
	resolution = 800,
	filename = generate.filename('EVOLVE_ctDNA', 'panel_cna_landscape_summarized','png')
	);

### SAVE SESSION INFO ##############################################################################
save.session.profile(generate.filename('SummarizeMOPS','SessionProfile','txt'));
