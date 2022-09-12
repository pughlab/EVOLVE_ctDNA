### plot_cnmops_profile_ccne1.R ####################################################################
# Use CN calls from panelCNmops to plot CN profile for CCNE1 to identify amplifications

### FUNCTIONS ######################################################################################
# function to generate a standardized filename
generate.filename <- function(project.stem, file.core, extension, include.date = TRUE) {

	# build up the filename
	file.name <- paste(project.stem, file.core, sep = '_');
	file.name <- paste(file.name, extension, sep = '.');

	if (include.date) {
		file.name <- paste(Sys.Date(), file.name, sep = '_');
		}

	return(file.name);
	}

# function to write session profile to file
save.session.profile <- function(file.name) {

	# open the file
	sink(file = file.name, split = FALSE);

	# write memory usage to file
	cat('### MEMORY USAGE ###############################################################');
	print(proc.time());

	# write sessionInfo to file
	cat("\n### SESSION INFO ###############################################################");
	print(sessionInfo());

	# close the file
	sink();

	}

### PREPARE SESSION ################################################################################
# import libraries
library(BoutrosLab.plotting.general);

setwd('/cluster/projects/ovgroup/projects/OV_Superset/EVOLVE/ctDNA/pipeline_suite/panelCNmops');

### READ DATA ######################################################################################
# get clinical data
clinical <- read.delim('/cluster/projects/ovgroup/projects/OV_Superset/EVOLVE/ctDNA/pipeline_suite/configs/2022-05-30_sample_info.txt')[,c('Patient','Sample','Tumor.Type')];

purity <- read.delim('../tumour_content/2022-07-13_EVOLVE_ctDNA_tp53_vafs.tsv');

# get cn.MOPS output
#cn.data <- read.delim('2022-07-19_EVOLVE_ctDNA__panelCN.mops_vs_diploid__cna_matrix.tsv');
#ratio.data <- read.delim('2022-07-19_EVOLVE_ctDNA__panelCN.mops_vs_diploid__ratio_matrix.tsv');
cn.data.long <- read.delim('2022-07-28_allUNIQUE_cohort_panelcn.mops_results.tsv');

# get cn.MOPS interval annotations
anno.data <- read.delim('PanelOfNormals/formatted_countWindows.bed');
colnames(anno.data) <- c('Chromosome','Start','End','Name','Gene','Exon');

### FORMAT DATA ####################################################################################
# format purity
#purity <- unique(purity[grepl('allUNIQUE', purity$Sample.ct),c('Patient','Sample','VAF.ct')]);
#colnames(purity) <- c('Patient','Sample','Purity');

# format CN data
cn.data.long$CN <- as.numeric(gsub('CN','',as.character(cn.data.long$CN)));
cn.data.long[which(cn.data.long$lowQual == 'lowQual'),]$CN <- NA;

cn.data <- reshape(
	cn.data.long[,c('Sample','Gene','Chr','Start','End','CN')],
	direction = 'wide',
	idvar = c('Chr','Start','End','Gene'),
	timevar = 'Sample'
	);

#colnames(cn.data)[3] <- 'Name';
#colnames(ratio.data)[3] <- 'Name';

cn.data <- merge(
	anno.data[,-4],
	cn.data,
	by.x = c('Chromosome','Start','End','Gene'),
	by.y = c('Chr','Start','End','Gene')
	);

#cn.data <- cn.data[,-3];
#cn.data$Chromosome <- factor(cn.data$Chromosome, levels = paste0('chr',c(1:22,'X','Y')));
cn.data <- cn.data[order(cn.data$Chromosome, cn.data$Start, cn.data$End),];

#ratio.data <- merge(anno.data, ratio.data, all.x = TRUE);
#ratio.data <- ratio.data[,-3];
#ratio.data$Chromosome <- factor(ratio.data$Chromosome, levels = paste0('chr',c(1:22,'X','Y')));
#ratio.data <- ratio.data[order(ratio.data$Chromosome, ratio.data$Start, ratio.data$End),];

# only plot CCNE1 for now
#plot.cn.data <- cn.data[which(cn.data$Gene == 'CCNE1'),c(1:5, grep('allUNIQUE',colnames(cn.data)))];
#plot.ratio.data <- ratio.data[which(ratio.data$Gene == 'CCNE1'),colnames(plot.cn.data)];

#colnames(plot.cn.data) <- gsub('_allUNIQUE','',gsub('\\.','-',gsub('_ctDNA','',colnames(plot.cn.data))));
#colnames(plot.ratio.data) <- colnames(plot.cn.data);

plot.cn.data <- cn.data;
colnames(plot.cn.data) <- gsub('CN-','',gsub('\\.','-',gsub('_ctDNA','',colnames(plot.cn.data))));

rownames(plot.cn.data) <- paste0(plot.cn.data$Exon, '_', plot.cn.data$Start);
#rownames(plot.ratio.data) <- paste0(plot.ratio.data$Exon, '_', plot.ratio.data$Start);

# format clinical for sample ordering
clinical$Timepoint <- 0;
clinical[grepl('^C', clinical$Tumor.Type),]$Timepoint <- as.numeric(
	gsub('^C','',gsub('D1','', clinical[grepl('^C', clinical$Tumor.Type),]$Tumor.Type))
	);
clinical[which(clinical$Tumor.Type == 'EOT'),]$Timepoint <- 35;

clinical <- merge(clinical, purity, all.x = TRUE);

clinical <- unique(clinical[order(clinical$Patient, clinical$Timepoint),]);

clinical$Sample <- gsub('_ctDNA','',clinical$Sample);

if (any(!clinical$Sample %in% colnames(plot.cn.data))) {
	clinical <- clinical[which(clinical$Sample %in% colnames(plot.cn.data)),];
	}

plot.cn.data <- data.frame(t(plot.cn.data[,clinical$Sample]));
#plot.ratio.data <- data.frame(t(plot.ratio.data[,clinical$Sample]));

# remove mostly duplicated interval
plot.cn.data <- plot.cn.data[,-1];
#plot.ratio.data <- plot.ratio.data[,-1];

patient.breaks <- get.line.breaks(rev(clinical$Patient));

# get overall call (20% intervals have a gain)
overall.cn <- rev(apply(
	plot.cn.data,
	1,
	function(i) {
		gains <- length(i[which(i > 2)]);
		losses <- length(i[which(i < 2)]);
		total <- length(i[!is.na(i)]);
		if (gains/total > 0.2) { return(2); }
		else if (losses/total > 0.2) { return(0); }
		else { return(1); }
		}
	));

high.tc <- rep(0, nrow(clinical));
high.tc[which(clinical$Purity > 0.05)] <- 1;
high.tc <- rev(high.tc);

# create heatmap
create.heatmap(
	plot.cn.data-2,
	cluster.dimensions = 'none',
	same.as.matrix = TRUE,
	covariates = list(
		rect = list(
			col = 'black',
			#fill = c('blue','#CEB1FF','white','#FF9E81','red')[match(overall.cn, c(0,1,2,3,4))],
			fill = c('blue','white','red')[match(overall.cn, c(0,1,2))],
			lwd = 0
#			),
#		rect = list(
#			col = 'black',
#			fill = c('white','black')[match(high.tc, c(0,1))],
#			lwd = 0
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

create.heatmap(
	plot.ratio.data,
	cluster.dimensions = 'none',
	same.as.matrix = TRUE,
	yaxis.lab = NA,
	xaxis.lab = NA,
	yaxis.cex = 0.6,
	xaxis.cex = 0.6,
	xaxis.tck = 0,
	yaxis.tck = 0,
	xlab.label = expression('CCNE1 Exon Copy-Number (log'['2']*'ratio)'),
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
	at = seq(-2,2,0.1),
	colour.scheme = c('blue','white','red'),
	colourkey.labels.at = seq(-2,2,1),
	colourkey.labels = seq(-2,2,1),
	colourkey.cex = 1,
	height = 11,
	width = 8,
	resolution = 400,
	filename = generate.filename('EVOLVE_ctDNA', 'panel_cna_ratio_landscape','png')
	);
