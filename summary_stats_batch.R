### summary_statistics.R ##########################################################################
# Perform survival statistics on clearance + clinical features

### PREPARE SESSION ################################################################################
# import libraries
library(BoutrosLab.plotting.general);

source('/cluster/home/sprokope/git/analysis/helper_functions/session.functions.R');

working.dir <- '/cluster/projects/ovgroup/projects/OV_Superset/EVOLVE/ctDNA/pipeline_suite/Ensemble_calls/tumour_content';

setwd(working.dir);

### READ DATA ######################################################################################
# find data files
load('../../2022-06-28_EVOLVE_ctDNA_clinicalCovariates.RData');

vaf.data <- read.delim('2022-07-28_EVOLVE_ctDNA__estimated_tumour_content.tsv');

# get patient order (only if plotting)
patient.order <- as.character(read.delim('../../tumour_content/patient_order_match_timeline.txt', header = F)$V1);

# get clinical data
clinical <- read.delim('../../configs/2022-07-27_sample_info_with_batch.txt');
#clinical <- clinical[,c('Patient','Sample','Tumor.Type')];

### FORMAT DATA ####################################################################################
# format data
vaf.data <- merge(
	unique(clinical),
	vaf.data[,c('Patient.ID','Sample','Final')],
	by.x = c('Patient','Sample'),
	by.y = c('Patient.ID','Sample'),
	all = TRUE
	);

vaf.data$Patient <- factor(vaf.data$Patient, levels = rev(patient.order));

vaf.reshaped <- reshape(
	vaf.data[,c('Patient','Group','Final')],
	direction = 'wide',
	timevar = 'Group',
	idvar = 'Patient'
	);
colnames(vaf.reshaped) <- gsub('Final\\.','',colnames(vaf.reshaped));

# merge mutation and clinical info
master.matrix <- merge(
	phenodata,
	vaf.reshaped,
	by.x = 'Patient.ID',
	by.y = 'Patient'
	);

master.matrix$N.cycles <- 0;
for (i in 1:nrow(master.matrix)) {
	patient <- master.matrix[i,]$Patient.ID;
	master.matrix[i,]$N.cycles <- max(vaf.data[which(vaf.data$Patient == patient),]$Timepoint);
	}

master.matrix$Best.Response <- factor(
	master.matrix$Best.Response,
	levels = c('Stable Disease','Partial Response','Progressive Disease (Objective)','Inevaluable'),
	labels = c('Stable Disease','Partial Response','Progressive Disease','Inevaluable')
	);

### STATISTICS #####################################################################################
# run simple stats
kruskal.test(x = master.matrix$baseline, g = master.matrix$COHORT);
kruskal.test(x = master.matrix$baseline, g = master.matrix$Best.Response);
kruskal.test(x = master.matrix$baseline, g = master.matrix$Age.cat);
kruskal.test(x = master.matrix$baseline, g = master.matrix$Ethnic.cat);

cor.test(x = master.matrix$baseline, y = master.matrix$AGE, method = 'spearman', use = 'pairwise');
cor.test(x = master.matrix$baseline, y = master.matrix$N.cycles, method = 'spearman', use = 'pairwise');

kruskal.test(x = master.matrix$on.trial, g = master.matrix$COHORT);
kruskal.test(x = master.matrix$on.trial, g = master.matrix$Best.Response);
kruskal.test(x = master.matrix$on.trial, g = master.matrix$Age.cat);
kruskal.test(x = master.matrix$on.trial, g = master.matrix$Ethnic.cat);

cor.test(x = master.matrix$on.trial, y = master.matrix$AGE, method = 'spearman', use = 'pairwise');

kruskal.test(x = master.matrix$EOT, g = master.matrix$COHORT);
kruskal.test(x = master.matrix$EOT, g = master.matrix$Best.Response);
kruskal.test(x = master.matrix$EOT, g = master.matrix$Age.cat);
kruskal.test(x = master.matrix$EOT, g = master.matrix$Ethnic.cat);

cor.test(x = master.matrix$EOT, y = master.matrix$AGE, method = 'spearman', use = 'pairwise');

# make some plots
plot.objects <- list();
for (i in c('baseline','on.trial','EOT')) {

	for (j in c('COHORT','Best.Response','AGE','N.cycles')) {
	
		ylab <- if (j == 'COHORT') { gsub('\\.',' ', i) } else { NULL }

		ylimits <- if (i == 'EOT') { c(-0.02, 0.5) } else if (i == 'on.trial') { c(-0.01, 0.2)
			} else { c(-0.05, 0.8) }
		yat <- if (i == 'EOT') { seq(0, 0.5, 0.1) } else if (i == 'on.trial') { seq(0, 0.2, 0.05)
			} else { seq(0, 0.8, 0.2) }

		xaxis.labels <- if (i == 'EOT') { NA } else { rep('', 10) }

		xlab <- if (i == 'EOT') {
			if (j == 'COHORT') {
				'Cohort' } else if (j == 'AGE') { 
				'Age' } else if (j == 'N.cycles') {
				'N cycles' } else { 'Best Response' } 
			} else { NULL }

		if (j %in% c('COHORT','Best.Response')) {
			plot.objects[[length(plot.objects)+1]] <- create.boxplot(
				get(i) ~ get(j),
				master.matrix,
				xlab.label = xlab,
				xaxis.lab = xaxis.labels,
				xaxis.rot = 90,
				xaxis.fontface = 'plain',
				yaxis.fontface = 'plain',
				xlab.cex = 1,
				ylab.cex = 1,
				xaxis.tck = c(0.5,0),
				yaxis.tck = c(0.5,0),
				xaxis.cex = 0.8,
				yaxis.cex = 0.8,
				ylimits = ylimits,
				yat = yat,
				ylab.label = ylab
				);
			} else {
			plot.objects[[length(plot.objects)+1]] <- create.scatterplot(
				get(i) ~ get(j),
				master.matrix,
				xlab.label = xlab,
				xaxis.lab = xaxis.labels,
				xaxis.fontface = 'plain',
				yaxis.fontface = 'plain',
				xlab.cex = 1,
				ylab.cex = 1,
				xaxis.tck = c(0.5,0),
				yaxis.tck = c(0.5,0),
				xaxis.cex = 0.8,
				yaxis.cex = 0.8,
				ylimits = ylimits,
				yat = yat,
				ylab.label = ylab
				);
			}
		}
	}

create.multipanelplot(
	plot.objects = plot.objects,
	layout.width = 4,
	layout.height = 3,
	plot.objects.widths = c(1,1,1,1),
	plot.objects.heights = c(1,1,1.6),
	left.legend.padding = 0,
	right.legend.padding = 0,
	bottom.legend.padding = 0,
	top.legend.padding = 0,
	ylab.axis.padding = 1,
	xlab.axis.padding = 3,
	y.spacing = -2,
	height = 8,
	width = 11,
	resolution = 800,
	filename = generate.filename('EVOLVE_ctDNA', 'tumourContent_stats','png')
	);

# make some plots in log10 space
plot.objects <- list();
for (i in c('baseline','on.trial','EOT')) {

	for (j in c('COHORT','Best.Response','AGE','N.cycles')) {
	
		ylab <- if (j == 'COHORT') {
			sub('EOT','time of progression', gsub('\\.',' ', i)) } else { NULL }

		ylimits <- if (i == 'EOT') { c(-3, log10(0.6)) } else if (i == 'on.trial') {
			c(-3, log10(0.2)) } else {
			c(-3, log10(0.8)) }
		yat <- if (i == 'EOT') { log10(c(0.001, 0.01, 0.1, 0.2, 0.4)) } else if (i == 'on.trial') {
			log10(c(0.001, 0.01, 0.1, 0.2)) } else {
			log10(c(0.001, 0.01, 0.1, 0.2, 0.4, 0.8)) }
		yaxis.labels <- if (i == 'EOT') {
			c(0.001, 0.01, 0.1, 0.2, 0.4) } else if (i == 'on.trial') {
			c(0.001, 0.01, 0.1, 0.2) } else {
			c(0.001, 0.01, 0.1, 0.2, 0.4, 0.8) }

		xaxis.labels <- if (i == 'EOT') { NA } else { rep('', 10) }

		xlab <- if (i == 'EOT') {
			if (j == 'COHORT') {
				'Cohort' } else if (j == 'AGE') { 
				'Age' } else if (j == 'N.cycles') {
				'N cycles' } else { 'Best Response' } 
			} else { NULL }

		if (j %in% c('COHORT','Best.Response')) {
			plot.objects[[length(plot.objects)+1]] <- create.boxplot(
				log10(get(i)) ~ get(j),
				master.matrix,
				add.stripplot = TRUE,
				xlab.label = xlab,
				xaxis.lab = xaxis.labels,
				xaxis.rot = 90,
				xaxis.fontface = 'plain',
				yaxis.fontface = 'plain',
				xlab.cex = 1,
				ylab.cex = 1,
				xaxis.tck = c(0.5,0),
				yaxis.tck = c(0.5,0),
				xaxis.cex = 0.8,
				yaxis.cex = 0.8,
				ylimits = ylimits,
				yat = yat,
				yaxis.lab = yaxis.labels,
				ylab.label = ylab
				);
			} else {
			plot.objects[[length(plot.objects)+1]] <- create.scatterplot(
				log10(get(i)) ~ get(j),
				master.matrix,
				alpha = 0.6,
				xlab.label = xlab,
				xaxis.lab = xaxis.labels,
				xaxis.fontface = 'plain',
				yaxis.fontface = 'plain',
				xlab.cex = 1,
				ylab.cex = 1,
				xaxis.tck = c(0.5,0),
				yaxis.tck = c(0.5,0),
				xaxis.cex = 0.8,
				yaxis.cex = 0.8,
				ylimits = ylimits,
				yat = yat,
				yaxis.lab = yaxis.labels,
				ylab.label = ylab
				);
			}
		}
	}

create.multipanelplot(
	plot.objects = plot.objects,
	layout.width = 4,
	layout.height = 3,
	plot.objects.widths = c(1,1,1,1),
	plot.objects.heights = c(1,1,1.6),
	left.legend.padding = 0,
	right.legend.padding = 0,
	bottom.legend.padding = 0,
	top.legend.padding = 0,
	ylab.axis.padding = 1,
	xlab.axis.padding = 3,
	y.spacing = -2,
	height = 8,
	width = 11,
	resolution = 800,
	filename = generate.filename('EVOLVE_ctDNA', 'tumourContent_stats_logged','png')
	);

### SAVE SESSION INFO ##############################################################################
save.session.profile(generate.filename('SummaryStats','SessionProfile','txt'));
