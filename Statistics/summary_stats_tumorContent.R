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
load('../../2022-09-06_EVOLVE_ctDNA_clinicalCovariates.RData');

vaf.data <- read.delim('2022-09-06_EVOLVE_ctDNA__estimated_tumour_content.tsv');

# get clinical data
#clinical <- read.delim('../../configs/2022-09-06_sample_info_with_batch.txt');
#clinical$Group <- factor(clinical$Group, levels = c('baseline','on.trial','EOT'));

### FORMAT DATA ####################################################################################
# format data
#vaf.data <- merge(
#	unique(clinical[,c('Patient','Sample','Timepoint','Group')]),
#	vaf.data[,c('Patient.ID','Sample','Final')],
#	by.x = c('Patient','Sample'),
#	by.y = c('Patient.ID','Sample'),
#	all = TRUE
#	);

vaf.reshaped <- reshape(
	vaf.data[,c('Patient.ID','Group.ctdna','Final')],
	direction = 'wide',
	timevar = 'Group.ctdna',
	idvar = 'Patient.ID'
	);
colnames(vaf.reshaped) <- gsub('Final\\.','',colnames(vaf.reshaped));

# merge mutation and clinical info
master.matrix <- merge(
	clinical,
	vaf.reshaped,
	by = 'Patient.ID'
	);

master.matrix$N.cycles <- 0;
for (i in 1:nrow(master.matrix)) {
	patient <- as.character(master.matrix[i,]$Patient.ID);
	master.matrix[i,]$N.cycles <- max(vaf.data[which(vaf.data$Patient == patient),]$Timepoint);
	}

master.matrix$Best.Response <- factor(
	master.matrix$Best.Response,
	levels = c('Stable Disease','Partial Response','Progressive Disease (Objective)','Inevaluable'),
	labels = c('Stable Disease','Partial Response','Progressive Disease','Inevaluable')
	);

master.matrix$COHORT <- factor(
	master.matrix$COHORT,
	levels = c('1. Platinum Sensitive', '2. Platinum Resistant', '3. Exploratory'),
	labels = c('Platinum Sensitive', 'Platinum Resistant', 'Exploratory')
	);

master.matrix$BRCA.cat <- factor(
	master.matrix$BRCA.cat,
	levels = c('WT','BRCA1','BRCA2')
	);

### STATISTICS #####################################################################################
# run simple stats
kruskal.test(x = master.matrix$baseline, g = master.matrix$COHORT);
kruskal.test(x = master.matrix$baseline, g = master.matrix$Best.Response);
kruskal.test(x = master.matrix$baseline, g = master.matrix$Age.cat);
kruskal.test(x = master.matrix$baseline, g = master.matrix$BRCA.cat);

cor.test(x = master.matrix$baseline, y = master.matrix$AGE, method = 'spearman', use = 'pairwise');
cor.test(x = master.matrix$baseline, y = master.matrix$N.cycles, method = 'spearman', use = 'pairwise');

kruskal.test(x = master.matrix$on.trial, g = master.matrix$COHORT);
kruskal.test(x = master.matrix$on.trial, g = master.matrix$Best.Response);
kruskal.test(x = master.matrix$on.trial, g = master.matrix$Age.cat);
kruskal.test(x = master.matrix$on.trial, g = master.matrix$BRCA.cat);

cor.test(x = master.matrix$on.trial, y = master.matrix$AGE, method = 'spearman', use = 'pairwise');

kruskal.test(x = master.matrix$EOT, g = master.matrix$COHORT);
kruskal.test(x = master.matrix$EOT, g = master.matrix$Best.Response);
kruskal.test(x = master.matrix$EOT, g = master.matrix$Age.cat);
kruskal.test(x = master.matrix$EOT, g = master.matrix$BRCA.cat);

cor.test(x = master.matrix$EOT, y = master.matrix$AGE, method = 'spearman', use = 'pairwise');

# make some plots
plot.objects <- list();
for (i in c('baseline','on.trial','EOT')) {

	for (j in c('COHORT','BRCA.cat','Best.Response','AGE','N.cycles')) {
	
		ylab <- if (j == 'COHORT') { gsub('\\.',' ', i) } else { NULL }
		if (i == 'EOT' & j == 'COHORT') { ylab <- 'end-of-treatment'; }

		ylimits <- if (i == 'EOT') { c(-0.02, 0.5) } else if (i == 'on.trial') { c(-0.01, 0.2)
			} else { c(-0.05, 0.8) }
		yat <- if (i == 'EOT') { seq(0, 0.5, 0.1) } else if (i == 'on.trial') { seq(0, 0.2, 0.05)
			} else { seq(0, 0.8, 0.2) }

		xaxis.labels <- if (i == 'EOT') { NA } else { rep('', 10) }

		xlab <- if (i == 'EOT') {
			if (j == 'COHORT') {
				'Cohort' } else if (j == 'BRCA.cat') {
				'BRCA Status' } else if (j == 'AGE') { 
				'Age' } else if (j == 'N.cycles') {
				'N cycles' } else { 'Best Response' } 
			} else { NULL }

		# get statistics
		if (j %in% c('COHORT','BRCA.cat','Best.Response')) {
			stats <- display.statistical.result(
				kruskal.test(x = master.matrix[,i], g = master.matrix[,j])$p.value,
				statistic.type = 'p', symbol = ' = ');
			p.key <- list(text = list(lab = stats),corner = c(0,1));
			} else {
			cor.stats <- as.numeric(cor.test(x = master.matrix[,i], y = master.matrix[,j],
				method = 'spearman', use = 'pairwise')[c('estimate','p.value')]);
			stats <- c(
				paste0('rho = ', round(cor.stats[1],2)),
				display.statistical.result(cor.stats[2], statistic.type = 'p', symbol = ' = ')
				);		
			p.key <- list(text = list(lab = stats),corner = c(1,1), x = 1, y = 1);
			}

		if (j %in% c('COHORT','BRCA.cat','Best.Response')) {
			plot.objects[[length(plot.objects)+1]] <- create.boxplot(
				get(i) ~ get(j),
				master.matrix,
				add.stripplot = TRUE,
				xlab.label = xlab,
				xaxis.lab = xaxis.labels,
				xaxis.rot = 90,
				xaxis.fontface = 'plain',
				yaxis.fontface = 'plain',
				xlab.cex = 1,
				ylab.cex = 1,
				xaxis.tck = if (i == 'EOT') { c(0.5,0) } else { c(0,0) },
				yaxis.tck = c(0.5,0),
				xaxis.cex = 0.8,
				yaxis.cex = 0.8,
				ylimits = ylimits,
				yat = yat,
				ylab.label = ylab,
				legend = list(inside = list(fun = draw.key, args = list(key = p.key)))
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
				ylab.label = ylab,
				legend = list(inside = list(fun = draw.key, args = list(key = p.key)))
				);
			}
		}
	}

create.multipanelplot(
	plot.objects = plot.objects,
	layout.width = 5,
	layout.height = 3,
	plot.objects.widths = c(1.08,1,1,1,1),
	plot.objects.heights = c(1,1,1.7),
	left.legend.padding = 0,
	right.legend.padding = 0,
	bottom.legend.padding = 0,
	top.legend.padding = 0,
	ylab.axis.padding = 1,
	xlab.axis.padding = 2,
	y.spacing = -2,
	height = 8,
	width = 12,
	resolution = 800,
	filename = generate.filename('EVOLVE_ctDNA', 'tumourContent_stats','png')
	);

# make some plots in log10 space
plot.objects <- list();
for (i in c('baseline','on.trial','EOT')) {

	for (j in c('COHORT','BRCA.cat','Best.Response','AGE','N.cycles')) {
	
		ylab <- if (j == 'COHORT') { gsub('\\.',' ', i) } else { NULL }
		if (i == 'EOT' & j == 'COHORT') { ylab <- 'end-of-treatment'; }
	
		ylimits <- if (i == 'EOT') { c(-3, log10(1)) } else if (i == 'on.trial') {
			c(-3, log10(0.4)) } else {
			c(-3, log10(10)) }
		yat <- if (i == 'EOT') { log10(c(0.001, 0.01, 0.1, 0.2, 0.4)) } else if (i == 'on.trial') {
			log10(c(0.001, 0.01, 0.1, 0.2, 0.4)) } else {
			log10(c(0.001, 0.01, 0.1, 0.2, 0.4, 0.8)) }
		yaxis.labels <- if (i == 'EOT') {
			c(0.001, 0.01, 0.1, 0.2, 0.4) } else if (i == 'on.trial') {
			c(0.001, 0.01, 0.1, 0.2, 0.4) } else {
			c(0.001, 0.01, 0.1, 0.2, 0.4, 0.8) }

		xaxis.labels <- if (i == 'EOT') { NA } else { rep('', 10) }

		xlab <- if (i == 'EOT') {
			if (j == 'COHORT') {
				'Cohort' } else if (j == 'BRCA.cat') {
				'BRCA Status' } else if (j == 'AGE') { 
				'Age' } else if (j == 'N.cycles') {
				'N cycles' } else { 'Best Response' } 
			} else { NULL }

		# get statistics
		if (j %in% c('COHORT','BRCA.cat','Best.Response')) {
			stats <- display.statistical.result(
				kruskal.test(x = master.matrix[,i], g = master.matrix[,j])$p.value,
				statistic.type = 'p', symbol = ' = ');
			p.key <- list(text = list(lab = stats),corner = c(0,1));
			} else {
			cor.stats <- as.numeric(cor.test(x = master.matrix[,i], y = master.matrix[,j],
				method = 'spearman', use = 'pairwise')[c('estimate','p.value')]);
			stats <- c(
				paste0('rho = ', round(cor.stats[1],2)),
				display.statistical.result(cor.stats[2], statistic.type = 'p', symbol = ' = ')
				);		
			p.key <- list(text = list(lab = stats),corner = c(1,1), x = 1, y = 1);
			}

		if (j %in% c('COHORT','BRCA.cat','Best.Response')) {
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
				xaxis.tck = if (i == 'EOT') { c(0.5,0) } else { c(0,0) },
				yaxis.tck = c(0.5,0),
				xaxis.cex = 0.8,
				yaxis.cex = 0.8,
				ylimits = ylimits,
				yat = yat,
				yaxis.lab = yaxis.labels,
				ylab.label = ylab,
				legend = list(inside = list(fun = draw.key, args = list(key = p.key)))
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
				xaxis.tck = if (i == 'EOT') { c(0.5,0) } else { c(0,0) },
				yaxis.tck = c(0.5,0),
				xaxis.cex = 0.8,
				yaxis.cex = 0.8,
				ylimits = ylimits,
				yat = yat,
				yaxis.lab = yaxis.labels,
				ylab.label = ylab,
				legend = list(
					inside = list(fun = draw.key, args = list(key = p.key),
					x = if (j == 'AGE') { 0.02 } else { 0.3 }))
				);
			}
		}
	}

create.multipanelplot(
	plot.objects = plot.objects,
	layout.width = 5,
	layout.height = 3,
	plot.objects.widths = c(1.08,1,1,1,1),
	plot.objects.heights = c(1,1,1.6),
	top.padding = 0,
	left.legend.padding = 0,
	right.legend.padding = 0,
	bottom.legend.padding = 0,
	top.legend.padding = 0,
	ylab.axis.padding = 1,
	xlab.axis.padding = 2,
	y.spacing = -2,
	height = 8,
	width = 12,
	resolution = 800,
	filename = generate.filename('EVOLVE_ctDNA', 'tumourContent_stats_logged','png')
	);

### SAVE SESSION INFO ##############################################################################
save.session.profile(generate.filename('SummaryStats','SessionProfile','txt'));
