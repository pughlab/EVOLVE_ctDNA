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

### FORMAT DATA ####################################################################################
# format data
vaf.reshaped <- reshape(
	vaf.data[,c('Patient.ID','Group.ctdna','Final')],
	direction = 'wide',
	timevar = 'Group.ctdna',
	idvar = 'Patient.ID'
	);
colnames(vaf.reshaped) <- gsub('Final\\.','',colnames(vaf.reshaped));

# calculate clearance
vaf.reshaped$delta <- apply(
	vaf.reshaped[,c('baseline','on.trial','EOT')],
	1,
	function(i) {
		i <- na.omit(i);
		if (length(i) < 2) { return(NA) } else { return( 100* ( (i[2]-i[1]) / i[1]) ) }
		}
	);

# merge mutation and clinical info
master.matrix <- merge(
	clinical,
	vaf.reshaped,
	by = 'Patient.ID'
	);

master.matrix[,c('delta.cycles','N.cycles')] <- 0;
for (i in 1:nrow(master.matrix)) {
	patient <- as.character(master.matrix[i,]$Patient.ID);
	master.matrix[i,]$delta.cycles <- sort(vaf.data[which(vaf.data$Patient == patient),]$Timepoint)[2];
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
response.p <- kruskal.test(x = master.matrix$delta, g = master.matrix$Response.cat)$p.value;
cohort.p <- kruskal.test(x = master.matrix$delta, g = master.matrix$COHORT)$p.value;
brca.p <- kruskal.test(x = master.matrix$delta, g = master.matrix$BRCA.cat)$p.value;

age.stats <- as.numeric(
	cor.test(x = master.matrix$delta, y = master.matrix$AGE, method = 'spearman', use = 'pairwise')
	[c('estimate','p.value')]);
cycle.stats <- as.numeric(
	cor.test(x = master.matrix$delta, y = master.matrix$N.cycles, method = 'spearman', use = 'pairwise')
	[c('estimate','p.value')]);

cohort.key <- list(text = list(lab = display.statistical.result(cohort.p, statistic.type = 'p', symbol = ' = ')), corner = c(0,1));
response.key <- list(text = list(lab = display.statistical.result(response.p, statistic.type = 'p', symbol = ' = ')), corner = c(0,1));
brca.key <- list(text = list(lab = display.statistical.result(brca.p, statistic.type = 'p', symbol = ' = ')), corner = c(0,1));

age.key <- list(text = list(lab = c(
	paste0('rho = ', round(age.stats[1],2)),
	display.statistical.result(age.stats[2], statistic.type = 'p', symbol = ' = ')
	)), corner = c(1,1), x = 1, y = 1);
cycle.key <- list(text = list(lab = c(
	paste0('rho = ', round(cycle.stats[1],2)),
	display.statistical.result(cycle.stats[2], statistic.type = 'p', symbol = ' = ')
	)), corner = c(1,1), x = 1, y = 1);

# make some plots
plot.objects <- list();
plot.objects[[1]] <- create.boxplot(
	delta ~ COHORT,
	master.matrix,
	add.stripplot = TRUE,
	xlab.label = 'Cohort',
	xaxis.lab = NA,
	xaxis.rot = 90,
	xaxis.fontface = 'plain',
	yaxis.fontface = 'plain',
	xlab.cex = 1,
	ylab.cex = 1,
	xaxis.tck = c(0.5,0),
	yaxis.tck = c(0.5,0),
	xaxis.cex = 0.8,
	yaxis.cex = 0.8,
	ylimits = c(-110,350),
	yat = seq(-100,300,100),
	ylab.label = expression(bold(Delta * 'ctDNA')),
	legend = list(inside = list(fun = draw.key, args = list(key = cohort.key)))
	);

plot.objects[[2]] <- create.boxplot(
	delta ~ BRCA.cat,
	master.matrix,
	add.stripplot = TRUE,
	xlab.label = 'BRCA Status',
	xaxis.lab = NA,
	xaxis.rot = 90,
	xaxis.fontface = 'plain',
	yaxis.fontface = 'plain',
	xaxis.tck = c(0.5,0),
	yaxis.tck = c(0.5,0),
	xlab.cex = 1,
	ylab.cex = 1,
	xaxis.cex = 0.8,
	yaxis.cex = 0.8,
	ylimits = c(-110,350),
	yat = seq(-100,300,100),
	legend = list(inside = list(fun = draw.key, args = list(key = brca.key))),
	ylab.label = NULL
	);

plot.objects[[3]] <- create.boxplot(
	delta ~ Best.Response,
	master.matrix,
	add.stripplot = TRUE,
	xlab.label = 'Best Response',
	xaxis.lab = NA,
	xaxis.rot = 90,
	xaxis.fontface = 'plain',
	yaxis.fontface = 'plain',
	xaxis.tck = c(0.5,0),
	yaxis.tck = c(0.5,0),
	xlab.cex = 1,
	ylab.cex = 1,
	xaxis.cex = 0.8,
	yaxis.cex = 0.8,
	ylimits = c(-110,350),
	yat = seq(-100,300,100),
	legend = list(inside = list(fun = draw.key, args = list(key = response.key))),
	ylab.label = NULL
	);

plot.objects[[4]] <- create.scatterplot(
	delta ~ AGE,
	master.matrix,
	xlab.label = 'Age',
	xaxis.lab = NA,
	xaxis.fontface = 'plain',
	yaxis.fontface = 'plain',
	xaxis.tck = c(0.5,0),
	yaxis.tck = c(0.5,0),
	xlab.cex = 1,
	ylab.cex = 1,
	xaxis.cex = 0.8,
	yaxis.cex = 0.8,
	ylimits = c(-110,350),
	yat = seq(-100,300,100),
	legend = list(inside = list(fun = draw.key, args = list(key = age.key))),
	ylab.label = NULL
	);

plot.objects[[5]] <- create.scatterplot(
	delta ~ N.cycles,
	master.matrix,
	xlab.label = 'N cycles',
	xaxis.lab = NA,
	xaxis.fontface = 'plain',
	yaxis.fontface = 'plain',
	xaxis.tck = c(0.5,0),
	yaxis.tck = c(0.5,0),
	xlab.cex = 1,
	ylab.cex = 1,
	xaxis.cex = 0.8,
	yaxis.cex = 0.8,
	ylimits = c(-110,350),
	yat = seq(-100,300,100),
	legend = list(inside = list(fun = draw.key, args = list(key = cycle.key))),
	ylab.label = NULL
	);

create.multipanelplot(
	plot.objects = plot.objects,
	layout.width = 5,
	layout.height = 1,
	plot.objects.widths = c(1.08,1,1,1,1),
	left.legend.padding = 0,
	right.legend.padding = 0,
	bottom.legend.padding = 0,
	top.legend.padding = 0,
	ylab.axis.padding = 1,
	xlab.axis.padding = 2,
	height = 4.5,
	width = 12,
	resolution = 800,
	filename = generate.filename('EVOLVE_ctDNA', 'tumourClearance_stats','png')
	);

write.table(
	master.matrix,
	file = generate.filename('EVOLVE_ctDNA','tumourClearance_data','tsv'),
	row.names = FALSE,
	col.names = TRUE,
	sep = '\t'
	);

### SAVE SESSION INFO ##############################################################################
save.session.profile(generate.filename('SummaryStats','SessionProfile','txt'));
