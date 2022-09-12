### summary_statistics.R ##########################################################################
# Perform survival statistics on clearance + clinical features

### PREPARE SESSION ################################################################################
# import libraries
library(BoutrosLab.plotting.general);

source('/cluster/home/sprokope/git/analysis/helper_functions/session.functions.R');

### READ DATA ######################################################################################
# find data files
load('../2022-06-28_EVOLVE_ctDNA_clinicalCovariates.RData');

clearance.data <- read.delim('2022-07-13_EVOLVE_ctDNA_tumour_clearance.tsv');
vaf.data <- read.delim('2022-07-13_EVOLVE_ctDNA_tp53_vafs.tsv');

# get patient order (only if plotting)
patient.order <- as.character(read.delim('patient_order_match_timeline.txt', header = F)$V1);

### FORMAT DATA ####################################################################################
# format data
clearance.data$Patient <- factor(clearance.data$Patient, levels = rev(patient.order));

vaf.data$Group <- 'eot';
vaf.data[which(vaf.data$Timepoint == 2),]$Group <- 'C2';
vaf.data[which(vaf.data$Timepoint == 1),]$Group <- 'baseline';
vaf.data[which(vaf.data$Timepoint == 0),]$Group <- 'baseline';
vaf.data$Group <- factor(vaf.data$Group, levels = c('baseline','C2','eot'), labels = c('baseline','on.trial','EOT'));

vaf.reshaped <- reshape(
	unique(vaf.data[grepl('allUNIQUE', vaf.data$Sample.ct),c('Patient','Group','VAF.ct')]),
	direction = 'wide',
	timevar = 'Group',
	idvar = 'Patient'
	);
colnames(vaf.reshaped) <- gsub('VAF\\.ct\\.','',colnames(vaf.reshaped));

# merge mutation and clinical info
master.matrix <- merge(
	phenodata,
	clearance.data[,c('Patient','delta')],
	by.x = 'Patient.ID',
	by.y = 'Patient'
	);

master.matrix$N.cycles <- 0;
for (i in 1:nrow(master.matrix)) {
	patient <- as.character(master.matrix[i,]$Patient.ID);
	master.matrix[i,]$N.cycles <- max(vaf.data[which(vaf.data$Patient == patient),]$Timepoint);
	}

### STATISTICS #####################################################################################
# run simple stats
kruskal.test(x = master.matrix$delta, g = master.matrix$Response.cat);
kruskal.test(x = master.matrix$delta, g = master.matrix$COHORT);

cor.test(x = master.matrix$delta, y = master.matrix$AGE, method = 'spearman', use = 'pairwise');
cor.test(x = master.matrix$delta, y = master.matrix$N.cycles, method = 'spearman', use = 'pairwise');

master.matrix$Response.cat <- factor(master.matrix$Response.cat, levels = c('SD','PR','PD','IE'));

# make some plots
plot.objects <- list();
plot.objects[[length(plot.objects)+1]] <- create.boxplot(
	delta ~ COHORT,
	master.matrix,
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
	ylab.label = expression(Delta * 'ctDNA')
	);

plot.objects[[length(plot.objects)+1]] <- create.boxplot(
	delta ~ Response.cat,
	master.matrix,
	xlab.label = 'Best Response',
	xaxis.lab = NA,
	xaxis.fontface = 'plain',
	yaxis.fontface = 'plain',
	xaxis.tck = c(0.5,0),
	yaxis.tck = c(0.5,0),
	xlab.cex = 1,
	ylab.cex = 1,
	xaxis.cex = 0.8,
	yaxis.cex = 0.8,
	ylab.label = NULL
	);

plot.objects[[length(plot.objects)+1]] <- create.scatterplot(
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
	ylab.label = NULL
	);

plot.objects[[length(plot.objects)+1]] <- create.scatterplot(
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
	ylab.label = NULL
	);

create.multipanelplot(
	plot.objects = plot.objects,
	layout.width = 4,
	layout.height = 1,
	plot.objects.widths = c(1,1,1,1),
	left.legend.padding = 0,
	right.legend.padding = 0,
	bottom.legend.padding = 0,
	top.legend.padding = 0,
	height = 5,
	width = 8,
	resolution = 800,
	filename = generate.filename('EVOLVE_ctDNA', 'tumourClearance_stats','png')
	);

### SAVE SESSION INFO ##############################################################################
save.session.profile(generate.filename('SummaryStats','SessionProfile','txt'));
