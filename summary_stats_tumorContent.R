### summary_statistics.R ##########################################################################
# Perform survival statistics on clearance + clinical features

### PREPARE SESSION ################################################################################
# import libraries
library(BoutrosLab.plotting.general);

source('/cluster/home/sprokope/git/analysis/helper_functions/session.functions.R');

working.dir <- '/cluster/projects/ovgroup/projects/OV_Superset/EVOLVE/ctDNA/pipeline_suite/tumour_content';

setwd(working.dir);

### READ DATA ######################################################################################
# find data files
load('../../2022-06-28_EVOLVE_ctDNA_clinicalCovariates.RData');

vaf.data <- read.delim('2022-07-13_EVOLVE_ctDNA_tp53_vafs.tsv');

# get patient order (only if plotting)
patient.order <- as.character(read.delim('patient_order_match_timeline.txt', header = F)$V1);

### FORMAT DATA ####################################################################################
# format data
vaf.data$Patient <- factor(vaf.data$Patient, levels = rev(patient.order));

vaf.data$Group <- 'eot';
vaf.data[which(vaf.data$Timepoint == 2),]$Group <- 'C2';
vaf.data[which(vaf.data$Timepoint == 1),]$Group <- 'baseline';
vaf.data[which(vaf.data$Timepoint == 0),]$Group <- 'baseline';
vaf.data$Group <- factor(vaf.data$Group, levels = c('baseline','C2','eot'), labels = c('baseline','on.trial','EOT'));

vaf.reshaped <- reshape(
	vaf.data[grepl('allUNIQUE', vaf.data$Sample.ct),c('Patient','Group','VAF.ct')],
	direction = 'wide',
	timevar = 'Group',
	idvar = 'Patient'
	);
colnames(vaf.reshaped) <- gsub('VAF\\.ct\\.','',colnames(vaf.reshaped));

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

	for (j in c('COHORT','Best.Response','AGE')) {
	
		ylab <- if (j == 'COHORT') { gsub('\\.',' ', i) } else { NULL }

		xaxis.labels <- if (i == 'EOT') { NA } else { rep('', 4) }

		xlab <- if (i == 'EOT') {
			if (j == 'COHORT') {
				'Cohort' } else if (j == 'AGE') { 
				'Age' } else { 'Best Response' } 
			} else { NULL }

		if (j != 'AGE') {
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
				xaxis.tck = c(0.5,0),
				yaxis.tck = c(0.5,0),
				xlab.cex = 1,
				ylab.cex = 1,
				xaxis.cex = 0.8,
				yaxis.cex = 0.8,
				ylab.label = ylab
				);
			}
		}
	}

create.multipanelplot(
	plot.objects = plot.objects,
	layout.width = 3,
	layout.height = 3,
	plot.objects.widths = c(1,1,1),
	plot.objects.heights = c(1,1,2),
	left.legend.padding = 0,
	right.legend.padding = 0,
	bottom.legend.padding = 0,
	top.legend.padding = 0,
	height = 8,
	width = 8,
	resolution = 800,
	filename = generate.filename('EVOLVE_ctDNA', 'tumourContent_stats','png')
	);

### SAVE SESSION INFO ##############################################################################
save.session.profile(generate.filename('SummaryStats','SessionProfile','txt'));
