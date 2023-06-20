### plot_tumour_content_and_clearance.R ############################################################
# Use allele frequencies from force-called positions

### PREPARE SESSION ################################################################################
# import libraries
library(BoutrosLab.plotting.general);

source('/cluster/home/sprokope/git/analysis/helper_functions/session.functions.R');

# move to working directory
setwd('/cluster/projects/ovgroup/projects/OV_Superset/EVOLVE/ctDNA/pipeline_suite/tumour_content');

load('../2022-06-28_EVOLVE_ctDNA_clinicalCovariates.RData');

### READ DATA ######################################################################################
# get clinical data
clinical <- read.delim('../configs/2022-05-30_sample_info.txt')[,c('Patient','Sample','Tumor.Type')];

clinical$Timepoint <- 0;
clinical[grepl('^C', clinical$Tumor.Type),]$Timepoint <- as.numeric(
	gsub('^C','',gsub('D1','', clinical[grepl('^C', clinical$Tumor.Type),]$Tumor.Type))
	);
clinical[which(clinical$Tumor.Type == 'EOT'),]$Timepoint <- 35;

clinical <- unique(clinical[order(clinical$Patient, clinical$Sample, clinical$Timepoint),]);

# collect list of all samples
all.samples <- as.character(unique(clinical$Sample));

# get bam-readcount output
readcount.files <- list.files(pattern = 'readcounts.tsv', recursive = TRUE);

my.data <- list();
for (i in readcount.files) {
	smp <- sub('_readcounts.tsv','',basename(i));
	my.data[[smp]] <- read.delim(i, header = FALSE)[,c(1:4,6:9)];
	}

# find original positions
known.mutations <- read.delim('EVOLVE_panel_covered_mutations.bed.new', header = FALSE);
colnames(known.mutations) <- c('Chromosome','Start','End','Sample','VAF');

known.mutations$Patient <- substr(known.mutations$Sample, 0, 11);
dup.idx <- which(known.mutations$Patient == 'EVO-009-009' & known.mutations$Start == 7674229);
known.mutations <- known.mutations[-dup.idx,];

### FORMAT DATA ####################################################################################
# format readcount output
combined.data <- do.call(rbind, my.data);
combined.data$Sample <- sapply(rownames(combined.data),function(i) { unlist(strsplit(i,'\\.'))[1] } );
rownames(combined.data) <- 1:nrow(combined.data);
colnames(combined.data) <- c('Chromosome','Position','Ref','Depth','Base.A','Base.C','Base.G','Base.T','Sample');
combined.data$Patient <- substr(combined.data$Sample,0,11);

# information is incomplete for indels, so skip them
combined.data$Depth <- as.numeric(combined.data$Depth);

indel.data <- combined.data[is.na(combined.data$Depth),];
combined.data <- combined.data[!is.na(combined.data$Depth),];

# extract allele counts
combined.data$Count.A <- sapply(
	combined.data$Base.A, function(i) { unlist(strsplit(as.character(i),':'))[2] } );
combined.data$Count.C <- sapply(
	combined.data$Base.C, function(i) { unlist(strsplit(as.character(i),':'))[2] } );
combined.data$Count.G <- sapply(
	combined.data$Base.G, function(i) { unlist(strsplit(as.character(i),':'))[2] } );
combined.data$Count.T <- sapply(
	combined.data$Base.T, function(i) { unlist(strsplit(as.character(i),':'))[2] } );

combined.data$Count.A <- as.numeric(combined.data$Count.A);
combined.data$Count.T <- as.numeric(combined.data$Count.T);
combined.data$Count.C <- as.numeric(combined.data$Count.C);
combined.data$Count.G <- as.numeric(combined.data$Count.G);

# calculate VAF
combined.data$VAF <- NA;

combined.data[which(combined.data$Ref == 'A'),]$VAF <- 1 - (combined.data[which(combined.data$Ref == 'A'),]$Count.A / combined.data[which(combined.data$Ref == 'A'),]$Depth);
combined.data[which(combined.data$Ref == 'T'),]$VAF <- 1 - (combined.data[which(combined.data$Ref == 'T'),]$Count.T / combined.data[which(combined.data$Ref == 'T'),]$Depth);
combined.data[which(combined.data$Ref == 'C'),]$VAF <- 1 - (combined.data[which(combined.data$Ref == 'C'),]$Count.C / combined.data[which(combined.data$Ref == 'C'),]$Depth);
combined.data[which(combined.data$Ref == 'G'),]$VAF <- 1 - (combined.data[which(combined.data$Ref == 'G'),]$Count.G / combined.data[which(combined.data$Ref == 'G'),]$Depth);

# indicate key fields
keep.fields <- c('Patient','Sample','Chromosome','Position','Depth','VAF');

vaf.data <- merge(
	known.mutations,
	combined.data[,keep.fields],
	by.x = c('Patient','Chromosome','Start'),
	by.y = c('Patient','Chromosome','Position'),
	suffixes = c('.wxs','.ct')
	);

vaf.data$Sample <- gsub('_allUNIQUE','', gsub('_DCS', '', gsub('_SSCS','', gsub('_DCS_SSCS','', vaf.data$Sample.ct))));

vaf.data <- merge(
	clinical,
	unique(vaf.data),
	all.x = TRUE
	);

write.table(
	vaf.data,
	file = generate.filename('EVOLVE_ctDNA','tp53_vafs','tsv'),
	row.names = FALSE,
	col.names = TRUE,
	sep = '\t'
	);

# trim to only 1 entry per sample (ie, allUNIQUE)
tmp.data <- unique(vaf.data[grepl('allUNIQUE', vaf.data$Sample.ct),!grepl('wxs',colnames(vaf.data))]);

# calcaluate ctDNA clearance
clearance.data <- data.frame(
	Patient = unique(clinical$Patient),
	T1.ID = NA,
	T2.ID = NA,
	T1 = NA,
	T2 = NA,
	delta = NA
	);

# function to get average percent change
get.delta.stats <- function(t1, t2) {
	if (length(t1) != length(t2)) { next; }
	tmp <- data.frame(T1 = t1, T2 = t2);
	tmp <- tmp[which(tmp$T1 > 0),];
	( 100 * (tmp$T2 - tmp$T1)/tmp$T1 );
	}

for (i in 1:nrow(clearance.data)) {
	smp <- clearance.data[i,]$Patient;
	tmp <- tmp.data[which(tmp.data$Patient == smp),];
	if (nrow(tmp) == 0) { next; }
	clearance.data[i,]$T1 <- min(tmp$Timepoint);
	clearance.data[i,]$T2 <- if (any(tmp$Timepoint == 2)) { 2 } else { max(tmp$Timepoint); }
	clearance.data[i,]$T1.ID <- as.character(tmp[which(tmp$Timepoint == clearance.data[i,]$T1),]$Sample);
	clearance.data[i,]$T2.ID <- as.character(tmp[which(tmp$Timepoint == clearance.data[i,]$T2),]$Sample);
	if (all(is.na(tmp$VAF.ct))) { next; }
	if (all(tmp[tmp$Timepoint %in% clearance.data[i,c('T1')],]$VAF.ct == 0)) { next; }

	# find clearance
	clearance.data[i,]$delta <- get.delta.stats(
		t1 = tmp.data[which(tmp.data$Sample == clearance.data[i,]$T1.ID),]$VAF.ct,
		t2 = tmp.data[which(tmp.data$Sample == clearance.data[i,]$T2.ID),]$VAF.ct
		);
	}

clearance.data[which(clearance.data$T1 == clearance.data$T2),]$delta <- NA;

write.table(
	clearance.data,
	file = generate.filename('EVOLVE_ctDNA','tumour_clearance','tsv'),
	row.names = FALSE,
	col.names = TRUE,
	sep = '\t'
	);

# for plotting
plot.data <- vaf.data;

clinical <- clinical[order(clinical$Patient, clinical$Timepoint),];
plot.data$Sample <- factor(plot.data$Sample, levels = clinical$Sample);

plot.data <- unique(plot.data[,!grepl('wxs', colnames(plot.data))]);
plot.data <- plot.data[grepl('allUNIQUE',plot.data$Sample.ct) | is.na(plot.data$Sample.ct),];
plot.data <- plot.data[order(plot.data$Sample),];

# fix the timing names
plot.data$Type <- 'eot';
plot.data[which(plot.data$Tumor.Type %in% c('Screening','C1D1')),]$Type <- 'baseline';
plot.data[which(plot.data$Tumor.Type == 'C2D1'),]$Type <- 'cycle2';

all.samples <- as.character(plot.data$Sample);

### SET-UP COVARIATES ##############################################################################
# set up covariates
type.colours <- rep('seagreen', length(all.samples));
type.colours[which(plot.data$Type == 'baseline')] <- 'skyblue';
type.colours[which(plot.data$Type == 'eot')] <- 'pink';

covariates.obj <- covariates.grob(
	covariates = list(
		rect = list(col = 'transparent', fill = type.colours, lwd = 1)
		),
	ord = 1:nrow(plot.data),
	side = 'top',
	grid.col = list(col = 'black', lwd = 1),
	col.lines = get.line.breaks(plot.data$Patient)-0.5,
	grid.border = list(col = 'black', lwd = 1)
	);

# make the plot legend (mutation type/consequence)
functional.legend <- list(
	legend = list(
		colours = c('skyblue','seagreen','pink'),
		labels = c('baseline','on-trial','end-of-trial')
		)
	);

clinical.legend <- legend.grob(
	legends = functional.legend,
	title.just = 'left',
	label.cex = 0.7,
	size = 1.2
	);

### PLOTTING #######################################################################################
# plot VAFs
create.scatterplot(
	VAF.ct*100 ~ Sample,
	plot.data,
	type = c('p','h'),
	xlab.label = NULL,
	ylab.label = 'VAF (%)',
	ylab.axis.padding = 2,
	xaxis.rot = 90,
	xaxis.cex = 0.8,
	xaxis.lab = gsub('_ctDNA','',plot.data$Sample),
	ylimits = c(-5,80),
	yat = seq(0,80,20),
	legend = list(
		inside = list(fun = covariates.obj, x = 0.5, y = 0.03),
		inside = list(fun = clinical.legend, x = 0.02, y = 0.98)
		),
	height = 5,
	width = 12,
	style = 'Nature',
	filename = generate.filename('EVOLVE_ctDNA','forced_allele_vafs','png')
	);

# plot clearance
timecourse.data <- merge(
	plot.data[,c('Patient','Sample','Timepoint')],
	phenodata,
	by.x = 'Patient',
	by.y = 'Patient.ID'
	);

timecourse.data$ctDNA <- NA;
for (smp in unique(timecourse.data$Patient)) {
	tmp <- plot.data[which(plot.data$Patient == smp),];
	if (all(is.na(tmp$VAF.ct))) { next; }
	deltas <- diff(tmp$VAF.ct)/tmp$VAF.ct[1:(nrow(tmp)-1)];
	timecourse.data[which(timecourse.data$Patient == smp),]$ctDNA <- c(0,deltas)*100;
	}

group.key1 <- list(
	points = list(pch = c(19,15,17,18), cex = 1),
	text = list(lab = covariate.colours$Response$labels, cex = 1),
	title = 'Best Response             ',
	cex.title = 1
	);

group.key2 <- list(
	lines = list(lwd = 1, size = 2, type = 'b', col = covariate.colours$Cohort$colours, pch = 19, cex = 0.5),
	text = list(lab = covariate.colours$Cohort$labels, cex = 1),
	title = 'Cohort                         ',
	cex.title = 1,
	divide = 1
	);

timecourse.data$Time <- 3;
timecourse.data[which(timecourse.data$Timepoint %in% c(0,1)),]$Time <- 1;
timecourse.data[which(timecourse.data$Timepoint == 2),]$Time <- 2;
timecourse.data$Time <- factor(timecourse.data$Time, levels = c(1,2,3), labels = c('baseline','on-trial','EOT'));

create.scatterplot(
	ctDNA ~ log10(timecourse.data$Timepoint + 0.1),
	timecourse.data,
	groups = timecourse.data$Patient,
	type = 'b',
	col = covariate.colours$Cohort$colours[match(timecourse.data[!duplicated(timecourse.data$Patient),]$Cohort.cat, covariate.colours$Cohort$levels)],
	cex = 1,
	alpha = 0.9,
	pch = c(19,15,17,18)[match(timecourse.data[!duplicated(timecourse.data$Patient),]$Response.cat, covariate.colours$Response$levels)],
	xlab.label = 'Timepoint',
	ylab.label = expression(Delta * 'ctDNA'),
	ylab.axis.padding = 2,
	ylimits = c(-100,500),
	yat = seq(-100,500,100),
	xlimits = log10(c(0.1,2.1,5.1,10.1,37.1)), #c(0,36),
	xat = log10(c(0.1,2.1,5.1,10.1,35.1)), #seq(0,35,5),
	xaxis.lab = c(0,2,5,10,'EOT'), #c(seq(0,30,5),'EOT'),
	abline.h = 0,
	abline.lty = 2,
	abline.col = 'grey70',
	legend = list(
		inside = list(fun = draw.key, args = list(key = group.key2), x = 0.01, y = 1.05),
		inside = list(fun = draw.key, args = list(key = group.key1), x = 0.01, y = 0.81)
		),
	height = 5,
	width = 6,
	style = 'Nature',
	filename = generate.filename('EVOLVE_ctDNA','forced_allele_percent_change_timecourse','png')
	);

### NEW 20220706 ###
patient.order <- as.character(read.delim('patient_order_match_timeline.txt', header = F)$V1);
plot.data$Patient <- factor(plot.data$Patient, levels = rev(patient.order));

plot.data$Group <- 'eot';
plot.data[which(plot.data$Timepoint == 2),]$Group <- 'C2';
plot.data[which(plot.data$Timepoint == 1),]$Group <- 'baseline';
plot.data[which(plot.data$Timepoint == 0),]$Group <- 'baseline';
plot.data$Group <- factor(plot.data$Group, levels = c('baseline','C2','eot'), labels = c('baseline','on-trial','EOT'));

# create barplot
x <- reshape(plot.data[,c('Patient','Group','VAF.ct')], direction = 'wide', timevar = 'Group', idvar = 'Patient');
colnames(x) <- gsub('VAF\\.ct\\.','',colnames(x));
colnames(x)[3] <- 'C2';

plot.objects <- list();
plot.objects[[1]] <- create.barplot(
	Patient ~ baseline,
	x,
	plot.horizontal = TRUE,
	col = 'grey90',
	yaxis.cex = 1,
	xaxis.cex = 1,
	xaxis.tck = c(0.5,0),
	yaxis.tck = c(0.5,0),
	xlab.label = NULL,
	ylab.label = NULL,
	style = 'Nature'
	);
plot.objects[[2]] <- create.barplot(
	Patient ~ C2,
	x,
	plot.horizontal = TRUE,
	col = 'grey90',
	yaxis.cex = 1,
	xaxis.cex = 1,
	xaxis.tck = c(0.5,0),
	yaxis.tck = c(0,0),
	xlab.label = NULL,
	ylab.label = NULL,
	style = 'Nature'
	);
plot.objects[[3]] <- create.barplot(
	Patient ~ EOT,
	x,
	plot.horizontal = TRUE,
	col = 'grey90',
	yaxis.cex = 1,
	xaxis.cex = 1,
	xaxis.tck = c(0.5,0),
	yaxis.tck = c(0,0),
	xlab.label = NULL,
	ylab.label = NULL,
	style = 'Nature'
	);

create.multiplot(
	plot.objects = plot.objects,
	panel.widths = c(1,1,1),
	xaxis.fontface = 'plain',
	yaxis.fontface = 'plain',
	xaxis.tck = c(0.5,0),
	yaxis.tck = c(0,0),
	xaxis.cex = 1,
	yaxis.cex = 0.8,
	axes.lwd = 1,
	plot.layout = c(3,1),
	remove.all.border.lines = FALSE,
	xlab.label = 'Estimated Tumour Fraction',
	xlab.cex = 1.5,
	print.new.legend = TRUE,
	legend = list(
		inside = list(fun = draw.key, args = list(key = list(text = list(c('baseline','on-trial','EOT'), cex = 1), columns = 3, between = 4)), x = 0.01, y = 1.04)
		),
	style = 'Nature',
	width = 5.5,
	filename = generate.filename('EVOLVE_ctDNA', 'tumourContent_plot','png')
	);

####################

### HEATMAP ###
# create heatmap showing exome and ctDNA TP53 VAFs
known.mutations$Group <- substr(known.mutations$Sample, 13, 14);

tmp <- reshape(
	known.mutations[,-4],
	direction = 'wide',
	timevar = 'Group',
	idvar = c('Chromosome','Start','End','Patient')
	);
rownames(tmp) <- 1:nrow(tmp);
known <- tmp;

ctdna <- vaf.data;
ctdna <- ctdna[order(ctdna$Timepoint),];
ctdna$Tumor.Type <- factor(ctdna$Tumor.Type, levels = unique(ctdna$Tumor.Type));

tmp <- reshape(
	unique(ctdna[grepl('allUNIQUE', ctdna$Sample.ct),c('Patient','Chromosome','Start','End','Tumor.Type','VAF.ct')]),
	direction = 'wide',
	idvar = c('Patient','Chromosome','Start','End'),
	timevar = 'Tumor.Type'
	);

# join for plotting
plot.data <- merge(known, tmp, by.x = colnames(known)[1:4], by.y = c('Chromosome','Start','End','Patient'));
rownames(plot.data) <- plot.data$Patient;

# calculate clearance
plot.data$delta <- apply(plot.data[,grepl('VAF.ct', colnames(plot.data))],1,function(i) { i <- na.omit(i); if (length(i) < 2) { return(0) } else { return( 100* ( (i[2]-i[1]) / i[1]) ) } } );

# sort data
plot.data <- plot.data[order(plot.data$delta, decreasing = FALSE),];
plot.data$Patient <- factor(plot.data$Patient, levels = as.character(plot.data$Patient));
plot.data$Order <- 1:nrow(plot.data);

tissue <- plot.data[,c('VAF.Ar','VAF.Bx','VAF.Pr')];
blood <- plot.data[,grepl('VAF.ct', colnames(plot.data))];

point.colours <- rep('cornflowerblue', nrow(plot.data));
point.colours[is.na(plot.data$delta)] <- NA;
point.colours[which(plot.data$delta > 0)] <- 'lightcoral';

# make plots
tissue.plot <- create.heatmap(
	t(plot.data[,c('VAF.Ar','VAF.Bx','VAF.Pr'),]),
	cluster.dimensions = 'none',
	yaxis.lab = NA,
	xaxis.lab.top = c('diagnostic','baseline','progression'),
	xaxis.top.cex = 1,
	xaxis.cex = 1,
	yaxis.cex = 0.8,
	xaxis.tck = 0,
	yaxis.tck = 0.2,
	xaxis.rot.top = 60,
	xaxis.fontface = 'plain',
	yaxis.fontface = 'plain',
	x.alternating = 2,
	axes.lwd = 1,
	grid.row = TRUE,
	force.grid.row = TRUE,
	row.colour = 'grey80',
	col.colour = 'grey80',
	row.lwd = 1, #if (length(all.samples) < 30) { 3 } else { 1 },
	col.lwd = 1, #if (length(all.samples) < 30) { 3 } else { 1 },
	grid.col = TRUE, #(length(all.samples) > 1),
	force.grid.col = TRUE, #(length(all.samples) > 1),
	print.colour.key = TRUE,
	colourkey.cex = 0.8,
	colourkey.labels.at = seq(0,1,0.5),
	colourkey.labels = c('0','0.5','1.0'),
	fill.colour = 'grey90',
	at = seq(0,1,0.01),
	colour.scheme = c('white','darkred'),
	cell.text = sprintf("%.2f", round(tissue[which(tissue > 0, arr.ind = T)],3)),
	text.cex = 0.6,
	text.col = sapply(
		tissue[which(tissue > 0, arr.ind = T)], function(i) { if (i > 0.5) { 'white' } else { 'black' } }),
	row.pos = which(tissue > 0, arr.ind = T)[,1],
	col.pos = which(tissue > 0, arr.ind = T)[,2]
	);

tissue.plot$xlab$label <- 'VAF';
tissue.plot$xlab$cex <- 1;
tissue.plot$xlab$fontface <- 'plain';

blood.plot <- create.heatmap(
	t(plot.data[,grepl('VAF.ct', colnames(plot.data))]),
	cluster.dimensions = 'none',
	yaxis.lab = rep('', nrow(plot.data)), 
	xaxis.lab.top = sub('Screening','baseline',gsub('VAF.ct.','',colnames(plot.data)[grepl('VAF.ct', colnames(plot.data))])),
	xaxis.top.cex = 1,
	xaxis.cex = 1,
	xaxis.tck = 0,
	yaxis.tck = 0,
	xaxis.rot.top = 60,
	xaxis.fontface = 'plain',
	x.alternating = 2,
	axes.lwd = 1,
	grid.row = TRUE,
	force.grid.row = TRUE,
	row.colour = 'grey80',
	col.colour = 'grey80',
	row.lwd = 1, #if (length(all.samples) < 30) { 3 } else { 1 },
	col.lwd = 1, #if (length(all.samples) < 30) { 3 } else { 1 },
	grid.col = TRUE, #(length(all.samples) > 1),
	force.grid.col = TRUE, #(length(all.samples) > 1),
	print.colour.key = TRUE,
	colourkey.cex = 0.8,
	colourkey.labels.at = seq(0,0.102,0.01),
	colourkey.labels = c(seq(0,0.09,0.01),'0.1+'),
	fill.colour = 'grey90',
	at = seq(0,0.101,0.001),
	colour.scheme = c('white','red','darkred'),
	colour.centering.value = 0.1,
	cell.text = sprintf("%.3f", round(blood[which(blood > 0, arr.ind = T)],3)),
	text.cex = 0.6,
	text.col = sapply(
		blood[which(blood > 0, arr.ind = T)], function(i) { if (i > 0.1) { 'white' } else { 'black' } }),
	row.pos = which(blood > 0, arr.ind = T)[,1],
	col.pos = which(blood > 0, arr.ind = T)[,2]
	);

blood.plot$xlab$label <- 'Estimated Tumor Fraction';
blood.plot$xlab$cex <- 1;
blood.plot$xlab$fontface <- 'plain';

# plot clearance
clearance.plot <- create.barplot(
	Order ~ delta,
	plot.data,
	col = point.colours[!is.na(point.colours)],
	ylab.label = NULL,
	xlab.label = expression(Delta * 'ctDNA'),
	xlab.cex = 1,
	xaxis.cex = 0.8,
	xlimits = c(-100,300),
	xat = c(-100,0,100,300),
	xaxis.fontface = 'plain',
	yaxis.lab = rep('',nrow(plot.data)),
	xaxis.tck = c(0.2,0),
	yaxis.tck = 0,
	plot.horizontal = TRUE,
	abline.v = 0
	);

clearance.plot$par.settings$axis.line$col <- 'transparent';

# Re-add bottom and left axes
clearance.plot$axis <- function(side, line.col = 'black', ...) {
	# Only draw axes on the bottom
	if (side %in% c('bottom')) {
		axis.default(side = side, line.col = 'black', ...);
		lims <- current.panel.limits();
		panel.abline(h = lims$ylim[1], v = lims$ylim[1]);
		}
	}

create.multipanelplot(
	plot.objects = list(tissue.plot, blood.plot, clearance.plot),
	layout.width = 3,
	layout.height = 1,
	plot.objects.widths = c(1.8,4,1.2),
	x.spacing = c(0.8,0),
	left.legend.padding = 0,
	right.legend.padding = 0,
	bottom.legend.padding = 1,
	top.legend.padding = 0,
	top.padding = 4,
	height = 6,
	width = 8,
	resolution = 800,
	filename = generate.filename('EVOLVE_ctDNA', 'tp53_VAFs','png')
	);

write.table(
	plot.data,
	file = generate.filename('EVOLVE_ctDNA','tumour_clearance_plotData','tsv'),
	row.names = FALSE,
	col.names = TRUE,
	sep = '\t'
	);

### SAVE SESSION INFO ##############################################################################
save.session.profile(generate.filename('tumourClearance','SessionProfile','txt'));
