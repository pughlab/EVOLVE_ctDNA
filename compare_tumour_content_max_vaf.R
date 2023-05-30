### compare_tumour_content_max_vaf.R ###############################################################
# Compare tumour content estimation (by TP53 mutation with bam-readcount) and by maximum VAF

### PREPARE SESSION ################################################################################
# import libraries
library(BoutrosLab.plotting.general);

source('/cluster/home/sprokope/git/analysis/helper_functions/session.functions.R');

input.dir <- '/cluster/projects/ovgroup/projects/OV_Superset/EVOLVE/ctDNA/pipeline_suite/Ensemble_calls';
output.dir <- 'tumour_content';

setwd(input.dir);

### READ DATA ######################################################################################
# get clinical data (for ctDNA)
timeline <- read.delim('../configs/2022-09-06_sample_info_with_batch.txt', stringsAsFactors = FALSE);

all.patients <- unique(timeline$Patient);

# get clinical data (for WXS)
phenodata <- read.delim('../configs/EVOLVE_clinical_data_updated2022.tsv', stringsAsFactors = FALSE);
phenodata <- phenodata[which(phenodata$Patient.ID %in% all.patients),c('Patient.ID','Sample.ID')];

phenodata$Sample.ID <- gsub('_', '-', phenodata$Sample.ID);
phenodata$Group <- substr(phenodata$Sample.ID, 13, 14);

# find data files
maf <- read.delim('2022-07-26_EVOLVE_ctdna_plus_exome__mutation_data_cbioportal.txt');

# read in estimated tumour content
purity <- read.delim('../tumour_content/2022-07-13_EVOLVE_ctDNA_tp53_vafs.tsv');

setwd(output.dir);

# clean up clinical data
clin.data <- merge(
	phenodata,
	timeline[,c('Patient','Sample','Timepoint','Group')],
	by.x = 'Patient.ID',
	by.y = 'Patient',
	suffixes = c('.wxs','.ctdna'),
	all = TRUE
	);

clin.data <- clin.data[-which(clin.data$Sample.ID == 'EVO-009-018-Bx'),];

### FORMAT DATA ####################################################################################
# format purity data
purity <- unique(purity[grepl('allUNIQUE', purity$Sample.ct),]);
colnames(purity)[12] <- 'Purity';

purity.formatted <- merge(
	clin.data,
	purity[,-4],
	by.x = c('Patient.ID','Sample.ID','Sample'),
	by.y = c('Patient','Sample.wxs','Sample'),
	all.x = TRUE
	);

purity.formatted <- purity.formatted[,c(colnames(clin.data),colnames(purity)[c(5:7,9,12)])];

# coverage in exome was too low to call confidently, but I suspect this is germline
maf[which(grepl('EVO-400-004', maf$Tumor_Sample_Barcode) & maf$Hugo_Symbol == 'BRCA2'),]$Mutation_Status <- 'germline';

# remove germline events
maf <- maf[which(maf$Mutation_Status != 'germline'),];

# format mutation data
maf$Patient <- substr(maf$Tumor_Sample_Barcode, 0, 11);
maf$t_vaf <- maf$t_alt_count / maf$t_depth;

keep.fields <- c('Patient','Hugo_Symbol','Chromosome','Start_Position','End_Position','Variant_Classification','Variant_Type','Tumor_Sample_Barcode','HGVSc','HGVSp_Short','t_vaf');

exome.idx <- !grepl('ctDNA', maf$Tumor_Sample_Barcode);
ctdna.idx <- grepl('ctDNA', maf$Tumor_Sample_Barcode);

merged <- merge(
	maf[exome.idx,keep.fields],
	maf[ctdna.idx,keep.fields],
	by = c('Patient','Hugo_Symbol','Chromosome','Start_Position','End_Position','Variant_Classification','Variant_Type','HGVSc','HGVSp_Short'),
	suffixes = c('.wxs','.ct'),
	all = TRUE
	);

# find maximum VAF for detected variants
purity.formatted$Max.VAF <- NA;
for (patient in unique(purity.formatted$Patient.ID)) {

	tmp.data <- maf[which(maf$Patient == patient & grepl('ctDNA', maf$Tumor_Sample_Barcode)),];
	if (nrow(tmp.data) == 0) { next; }

	tmp.data <- tmp.data[,keep.fields];

	# try to filter to variants present in all samples
	n.smps <- length(unique(tmp.data$Tumor_Sample_Barcode));
	tmp <- aggregate(
		Tumor_Sample_Barcode ~ Chromosome + Start_Position + End_Position,
		tmp.data,
		length
		);

	if (any(tmp$Tumor_Sample_Barcode == n.smps)) {
		tmp <- tmp[which(tmp$Tumor_Sample_Barcode == n.smps),];
		tmp.data <- merge(tmp.data, tmp[,c('Chromosome','Start_Position','End_Position')]);
		}

	# if there is a baseline sample, get the coordinates for the max VAF position
	if (any(grepl('Screening|C1D1', tmp.data$Tumor_Sample_Barcode))) {
		smp.data <- tmp.data[grepl('Screening|C1D1', tmp.data$Tumor_Sample_Barcode),];
		} else {
		# otherwise check all samples
		smp.data <- tmp.data;
		}

	smp.data <- smp.data[order(smp.data$t_vaf, decreasing = TRUE),];
	coordinates <- smp.data[1,c('Chromosome','Start_Position','End_Position')];

	tmp.data <- merge(tmp.data, coordinates);

	for (smp in tmp.data$Tumor_Sample_Barcode) {
		max.vaf <- tmp.data[which(tmp.data$Tumor_Sample_Barcode == smp),]$t_vaf[1];
		purity.formatted[which(purity.formatted$Sample == smp),]$Max.VAF <- max.vaf;
		}

	rm(tmp.data,smp.data,coordinates,max.vaf,smp,tmp,n.smps);
	}

# combine results
results <- unique(purity.formatted[,c('Patient.ID','Sample','Timepoint','Group.ctdna','Purity','Max.VAF')]);
results$Final <- results$Purity;
results[is.na(results$Purity),]$Final <- results[is.na(results$Purity),]$Max.VAF;

# save results
write.table(
	results,
	file = generate.filename('EVOLVE_ctDNA', '_estimated_tumour_content','tsv'),
	row.names = FALSE,
	col.names = TRUE,
	sep = '\t'
	);

### PLOT DATA ######################################################################################
# first, let's compare per-position VAFs between WXS and ctDNA
plot.data <- merge(
	merged[!is.na(merged$t_vaf.wxs) & !is.na(merged$t_vaf.ct),c(1,10:13)],
	clin.data,
	by.x = c('Patient','Tumor_Sample_Barcode.wxs','Tumor_Sample_Barcode.ct'),
	by.y = c('Patient.ID','Sample.ID','Sample'),
	all.x = TRUE
	);

# create plot
create.scatterplot(
	t_vaf.ct ~ t_vaf.wxs,
	plot.data,
	cex = c(0.5,0.5,1.5)[match(plot.data$Group.ctdna, c('EOT','on.trial','baseline'))],
	col = c('grey70','grey70','black')[match(plot.data$Group.ctdna, c('EOT','on.trial','baseline'))],
	alpha = 0.8,
	xlimits = c(0,1),
	xat = seq(0,1,0.2),
	ylimits = c(0,0.8),
	yat = seq(0,1,0.1),
	yaxis.tck = c(0.75,0),
	xaxis.tck = c(0.75,0),
	ylab.label = 'VAF (cfDNA)',
	xlab.label = 'VAF (exome)',
	ylab.cex = 1.5,
	xlab.cex = 1.5,
	yaxis.cex = 1,
	xaxis.cex = 1,
	yaxis.fontface = 'plain',
	xaxis.fontface = 'plain',
	key = get.corr.key(
		y = plot.data[which(plot.data$Group.ctdna == 'baseline'),]$t_vaf.ct,
		x = plot.data[which(plot.data$Group.ctdna == 'baseline'),]$t_vaf.wxs,
		y.pos = 0.95,
		x.pos = 0.05,
		key.cex = 1.5
		),
	filename = generate.filename('EVOLVE_ctDNA','_exome_vaf__vs__ctdna_vaf','png'),
	style = 'Nature'
	);

# next, let's compare purity using TP53 VAF and max VAF
results$Purity.mod <- results$Purity;
results[is.na(results$Purity),]$Purity.mod <- -0.07;
results$Max.VAF.mod <- results$Max.VAF;
results[is.na(results$Max.VAF),]$Max.VAF.mod <- -0.07;

# create plot
create.scatterplot(
	Purity.mod ~ Max.VAF.mod,
	results,
	alpha = 0.8,
	ylimits = c(-0.1,1),
	yat = c(-0.07, seq(0,1,0.2)),
	yaxis.lab = c('NA', '0','0.2','0.4','0.6','0.8','1.0'),
	xlimits = c(-0.1,1),
	xat = c(-0.07, seq(0,1,0.2)),
	xaxis.lab = c('NA', '0','0.2','0.4','0.6','0.8','1.0'),
	yaxis.tck = c(0.75,0),
	xaxis.tck = c(0.75,0),
	ylab.label = 'Estimated Tumour Fraction',
	xlab.label = 'Maximum VAF',
	ylab.cex = 1.5,
	xlab.cex = 1.5,
	yaxis.cex = 1,
	xaxis.cex = 1,
	abline.h = 0,
	abline.v = 0,
	abline.lty = 2,
	abline.lwd = 1,
	add.rectangle = TRUE,
	xleft.rectangle = c(-0.1,0),
	xright.rectangle = c(0,1),
	ytop.rectangle = c(1,0),
	ybottom.rectangle = c(0,-0.1),
	col.rectangle = c('grey50','grey50'),
	alpha.rectangle = 0.5,
	key = get.corr.key(
		x = results$Max.VAF,
		y = results$Purity,
		label.items = c('pearson','pearson.p'),
		x.pos = 0.16,
		y.pos = 0.95,
		key.cex = 1.5
		),
	file = generate.filename('EVOLVE_ctDNA','_tp53_vaf__vs__max_vaf','png'),
	style = 'Nature'
	);

# now, re-create heatmap showing exome and ctDNA VAFs
# find original positions
known.mutations <- read.delim('../../tumour_content/EVOLVE_panel_covered_mutations.bed.new', header = FALSE);
colnames(known.mutations) <- c('Chromosome','Start','End','Sample','VAF');

known.mutations <- merge(clin.data, known.mutations, by.x = 'Sample.ID', by.y = 'Sample', all.x = TRUE);

dup.idx <- which(known.mutations$Patient == 'EVO-009-009' & known.mutations$Start == 7674229);
known.mutations <- known.mutations[-dup.idx,];

# figure out which positions to plot
tmp <- merge(
	known.mutations[,c('Patient.ID','Sample.ID','Group.wxs','VAF','Sample','Timepoint','Group.ctdna')],
	results[,c('Patient.ID','Sample','Purity','Final')],
	by.x = c('Patient.ID','Sample'),
	by.y = c('Patient.ID','Sample'),
	all = TRUE
	);

# format for plotting
plot.data <- merge(
	reshape(unique(tmp[,c(1,4,5)]), direction = 'wide', idvar = 'Patient.ID', timevar = 'Group.wxs'),
	reshape(unique(tmp[,c(1,7,9)]), direction = 'wide', idvar = 'Patient.ID', timevar = 'Group.ctdna')
	);

colnames(plot.data) <- gsub('Final', 'VAF.ct', colnames(plot.data));
rownames(plot.data) <- plot.data$Patient.ID;
plot.data <- plot.data[,c('VAF.Ar','VAF.Bx','VAF.Pr','VAF.ct.baseline','VAF.ct.on.trial','VAF.ct.EOT')];

# calculate clearance
plot.data$delta <- apply(
	plot.data[,c('VAF.ct.baseline','VAF.ct.on.trial')],
	1,
	function(i) {
		i <- na.omit(i);
		if (length(i) < 2) { return(NA) } else { return( 100* ( (i[2]-i[1]) / i[1]) ) }
		}
	);

plot.data$max <- apply(
	plot.data[,grepl('VAF.ct', colnames(plot.data))],
	1,
	function(i) {
		if (all(is.na(i))) { return(NA) } else { max(i, na.rm = TRUE) }
		}	
	);

# sort data
plot.data <- plot.data[order(plot.data$delta, plot.data$max, decreasing = FALSE, na.last = FALSE),];
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

tissue.plot$xlab$label <- 'VAF (WES)';
tissue.plot$xlab$cex <- 1;
tissue.plot$xlab$fontface <- 'plain';

blood.plot <- create.heatmap(
	t(plot.data[,grepl('VAF.ct', colnames(plot.data))]),
	cluster.dimensions = 'none',
	yaxis.lab = rep('', nrow(plot.data)), 
	xaxis.lab.top = c('baseline','on-trial','end-of-treatment'),
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
	colourkey.labels.at = c(0, 0.02, 0.05, 0.1), #seq(0,0.102,0.02),
	colourkey.labels = c(0,0.02,0.05,'0.1+'),
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

blood.plot$xlab$label <- 'Estimated\nTumour Fraction';
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

# combine them
create.multipanelplot(
	plot.objects = list(tissue.plot, blood.plot, clearance.plot),
	layout.width = 3,
	layout.height = 1,
	plot.objects.widths = c(2,1.3,1.2),
	x.spacing = c(0.3,0),
	left.legend.padding = 0,
	right.legend.padding = 0,
	bottom.legend.padding = 0,
	top.legend.padding = 0,
	top.padding = 5,
	height = 6,
	width = 5,
	resolution = 800,
	filename = generate.filename('EVOLVE_ctDNA', '_estimated_tumour_contents','png')
	);

# create covariates
load('../../2022-09-06_EVOLVE_ctDNA_clinicalCovariates.RData');

plot.data.anno <- merge(clinical[,c(1,9:11,14)], plot.data, by.x = 'Patient.ID', by.y = 'row.names');
plot.data.anno <- plot.data.anno[order(plot.data.anno$Order),];

smp.covariates <- list(
	rect = list(
		col = 'transparent',
		fill = covariate.colours$Cohort$colours[match(plot.data.anno$Cohort.cat, covariate.colours$Cohort$levels)],
		lwd = 1
		),
	rect = list(
		col = 'transparent',
		fill = covariate.colours$BRCA$colours[match(plot.data.anno$BRCA.cat, covariate.colours$BRCA$levels)],
		lwd = 1
		),
	rect = list(
		col = 'transparent',
		fill = covariate.colours$Age$colours[match(plot.data.anno$Age.cat, covariate.colours$Age$levels)],
		lwd = 1
		),
	rect = list(
		col = 'transparent',
		fill = covariate.colours$Response$colours[match(plot.data.anno$Response.cat, covariate.colours$Response$levels)],
		lwd = 1
		)
	);

smp.covariate.grob <- covariates.grob(
	covariates = smp.covariates,
	ord = 1:nrow(plot.data.anno),
	side = 'right',
	size = 1,
	grid.col = list(col = 'white', lwd = 2),
	col.lines = 1:length(smp.covariates)
	);

# plot clearance separately
create.barplot(
	Order ~ delta,
	plot.data.anno,
	col = point.colours[!is.na(point.colours)],
	ylab.label = NULL,
	xlab.label = expression(Delta * 'ctDNA'),
	xlab.cex = 1.5,
	xaxis.cex = 1.2,
	yaxis.cex = 1.2,
	xlimits = c(-120,300),
	xat = c(-100,0,100,300),
	yaxis.lab = rep('', nrow(plot.data.anno)),
	xaxis.fontface = 'plain',
	yaxis.fontface = 'plain',
	xaxis.tck = c(1,0),
	yaxis.tck = 0,
	plot.horizontal = TRUE,
	legend = list(
		inside = list(fun = smp.covariate.grob, x = -0.3, y = 0.5)
		),
	left.padding = 8,
	height = 9,
	width = 4,
	style = 'Nature',
	filename = generate.filename('EVOLVE_ctDNA', '_estimated_tumour_clearance','png')
	);

# try a scatterplot per patient?
plot.data <- results[order(results$Patient.ID, results$Timepoint),];
plot.data$Group.ctdna <- factor(
	plot.data$Group.ctdna,
	levels = c('baseline','on.trial','EOT'),
	labels = c('baseline','on-trial','end-of-treatment')
	);

nv.key <- list(
	text = list(lab = 'no data', cex = 1)
	);

create.scatterplot(
	log10(Final) ~ Group.ctdna | Patient.ID,
	plot.data,
	type = 'b',
	xaxis.rot = 90,
	x.spacing = 1,
	y.spacing = 1,
	as.table = TRUE,
	xaxis.tck = c(1,0),
	yaxis.tck = c(1,0),
	xlab.label = NULL,
	ylab.label = expression('Estimated Tumour Fraction'),
	ylab.axis.padding = 2,
	yat = log10(c(0.001,0.01,0.1,1)),
	yaxis.lab = c(0.001,0.01,0.1,1),
	xaxis.cex = 1.2,
	yaxis.cex = 1.2,
	ylab.cex = 1.5,
	xaxis.fontface = 'plain',
	yaxis.fontface = 'plain',
	strip.fontface = 'plain',
	strip.cex = 0.9,
	legend = list(
		inside = list(fun = draw.key, args = list(key = nv.key), x = 0.187, y = 0.5),
		inside = list(fun = draw.key, args = list(key = nv.key), x = 0.187, y = 0.09),
		inside = list(fun = draw.key, args = list(key = nv.key), x = 0.53, y = 0.09)
		),
	height = 9,
	width = 9,
	filename = generate.filename('EVOLVE_ctDNA', '_estimated_tumour_fraction','png')
	);









### OTHER ####
exome.purity <- read.delim('/cluster/projects/ovgroup/projects/OV_Superset/Combined/EVOLVE/cBioportal/upload_20211130/sample_clinical_data.txt', skip = 4);
exome.purity <- exome.purity[,c('PATIENT_ID','SAMPLE_ID','CELLULARITY')];
colnames(exome.purity) <- c('Patient','Tumor_Sample_Barcode','Cellularity');

tmp1 <- merge(
	maf[exome.idx,keep.fields],
	exome.purity
	);

colnames(purity) <- c('Patient','Tumor_Sample_Barcode','Cellularity');
tmp2 <- merge(
	maf[ctdna.idx,keep.fields],
	purity,
	all.x = TRUE
	);

tmp1$ccf <- tmp1$t_vaf / tmp1$Cellularity;
tmp2$ccf <- tmp2$t_vaf / tmp2$Cellularity;

merged <- merge(
	tmp1[,c(keep.fields,'ccf')],
	tmp2[,c(keep.fields,'ccf')],
	by = c('Patient','Hugo_Symbol','Chromosome','Start_Position','End_Position','Variant_Classification','Variant_Type','HGVSc','HGVSp_Short'),
	suffixes = c('.wxs','.ct'),
	all = TRUE
	);

results <- merged[!is.na(merged$Tumor_Sample_Barcode.ct),];
cor(results[,c('ccf.wxs','ccf.ct')], use = 'pairwise');


tmp <- results[,c('Patient','Hugo_Symbol','Tumor_Sample_Barcode.wxs','ccf.wxs','Tumor_Sample_Barcode.ct','t_vaf.ct','ccf.ct')];

tmp$ccf.wxs <- round(tmp$ccf.wxs, 3);
tmp$ccf.ct <- round(tmp$ccf.ct, 3);
tmp$t_vaf.ct <- round(tmp$t_vaf.ct, 3);

tmp$Tumor_Sample_Barcode.ct <- gsub('_ctDNA','', tmp$Tumor_Sample_Barcode.ct)
colnames(tmp)[3] <- 'WXS';
colnames(tmp)[5] <- 'ctDNA';

version1 <- tmp;

### ALTERNATE ###
colnames(max.vaf) <- c('Patient','Tumor_Sample_Barcode','Cellularity');
tmp2 <- merge(
	maf[ctdna.idx,keep.fields],
	max.vaf,
	all.x = TRUE
	);

tmp2$ccf <- tmp2$t_vaf / tmp2$Cellularity;

merged <- merge(
	tmp1[,c(keep.fields,'ccf')],
	tmp2[,c(keep.fields,'ccf')],
	by = c('Patient','Hugo_Symbol','Chromosome','Start_Position','End_Position','Variant_Classification','Variant_Type','HGVSc','HGVSp_Short'),
	suffixes = c('.wxs','.ct'),
	all = TRUE
	);

results <- merged[!is.na(merged$Tumor_Sample_Barcode.ct),];
cor(results[,c('ccf.wxs','ccf.ct')], use = 'pairwise');

tmp <- results[,c('Patient','Hugo_Symbol','Tumor_Sample_Barcode.wxs','ccf.wxs','Tumor_Sample_Barcode.ct','t_vaf.ct','ccf.ct')];

tmp$ccf.wxs <- round(tmp$ccf.wxs, 3);
tmp$ccf.ct <- round(tmp$ccf.ct, 3);
tmp$t_vaf.ct <- round(tmp$t_vaf.ct, 3);

tmp$Tumor_Sample_Barcode.ct <- gsub('_ctDNA','', tmp$Tumor_Sample_Barcode.ct)
colnames(tmp)[3] <- 'WXS';
colnames(tmp)[5] <- 'ctDNA';

version2 <- tmp;
