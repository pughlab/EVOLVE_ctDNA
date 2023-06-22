### plot_mutation_summary.R ########################################################################
# Plot somatic mutation profile of ctDNA, indicated expected result (based on prior whole-exome
# sequencing of the matched tumour tissue) where possible

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
#library(GenomicRanges);

### VARIANT CODING
# 1 = missense, 2 = stop gain, 3 = stop loss, 4 = splicing, 5 = frameshift, 6 = in frame indel, 7 = tss
# 8 = RNA, 9 = other (up/downstream, UTR, intergenic, silent, intron), 10 = ITD
variant.codes <- data.frame(
	Classification = c("3'Flank", "5'Flank", "Intron", "RNA", "IGR", "3'UTR", "5'UTR", "Silent",
		"Missense_Mutation", "Splice_Region", "Splice_Site", "In_Frame_Del", "In_Frame_Ins",
		"Frame_Shift_Del", "Frame_Shift_Ins", "Nonsense_Mutation", "Nonstop_Mutation",
		"Translation_Start_Site", "ITD"),
	Group = c('other','other','other','RNA','other','other','other','other','missense',
		'splice_site','splice_site','in_frame_indel','in_frame_indel','frameshift_del',
		'frameshift_ins', 'nonsense', 'nonstop', 'tss', 'itd')
	);

# for these plots, we will ignore some variant types
variant.colours <- c('darkseagreen4','#9AA3F2','yellow','darkorange3','#F9B38E','grey50')
names(variant.colours) <- c('missense','nonsense','splicing','frameshift_ins','frameshift_del','noncoding');

variant.codes$Code <- c(1:6)[match(variant.codes$Group, c('missense','nonsense','splice_site','frameshift_ins','frameshift_del','other'))];

# list samples to exclude
smps.without.wxs <- c('EVO-009-014-Ar','EVO-009-018-Bx','EVO-400-004-Bx');
smps.without.ctdna <- c('EVO-009-025-Ar','EVO-009-025-Bx','EVO-400-001-Ar','EVO-400-001-Bx','EVO-400-002-Ar','EVO-400-002-Bx','EVO-400-006-Ar','EVO-400-006-Bx');

### READ DATA ######################################################################################
# get clinical covariates
load('/Users/sprokopec/git/EVOLVE_ctDNA/data/EVOLVE_ctDNA__clinical_timeline.RData');

# get mutation data
load('/Users/sprokopec/git/EVOLVE_ctDNA/data/EVOLVE_ctDNA__mutation_data.RData');

# get estimated ctDNA levels
tumour.content <- read.delim('/Users/sprokopec/git/EVOLVE_ctDNA/data/estimated_tumour_content.txt')

### FORMAT PLOT DATA ###############################################################################
# apply variant coding
mutation.data$Code <- variant.codes$Code[match(
	mutation.data$Variant_Classification,
	variant.codes$Classification)];

# reshape data for plotting
plot.data <- reshape(
	unique(mutation.data[!is.na(mutation.data$Sample),c('Hugo_Symbol','Sample','Code')]),
	direction = 'wide',
	timevar = 'Sample',
	idvar = 'Hugo_Symbol'
	);

rownames(plot.data) <- plot.data$Hugo_Symbol;
plot.data <- plot.data[,-1];
colnames(plot.data) <- gsub('Code\\.','',colnames(plot.data));

dot.data <- reshape( 
	unique(mutation.data[!is.na(mutation.data$Sample),c('Hugo_Symbol','Sample','Validated')]),
	direction = 'wide',
	timevar = 'Sample',
	idvar = 'Hugo_Symbol'
	);

rownames(dot.data) <- dot.data$Hugo_Symbol;
dot.data <- dot.data[,-1];
colnames(dot.data) <- gsub('Validated\\.','',colnames(dot.data));

# order data
all.patients <- as.character(unique(
	timeline[which(timeline$Sample != ''),]$Patient
	));

phenodata <- timeline[which(timeline$Sample != ''),];
phenodata$Patient <- factor(phenodata$Patient, levels = all.patients);
phenodata$Group <- factor(phenodata$Group, levels = c('baseline','on.trial','EOT'));
phenodata <- phenodata[order(phenodata$Patient, phenodata$Group),];

phenodata$ORDER <- NA;
for (patient in all.patients) {
	idx <- which(phenodata$Patient == patient);
	phenodata[idx,]$ORDER <- 1:length(idx);
	}

# fill in missing samples (no mutation does not mean not tested)
missing.samples <- setdiff(phenodata$Sample, colnames(plot.data));
plot.data[,missing.samples] <- 0;
dot.data[,missing.samples] <- 0;

plot.data[is.na(plot.data)] <- 0;
dot.data[is.na(dot.data)] <- 0;

# fill in missing mutations (in WXS but not ctDNA)
missing.muts <- unique(mutation.data[which(mutation.data$Validated == 1),c('Patient','Hugo_Symbol')]);

for (i in 1:nrow(missing.muts)) {

	gene <- missing.muts[i,]$Hugo_Symbol;
	patient <- missing.muts[i,]$Patient;
	smps <- colnames(dot.data)[grep(patient,colnames(dot.data))];

	if (all(dot.data[gene,smps] == 1)) { next; }
	idx <- which(dot.data[gene,smps] == 0);
	dot.data[gene,smps[idx]] <- 2;

	}

# indicate reversions
for (i in 1:nrow(reversion.data)) {
	patient <- reversion.data[i,]$Patient;
	gene <- reversion.data[i,]$Hugo_Symbol;
	dot.data[gene, grepl(patient, colnames(dot.data))] <- reversion.data[i,]$Dot.Code;
	plot.data[gene, grepl(patient, colnames(plot.data))] <- reversion.data[i,]$BG.Code;
 	}

# indicate CCNE1 amplifications
plot.data['CCNE1',] <- 0;
dot.data['CCNE1',] <- 0;

# based on original EVOLVE publication
wxs.with.ccne1.amps <- c('EVO-009-004','EVO-009-006','EVO-009-007','EVO-009-009','EVO-009-011','EVO-400-007','EVO-400-008');

# use panelCN.mops data
ctdna.with.ccne1.amps <- gsub('_ctDNA','',as.character(cna.data[which(cna.data$Call == 1),]$Sample));

# fill in findings
plot.data['CCNE1',ctdna.with.ccne1.amps] <- 7;
for (patient in wxs.with.ccne1.amps) {
	dot.data['CCNE1',grepl(patient, colnames(dot.data))] <- 2;
	}

dot.data['CCNE1',which(plot.data['CCNE1',] * dot.data['CCNE1',] > 0)] <- 1;

# sort data
gene.order <- c('TP53','BRCA1','BRCA2','PALB2','CCNE1');
plot.data <- plot.data[gene.order,phenodata$Sample];
dot.data <- dot.data[gene.order,phenodata$Sample];

# format purity estimates
tumour.content$Sample <- gsub('_ctDNA','',tumour.content$Sample);

# for cases with no TP53 mutation, use the maximum somatic VAF
tumour.content$Final <- tumour.content$Estimate;
na.idx <- which(is.na(tumour.content$Estimate));
tumour.content[na.idx,]$Final <- tumour.content[na.idx,]$Max.VAF;

purity.estimates <- merge(
	phenodata[,c('Patient','Sample','ORDER')],
	tumour.content[,c('Patient','Sample','Final')]
	);

purity.estimates <- purity.estimates[order(purity.estimates$Patient, purity.estimates$ORDER),];
purity.estimates$Sample <- factor(purity.estimates$Sample, levels = colnames(plot.data));

### COVARIATES #####################################################################################
# make the plot legend (mutation type/consequence)
functional.legend <- legend.grob(
	legends = list(
		legend = list(
			colours = variant.colours[1:4],
			labels = names(variant.colours)[1:4],
			title = 'Variant Type'
			),
		legend = list(
			colours = c(variant.colours[5:6],'red'),
			labels = c(names(variant.colours)[5:6], 'amplification')
			)
		),
	title.just = 'left',
	title.fontface = 'plain',
	label.cex = 0.8,
	title.cex = 0.9,
	layout = c(2,1),
	size = 1.5
	);

# make sample covariates
covariate.data <- merge(phenodata, clinical, by.x = 'Patient', by.y = 'Patient.ID', all.x = TRUE);

smp.covariates <- list(
	rect = list(
		col = 'transparent',
		fill = c('skyblue','seagreen','pink')[match(covariate.data$Group, c('baseline','on.trial','EOT'))],
		lwd = 1
		),
	rect = list(
		col = 'transparent',
		fill = covariate.colours$Response$colours[match(covariate.data$Response.cat, covariate.colours$Response$levels)],
		lwd = 1
		),
	rect = list(
		col = 'transparent',
		fill = covariate.colours$Age$colours[match(covariate.data$Age.cat, covariate.colours$Age$levels)],
		lwd = 1
		),
	rect = list(
		col = 'transparent',
		fill = covariate.colours$BRCA$colours[match(covariate.data$BRCA.cat, covariate.colours$BRCA$levels)],
		lwd = 1
		),
	rect = list(
		col = 'transparent',
		fill = covariate.colours$Tissue$colours[match(covariate.data$Oncotree.Code, covariate.colours$Tissue$levels)],
		lwd = 1
		),
	rect = list(
		col = 'transparent',
		fill = covariate.colours$Cohort$colours[match(covariate.data$Cohort.cat, covariate.colours$Cohort$levels)],
		lwd = 1
		)
	);

smp.covariate.grob <- covariates.grob(
	covariates = smp.covariates,
	ord = 1:nrow(covariate.data),
	side = 'top',
	size = 0.7,
	grid.row = list(col = 'white', lwd = 1),
	row.lines = 1:length(smp.covariates),
	grid.col = list(col = 'white', lwd = 1),
	col.lines = get.line.breaks(covariate.data$Patient)-0.5
	);

smp.legends <- list(
	legend = list(
		colours = covariate.colours$Cohort$colours,
		labels = covariate.colours$Cohort$labels,
		title = 'Cohort'
		),
	legend = list(
		colours = covariate.colours$Tissue$colours,
		labels = covariate.colours$Tissue$labels,
		title = 'Tissue'
		),
	legend = list(
		colours = covariate.colours$BRCA$colours,
		labels = covariate.colours$BRCA$labels,
		title = 'Germline BRCA'
		),
	legend = list(
		colours = covariate.colours$Age$colours,
		labels = covariate.colours$Age$labels,
		title = 'Age'
		),
	legend = list(
		colours = covariate.colours$Response$colours,
		labels = covariate.colours$Response$labels,
		title = 'Best Response'
		),
	legend = list(
		colours = c('skyblue','seagreen','pink'),
		labels = c('baseline','on-trial','end-of-treatment'),
		title = 'Timepoint'
		)
	);

smp.legend.grob <- legend.grob(
	legends = smp.legends,
	label.cex = 0.8,
	title.cex = 0.9,
	title.fontface = 'plain',
	title.just = 'left',
	size = 1.5,
	layout = c(length(smp.legends),1)
	);

covariate.key <- list(
	text = list(lab = c('Timepoint','Response','Age','gBRCA','Tissue','Cohort'),cex = 0.8,adj = 1)
	);

### PLOT DATA ######################################################################################
# create function to determine spot size
spot.size.vaf <- function(x) {
	sizes <- rep(0,length(x));
	sizes[which(x == 1)] <- 1.5;
	sizes[which(x == 2)] <- 1;
	sizes[which(x == 3)] <- 1;
	sizes[which(x == 4)] <- 1.5;
	return(sizes);
	}
spot.colour <- function(x) { 'black'; }
spot.shape <- function(x) {
	shapes <- rep(19, length(x));
	shapes[which(x == 1)] <- 22
	shapes[which(x == 2)] <- 0;
	shapes[which(x == 3)] <- 5;
	shapes[which(x == 4)] <- 18;
	return(shapes);
	}

dot.key <- list(
	points = list(
		pch = c(22,0,18,5),
		cex = c(1.5,1,1.5,1),
		col = 'black',
		fill = 'black'
		),
	text = list(
		lab = c(
			'expected, confirmed in ctDNA',
			'expected, not detected in ctDNA',
			'reversion expected, confirmed in ctDNA',
			'reversion expected, not detected in ctDNA'
			),
		cex = 0.8
		)
	);

# determine where to put lines
patient.splits <- rep('grey90',nrow(covariate.data));
patient.splits[get.line.breaks(covariate.data$Patient)+0.5] <- 'black';

# make the plot!
mutation.plot <- create.dotmap(
	x = dot.data,
	bg.data = plot.data,
	spot.size.function = spot.size.vaf,
	spot.colour.function = spot.colour,
	pch = spot.shape(unlist(dot.data)),
	colour.scheme = c('white',variant.colours,'red'),
	at = seq(-0.5,length(variant.colours)+2,1),
	legend = list(
		inside = list(fun = smp.covariate.grob, x = 0.5, y = -0.05),
		inside = list(fun = smp.legend.grob, x = -0.04, y = -0.54),
		inside = list(fun = functional.legend, x = 0.56, y = -0.54),
		inside = list(fun = draw.key, args = list(dot.key), x = 0.77, y = -0.61),
		inside = list(fun = draw.key, args = list(key = covariate.key), x = -0.06, y = -0.05)
		),
	top.padding = 1,
	bottom.padding = 19,
	left.padding = 1,
	xaxis.lab = rep('', ncol(dot.data)),
	yaxis.lab = c(
		expression(italic('TP53')),
		expression(italic('BRCA1')),
		expression(italic('BRCA2')),
		expression(italic('PALB2')),
		expression(italic('CCNE1'))
		),
	yaxis.cex = 1,
	xaxis.tck = 0,
	yaxis.tck = 0,
	yaxis.fontface = 'plain',
	na.spot.size = 0,
	bg.alpha = 1,
	lwd = 1,
	col.lwd = 1,
	row.lwd = 1,
	col.colour = patient.splits
	);

purity.estimates$Estimate <- log10(purity.estimates$Final*1000);

purity.plot <- create.scatterplot(
	Estimate ~ Sample,
	purity.estimates,
	type = c('h','p'),
	xaxis.lab = rep('',nrow(tumour.content)),
	xlab.label = NULL,
	ylab.label = 'ctDNA level',
	ylab.cex = 1.1,
	xaxis.tck = 0,
	ylimits = c(0,3),
	yat = seq(0,3,1), #log10(c(0.001,0.01,0.1,1)),
	yaxis.lab = c(0.001,0.01,0.1,1),
	yaxis.cex = 1,
	abline.h = c(1,2), # 1% and 10%
	abline.lty = 2,
	abline.lwd = 1,
	abline.col = c('grey50','red'),
	style = 'Nature'
	);

create.multipanelplot(
	plot.objects = list(purity.plot, mutation.plot),
	plot.objects.heights = c(1,2),
	plot.objects.widths = 1,
	layout.height = 2,
	layout.width = 1,
	y.spacing = -2,
	ylab.axis.padding = 0,
	left.legend.padding = 0,
	right.legend.padding = 0,
	top.legend.padding = 0,
	bottom.legend.padding = 4,
	top.padding = -1,
	bottom.padding = 6,
	height = 6,
	width = 15,
	resolution = 200,
	filename = generate.filename('EVOLVE_ctDNA', 'mutation_summary__Figure4','png')
	);

### SAVE SESSION INFO ##############################################################################
save.session.profile(generate.filename('Figure4','SessionProfile','txt'));
