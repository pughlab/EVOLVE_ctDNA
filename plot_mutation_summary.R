### plot_validated_snvs.R ##########################################################################
# Plot status of known (exome) mutations detected in ctDNA

### PREPARE SESSION ################################################################################
# import libraries
library(BoutrosLab.plotting.general);
library(GenomicRanges);

source('/cluster/home/sprokope/git/analysis/helper_functions/session.functions.R');

input.dir <- '/cluster/projects/ovgroup/projects/OV_Superset/EVOLVE/ctDNA/pipeline_suite/';
output.dir <- '/cluster/projects/ovgroup/projects/OV_Superset/EVOLVE/ctDNA/pipeline_suite/paper_figures/';

setwd(input.dir);

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
		'frameshift_ins', 'nonsense', 'nonstop', 'tss', 'itd'),
	Code = c(9, 9, 9, 8, 9, 9, 9, 9, 1, 4, 4, 7, 7, 6, 5, 2, 3, 7, 10)
	);

#variant.colours <- c('darkseagreen4','darkorchid4','#9AA3F2','yellow','darkorange3','#F9B38E','turquoise1','plum','grey50')
#names(variant.colours) <- c('missense','nonsense','nonstop','splicing','frameshift_ins','frameshift_del','tss','RNA','noncoding');

# for these plots, we will ignore some variant types
variant.colours <- c('darkseagreen4','#9AA3F2','yellow','darkorange3','#F9B38E','grey50')
names(variant.colours) <- c('missense','nonsense','splicing','frameshift_ins','frameshift_del','noncoding');

variant.codes$NEW.CODE <- c(1:6)[match(variant.codes$Group, c('missense','nonsense','splice_site','frameshift_ins','frameshift_del','other'))];

# list samples to exclude
smps.without.wxs <- c('EVO-009-014-Ar','EVO-009-018-Bx','EVO-400-004-Bx');
smps.without.ctdna <- c('EVO-009-025-Ar','EVO-009-025-Bx','EVO-400-001-Ar','EVO-400-001-Bx','EVO-400-002-Ar','EVO-400-002-Bx','EVO-400-006-Ar','EVO-400-006-Bx');

### READ DATA ######################################################################################
# get target regions
target_bed <- read.delim('CHARM-MMR_plus_EVOLVE_hg38.bed', header = FALSE)
colnames(target_bed) <- c('Chromosome','Start','End','ID','V5','Strand');

# get clinical covariates
load('2022-09-06_EVOLVE_ctDNA_clinicalCovariates.RData');

sample.info <- read.delim('configs/2022-09-06_sample_info_with_batch.txt');

# read in data
timeline <- read.delim(
	'cBioportal/EVOLVE_ctDNA_20220908/data_clinical_timeline_specimen.txt', stringsAsFactors = FALSE);
cbio.data <- read.delim('cBioportal/EVOLVE_ctDNA_20220908/mutation_data_extended.txt', stringsAsFactors = FALSE);

# get cnv data (for CCNE1 amplifications)
cnvkit <- read.delim('CNVKit/cnv_calls_with_purity__ci/2022-09-12_EVOLVE_ctDNA__cnvkit_calls.tsv');

# move to output directory
setwd(output.dir);

### FORMAT PLOT DATA ###############################################################################
cbio.data$Patient <- substr(cbio.data$Tumor_Sample_Barcode, 0, 11);

# get clinical data (for sorting)
timeline <- timeline[order(timeline$PATIENT_ID, timeline$START_DATE),c(1,5,2,7)];
timeline$SAMPLE_ID <- gsub('_ctDNA','', timeline$SAMPLE_ID);
timeline <- timeline[-which(timeline$SAMPLE_ID %in% smps.without.wxs),];
timeline <- timeline[-which(timeline$SAMPLE_ID %in% smps.without.ctdna),];

# indicate samples to keep (some exome do not have ctdna)
keep.exome.smps <- timeline[which(timeline$SPECIMEN_TYPE == 'TISSUE'),]$SAMPLE_ID;

# format target regions
target_bed$Start <- target_bed$Start - 100;
target_bed$End <- target_bed$End + 100;

target.gr <- makeGRangesFromDataFrame(
	target_bed[,1:3],
	ignore.strand = TRUE,
	starts.in.df.are.0based = TRUE
	);

evolve.gr <- makeGRangesFromDataFrame(
	cbio.data[which(cbio.data$Mutation_Status == 'somatic'),],
	seqnames.field = 'Chromosome',
	start.field = 'Start_Position',
	end.field = 'End_Position',
	ignore.strand = TRUE,
	keep.extra.columns = TRUE
	);

# filter data to target regions
cbio.filtered <- data.frame(subsetByOverlaps(evolve.gr, target.gr));
colnames(cbio.filtered)[1:3] <- c('Chromosome','Start_Position','End_Position');

key.fields <- c('Patient','Tumor_Sample_Barcode','Hugo_Symbol','Chromosome','Start_Position','End_Position','Variant_Classification','Variant_Type','HGVSc','HGVSp','HGVSp_Short','t_depth','t_ref_count','t_alt_count','n_depth','n_ref_count','n_alt_count');

# filter a few synonymous/intronic variants
cbio.filtered <- cbio.filtered[which(cbio.filtered$Consequence != 'synonymous_variant'),];
cbio.filtered <- cbio.filtered[which(cbio.filtered$Consequence != 'intron_variant'),];

# merge exome and ctDNA by position
exome.smps <- which(
	!is.na(cbio.filtered$Matched_Norm_Sample_Barcode) &
	cbio.filtered$Tumor_Sample_Barcode %in% keep.exome.smps
	);
ctdna.smps <- which(is.na(cbio.filtered$Matched_Norm_Sample_Barcode));

merged <- merge(
	cbio.filtered[exome.smps,key.fields],
	cbio.filtered[ctdna.smps,key.fields[1:14]],
	by = c('Patient','Hugo_Symbol','Chromosome','Start_Position','End_Position','Variant_Classification','Variant_Type','HGVSc','HGVSp','HGVSp_Short'),
	suffixes = c('.wxs','.ctdna'),
	all = TRUE
	);

# and only targeted genes
merged <- merged[which(merged$Hugo_Symbol %in% c('TP53','BRCA1','BRCA2','PALB2','ABCB1','CCNE1')),];

# calculate VAF
merged$t_vaf.wxs <- merged$t_alt_count.wxs / merged$t_depth.wxs
merged$t_vaf.ctdna <- merged$t_alt_count.ctdna / merged$t_depth.ctdna
merged[which(merged$t_vaf.ctdna < 0.01),]$t_vaf.ctdna <- NA;

# apply variant coding
merged$Code <- variant.codes$NEW.CODE[match(merged$Variant_Classification, variant.codes$Classification)];

# was the WXS mutation validated in ctDNA?
merged$Validated <- 0;
merged[which(merged$t_vaf.wxs > 0),]$Validated <- 1;

# sort data
merged <- merged[order(merged$Patient, merged$Hugo_Symbol, -merged$t_vaf.wxs, -merged$t_vaf.ctdna, na.last = TRUE),];

# remove redundant entries
to.remove <- intersect(
	which(duplicated(merged[,c('Patient','Hugo_Symbol','Tumor_Sample_Barcode.ctdna')])),
	which(is.na(merged$Tumor_Sample_Barcode.wxs))
	);

if (length(to.remove) > 0) {
	merged <- merged[-to.remove,];
	}

# fix sample names
merged$Tumor_Sample_Barcode.ctdna <- gsub('_ctDNA','',merged$Tumor_Sample_Barcode.ctdna);

# drop any ONLY ctDNA variant with VAF < 0.01
if (any(is.na(merged$t_vaf.ctdna) & is.na(merged$t_vaf.wxs))) {
	merged <- merged[-which(is.na(merged$t_vaf.ctdna) & is.na(merged$t_vaf.wxs)),];
	}

## format data for plotting
plot.data <- reshape(
	unique(merged[!is.na(merged$Tumor_Sample_Barcode.ctdna),c('Hugo_Symbol','Tumor_Sample_Barcode.ctdna','Code')]),
	direction = 'wide',
	timevar = 'Tumor_Sample_Barcode.ctdna',
	idvar = 'Hugo_Symbol'
	);

rownames(plot.data) <- plot.data$Hugo_Symbol;
plot.data <- plot.data[,-1];
colnames(plot.data) <- gsub('Code\\.','',colnames(plot.data));

dot.data <- reshape( 
	unique(merged[!is.na(merged$Tumor_Sample_Barcode.ctdna),c('Hugo_Symbol','Tumor_Sample_Barcode.ctdna','Validated')]),
	direction = 'wide',
	timevar = 'Tumor_Sample_Barcode.ctdna',
	idvar = 'Hugo_Symbol'
	);

rownames(dot.data) <- dot.data$Hugo_Symbol;
dot.data <- dot.data[,-1];
colnames(dot.data) <- gsub('Validated\\.','',colnames(dot.data));

# order data
all.patients <- unique(timeline[which(timeline$SPECIMEN_TYPE == 'BLOOD'),]$PATIENT_ID);

phenodata <- merge(
	clinical,
	sample.info,
	by.x = 'Patient.ID',
	by.y = 'Patient',
	all = TRUE
	);

phenodata <- phenodata[which(phenodata$Patient.ID %in% all.patients),];
phenodata$Patient.ID <- factor(phenodata$Patient.ID, levels = unique(timeline$PATIENT_ID));
phenodata$Sample <- gsub('_ctDNA','',phenodata$Sample);
phenodata$Group <- factor(phenodata$Group, levels = c('baseline','on.trial','EOT'));
phenodata <- phenodata[order(phenodata$Patient.ID, phenodata$Group),];

phenodata$ORDER <- NA;
for (patient in all.patients) {
	idx <- which(phenodata$Patient.ID == patient);
	phenodata[idx,]$ORDER <- 1:length(idx);
	}

# fill in missing samples (no mutation does not mean not tested)
missing.samples <- setdiff(phenodata$Sample, colnames(plot.data));
plot.data[,missing.samples] <- 0;
dot.data[,missing.samples] <- 0;

plot.data[is.na(plot.data)] <- 0;
dot.data[is.na(dot.data)] <- 0;

# fill in missing mutations (in WXS but not ctDNA)
missing.muts <- unique(merged[which(merged$Validated == 1),c('Patient','Hugo_Symbol')]);

for (i in 1:nrow(missing.muts)) {

	gene <- missing.muts[i,]$Hugo_Symbol;
	patient <- missing.muts[i,]$Patient;
	smps <- colnames(dot.data)[grep(patient,colnames(dot.data))];

	if (all(dot.data[gene,smps] == 1)) { next; }
	idx <- which(dot.data[gene,smps] == 0);
	dot.data[gene,smps[idx]] <- 2;

	}

# indicate reversions
#dot.data['BRCA1',grepl('EVO-009-001', colnames(dot.data))] <- 2; # somatic snp + somatic DEL
dot.data['BRCA1',grepl('EVO-009-003', colnames(dot.data))] <- 3; # germline snp + somatic snp
dot.data['BRCA1',grepl('EVO-009-013', colnames(dot.data))] <- 3; # germline snp + somatic indel
dot.data['BRCA1',grepl('EVO-009-023', colnames(dot.data))] <- 3; # germline snp + somatic indel

plot.data['BRCA2',grepl('EVO-009-006', colnames(plot.data))] <- 6; # triallelic
dot.data['BRCA2',grepl('EVO-009-006', colnames(dot.data))] <- 4; # triallelic
 
# indicate CCNE1 amplifications
plot.data['CCNE1',] <- 0;
dot.data['CCNE1',] <- 0;

wxs.with.ccne1.amps <- c('EVO-009-002','EVO-009-004','EVO-009-006','EVO-009-007','EVO-009-008','EVO-009-009','EVO-009-011');
#ctdna.with.ccne1.amps <- c('EVO-009-004_C2D1','EVO-009-007_C2D1','EVO-009-013_Screening','EVO-009-013_C2D1','EVO-009-017_C2D1','EVO-009-020_C2D1','EVO-009-023_C2D1','EVO-009-024_C2D1','EVO-400-003_C2D1');

# use CNVKit data
ccne1.amps <- cnvkit[which(cnvkit$gene == 'CCNE1' & cnvkit$cn >= 4),];
ctdna.with.ccne1.amps <- gsub('_ctDNA','',as.character(ccne1.amps$Sample));

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

# fix coding (for heatmap colours)
#plot.data[which(plot.data == 9, arr.ind = TRUE)] <- 8;

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
covariate.data <- phenodata;

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
	col.lines = get.line.breaks(covariate.data$Patient.ID)-0.5
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
patient.splits[get.line.breaks(covariate.data$Patient.ID)+0.5] <- 'black';

# make the plot!
create.dotmap(
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
	col.colour = patient.splits,
	height = 4.4,
	width = 15,
	resolution = 1200,
	filename = generate.filename('EVOLVE_ctDNA', 'SNV_validation_heatmap','png')
	);

### SAVE SESSION INFO ##############################################################################
save.session.profile(generate.filename('SNVplot','SessionProfile','txt'));
