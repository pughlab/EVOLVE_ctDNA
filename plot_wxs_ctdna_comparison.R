### plot_validated_snvs.R ##########################################################################
# Plot status of known (exome) mutations detected in ctDNA

### PREPARE SESSION ################################################################################
# import libraries
library(BoutrosLab.plotting.general);
library(GenomicRanges);

source('/cluster/home/sprokope/git/analysis/helper_functions/session.functions.R');

input.dir <- '/cluster/projects/ovgroup/projects/OV_Superset/EVOLVE/ctDNA/pipeline_suite/';
output.dir <- '/cluster/projects/ovgroup/projects/OV_Superset/EVOLVE/ctDNA/pipeline_suite/Ensemble_calls/';

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

variant.colours <- c('darkseagreen4','darkorchid4','#9AA3F2','yellow','darkorange3','#F9B38E','turquoise1','plum','grey50')
names(variant.colours) <- c('missense','nonsense','nonstop','splicing','frameshift_ins','frameshift_del','tss','RNA','noncoding');

# for these plots, we will ignore some variant types
variant.colours <- variant.colours[c(1:7,9)];

# list samples to exclude
smps.without.wxs <- c('EVO-009-014-Ar','EVO-009-018-Bx','EVO-400-004-Bx');
smps.without.ctdna <- c('EVO-009-025-Ar','EVO-009-025-Bx','EVO-400-001-Ar','EVO-400-001-Bx','EVO-400-002-Ar','EVO-400-002-Bx','EVO-400-006-Ar','EVO-400-006-Bx');

### READ DATA ######################################################################################
# get target regions
target_bed <- read.delim('CHARM-MMR_plus_EVOLVE_hg38.bed', header = FALSE)
colnames(target_bed) <- c('Chromosome','Start','End','ID','V5','Strand');

# get clinical covariates
load('2022-06-28_EVOLVE_ctDNA_clinicalCovariates.RData');

# read in data
timeline <- read.delim('cBioportal/EVOLVE_ctDNA_20220726/data_clinical_timeline_specimen.txt');
cna.data <- read.delim('cBioportal/EVOLVE_ctDNA_20220726/discrete_cna_data.txt');
cbio.data <- read.delim('cBioportal/EVOLVE_ctDNA_20220726/mutation_data_extended.txt');

setwd(output.dir);

ccne1.data <- read.delim('');

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

# calculate VAF
merged$t_vaf.wxs <- merged$t_alt_count.wxs / merged$t_depth.wxs
merged$t_vaf.ctdna <- merged$t_alt_count.ctdna / merged$t_depth.ctdna

###
df <- merged[,c('Hugo_Symbol','HGVSp_Short','Tumor_Sample_Barcode.wxs','Tumor_Sample_Barcode.ctdna','t_vaf.wxs','t_vaf.ctdna')];
colnames(df) <- c('Gene','HGVSp','WXS','ctDNA','VAF.wxs','VAF.ct')

df$ctDNA <- gsub('_ctDNA','',df$ctDNA)
df$VAF.wxs <- round(df$VAF.wxs,3)
df$VAF.ct <- round(df$VAF.ct,3)
###

# apply variant coding
merged$Code <- variant.codes$Code[match(merged$Variant_Classification, variant.codes$Classification)];

# for the purposes of this plot, only include WXS mutations
merged <- merged[!is.na(merged$Tumor_Sample_Barcode.wxs),];
merged <- merged[order(merged$Patient, merged$Hugo_Symbol, -merged$t_vaf.wxs),];

to.remove <- intersect(
	which(duplicated(merged[,c('Patient','Hugo_Symbol','Tumor_Sample_Barcode.wxs')])),
	which(is.na(merged$Tumor_Sample_Barcode.ctdna))
	);

if (length(to.remove) > 0) {
	merged <- merged[-to.remove,];
	}

# and only targeted genes
merged <- merged[which(merged$Hugo_Symbol %in% c('TP53','BRCA1','CCNE1')),];

# fix sample names
merged$Tumor_Sample_Barcode.ctdna <- gsub('_ctDNA','',merged$Tumor_Sample_Barcode.ctdna);

# make the plot legend (mutation type/consequence)
functional.legend <- legend.grob(
	legends = list(
		legend = list(
			colours = variant.colours[1:4],
			labels = names(variant.colours)[1:4],
			title = 'Variant Type'
			),
		legend = list(
			colours = variant.colours[5:8],
			labels = names(variant.colours)[5:8]
			),
		legend = list(
			colours = c('red','white','blue'),
			labels = c('amplification','neutral','deep deletion')
			)
		),
	title.just = 'left',
	title.fontface = 'plain',
	label.cex = 0.7,
	title.cex = 0.7,
	layout = c(3,1),
	size = 1
	);

# format plot data
plot.data <- merge(
	reshape(
		unique(merged[,c('Hugo_Symbol','Tumor_Sample_Barcode.wxs','Code')]),
		direction = 'wide', timevar = 'Tumor_Sample_Barcode.wxs', idvar = 'Hugo_Symbol'),
	reshape(
		unique(merged[!is.na(merged$Tumor_Sample_Barcode.ctdna),c('Hugo_Symbol','Tumor_Sample_Barcode.ctdna','Code')]),
		direction = 'wide', timevar = 'Tumor_Sample_Barcode.ctdna', idvar = 'Hugo_Symbol'),
	by = 'Hugo_Symbol',
	all = TRUE
	);

rownames(plot.data) <- plot.data$Hugo_Symbol;
plot.data <- plot.data[,-1];
colnames(plot.data) <- gsub('Code\\.','', colnames(plot.data));

vaf.data <- merge(
	reshape(
		unique(merged[,c('Hugo_Symbol','Tumor_Sample_Barcode.wxs','t_vaf.wxs')]),
		direction = 'wide', timevar = 'Tumor_Sample_Barcode.wxs', idvar = 'Hugo_Symbol'),
	reshape(
		unique(merged[!is.na(merged$Tumor_Sample_Barcode.ctdna),c('Hugo_Symbol','Tumor_Sample_Barcode.ctdna','t_vaf.ctdna')]),
		direction = 'wide', timevar = 'Tumor_Sample_Barcode.ctdna', idvar = 'Hugo_Symbol'),
	by = 'Hugo_Symbol',
	all = TRUE
	);

rownames(vaf.data) <- vaf.data$Hugo_Symbol;
vaf.data <- vaf.data[,-1];
colnames(vaf.data) <- gsub('t_vaf.ctdna.','', gsub('t_vaf.wxs.','',colnames(vaf.data)));

#plot.data['BRCA2',] <- 0;
#plot.data['PALB2',] <- 0;
plot.data['CCNE1',] <- 0;
#vaf.data['BRCA2',] <- 0;
#vaf.data['PALB2',] <- 0;
vaf.data['CCNE1',] <- 0;

# order data
all.patients <- unique(timeline[which(timeline$SPECIMEN_TYPE == 'BLOOD'),]$PATIENT_ID);

clinical <- timeline[which(timeline$PATIENT_ID %in% all.patients),];
clinical$SPECIMEN_TYPE <- factor(clinical$SPECIMEN_TYPE, levels = c('TISSUE','BLOOD'));
clinical <- clinical[order(clinical$PATIENT_ID, clinical$SPECIMEN_TYPE, clinical$START_DATE),];

clinical$ORDER <- NA
for (patient in all.patients) {
	idx <- which(clinical$PATIENT_ID == patient);
	clinical[idx,]$ORDER <- 1:length(idx);
	}

# fill in missing samples (no mutation does not mean not tested)
missing.samples <- setdiff(clinical$SAMPLE_ID, colnames(vaf.data));
vaf.data[,missing.samples] <- 0;
plot.data[,missing.samples] <- 0;

vaf.data <- vaf.data[c('TP53','BRCA1','CCNE1'),clinical$SAMPLE_ID];
plot.data <- plot.data[c('TP53','BRCA1','CCNE1'),clinical$SAMPLE_ID];
vaf.data[is.na(vaf.data)] <- 0;
plot.data[is.na(plot.data)] <- 0;

plot.data[which(plot.data == 9, arr.ind = TRUE)] <- 8;

# add in CN status
gene.cn <- cna.data[which(cna.data$Hugo_Symbol == 'CCNE1'),];
rownames(gene.cn) <- gene.cn$Hugo_Symbol;
colnames(gene.cn) <- gsub('\\.','-',gsub('_ctDNA','', colnames(gene.cn)));
gene.cn <- gene.cn[,colnames(plot.data)];
gene.cn['TP53',] <- 0;
gene.cn['BRCA1',] <- 0;

gene.cn[which(gene.cn == 1, arr.ind = TRUE)] <- 0;
gene.cn[which(gene.cn == -1, arr.ind = TRUE)] <- 0;

# recode CNs
gene.cn[which(gene.cn == -2, arr.ind = TRUE)] <- 9;
gene.cn[which(gene.cn == 2, arr.ind = TRUE)] <- 10;

top.data <- gene.cn;
mid.data <- plot.data;
bottom.data <- gene.cn;

rownames(top.data) <- paste0(rownames(top.data),'.1');
rownames(mid.data) <- paste0(rownames(mid.data),'.2');
rownames(bottom.data) <- paste0(rownames(bottom.data),'.3');

new.bg.data <- rbind(top.data, mid.data, bottom.data);

gene.order <- data.frame(
	ID = rownames(new.bg.data),
	Gene = sapply(rownames(new.bg.data), function(i) { unlist(strsplit(i,'\\.'))[1] } )
	);
gene.order$Gene <- factor(gene.order$Gene, levels = c('TP53','BRCA1','CCNE1'));

gene.order <- gene.order[order(gene.order$Gene, gene.order$ID),];
gene.order$ID <- as.character(gene.order$ID);

top.data[which(top.data != 0, arr.ind = TRUE)] <- 0;
mid.data <- vaf.data;
rownames(mid.data) <- paste0(rownames(mid.data),'.2');
bottom.data[which(bottom.data != 0, arr.ind = TRUE)] <- 0;

new.vaf.data <- rbind(top.data, mid.data, bottom.data);

# fill in gaps
for (gene in levels(gene.order$Gene)) {
	for (i in colnames(new.bg.data)) {
		cn <- new.bg.data[paste0(gene, '.1'),i];
		snv <- new.bg.data[paste0(gene, '.2'),i];
		if (snv == 0) { new.bg.data[paste0(gene, '.2'),i] <- cn; }
		}
	}

### COVARIATES #####################################################################################
# make sample covariates
covariate.data <- merge(clinical, phenodata, by.x = 'PATIENT_ID', by.y = 'Patient.ID');
covariate.data$SAMPLE_ID <- factor(covariate.data$SAMPLE_ID, levels = clinical$SAMPLE_ID);
covariate.data <- covariate.data[order(covariate.data$SAMPLE_ID),];

covariate.data$Time <- NA;
covariate.data[grepl('Ar', covariate.data$SAMPLE_ID),]$Time <- 1;
covariate.data[grepl('Bx|Screening|C1', covariate.data$SAMPLE_ID),]$Time <- 2;
covariate.data[grepl('Pr|EOT', covariate.data$SAMPLE_ID),]$Time <- 4;
for (patient in unique(covariate.data$PATIENT_ID)) {
	idx <- which(covariate.data$PATIENT_ID == patient);
	max.idx <- which(covariate.data$ORDER == max(covariate.data[idx,]$ORDER));
	covariate.data[intersect(idx, max.idx),]$Time <- 4;
	}
covariate.data[grepl('C2D1', covariate.data$SAMPLE_ID),]$Time <- 3;

smp.covariates <- list(
	rect = list(
		col = 'transparent',
		fill = default.colours(6,'seq.blue')[c(1,3,4,6)][match(covariate.data$Time, 1:4)],
		lwd = 1
		),
	rect = list(
		col = 'transparent',
		fill = c('tan','rosybrown')[match(covariate.data$SPECIMEN_TYPE, c('BLOOD','TISSUE'))],
		lwd = 1
		),
	rect = list(
		col = 'transparent',
		fill = covariate.colours$Response$colours[match(covariate.data$Response, covariate.colours$Response$levels)],
		lwd = 1
		),
	rect = list(
		col = 'transparent',
		fill = covariate.colours$Age$colours[match(covariate.data$Age.cat, covariate.colours$Age$levels)],
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
	size = 0.6,
	grid.row = list(col = 'white', lwd = 1),
	row.lines = 1:length(smp.covariates),
	grid.col = list(col = 'white', lwd = 1),
	col.lines = get.line.breaks(covariate.data$PATIENT_ID)-0.5
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
		colours = c('tan','rosybrown'),
		labels = c('plasma','tissue'),
		title = 'Sample Type'
		),
	legend = list(
		colours = default.colours(6, 'seq.blue')[c(1,3,4,6)],
		labels = c('diagnosis','baseline','cycle 2','time of progression'),
		title = 'Timepoint'
		)
	);

smp.legend.grob <- legend.grob(
	legends = smp.legends,
	label.cex = 0.7,
	title.cex = 0.7,
	title.fontface = 'plain',
	title.just = 'left',
	size = 1,
	layout = c(length(smp.legends),1)
	);

### PLOT DATA ######################################################################################
# create function to determine spot size
spot.size.vaf <- function(x) { abs(x); }
spot.colour <- function(x) { 'black'; }

# create a key to describe VAF
dot.key <- list(
	points = list(
		col = 'black',
		cex = spot.size.vaf(c(1,0.8,0.5,0.2,0.1,0)),
		pch = 19
		),
	text = list(
		lab = c('1.0','0.8','0.5','0.2','0.1','0'),
		cex = 0.7
		),
	padding.text = 1.05,
	title = 'VAF',
	cex.title = 0.8
	);

# determine where to put lines
patient.splits <- rep('grey90',nrow(clinical));
patient.splits[get.line.breaks(clinical$PATIENT_ID)+0.5] <- 'black';

gene.splits <- rep('transparent', nrow(gene.order));
gene.splits[get.line.breaks(gene.order$Gene)+0.5] <- 'black';

# make the plot!
create.dotmap(
	x = new.vaf.data[gene.order$ID,],
	bg.data = new.bg.data[gene.order$ID,colnames(new.vaf.data)],
	spot.size.function = spot.size.vaf,
	spot.colour.function = spot.colour,
	pch = 19,
	colour.scheme = c('white',variant.colours,'blue','red'),
	at = seq(-0.5,length(variant.colours)+3,1),
	total.colours = length(variant.colours)+4,
	colourkey = FALSE,
	legend = list(
		inside = list(fun = smp.covariate.grob, x = 0.5, y = -0.05),
		inside = list(fun = draw.key, args = list(key = dot.key), x = -0.06, y = -0.4),
		inside = list(fun = functional.legend, x = 0, y = -0.6),
		inside = list(fun = smp.legend.grob, x = 0.5, y = -0.6)
		),
	top.padding = 1,
	bottom.padding = 15,
	xaxis.lab = rep('', ncol(vaf.data)),
	yaxis.lab = c(
		'', expression(italic('TP53')), '',
		'', expression(italic('BRCA1')), '',
		'', expression(italic('CCNE1')), ''
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
	row.colour = gene.splits,
	height = 3.5,
	width = 15,
	resolution = 200,
	filename = generate.filename('EVOLVE_ctDNA', 'SNV_validation_heatmap','png')
	);

### SAVE SESSION INFO ##############################################################################
save.session.profile(generate.filename('SNVplot','SessionProfile','txt'));
