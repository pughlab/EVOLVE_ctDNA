### plot_snv_summary.R ############################################################################
# Identify and plot recurrent SNVs.

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
library(GenomicRanges);

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
names(variant.colours) <- c('missense','nonsense','nonstop','splicing','frameshift_ins','frameshift_del','tss','RNA','other');

# for these plots, we will ignore some variant types
variant.colours <- variant.colours[c(1:6,9)];

### READ DATA ######################################################################################
# get data
input.data <- read.delim(
	'2021-08-11_EVOLVE_ctDNA_ensemble_oncokb_annotated_filtered.tsv',
	stringsAsFactors = FALSE
	);

clinical <- read.delim('sample_info.txt', stringsAsFactors = FALSE);

## For now, only include ALL_UNIQUE
clinical <- droplevels(clinical[which(clinical$BAM == 'allUNIQUE'),]);

clinical[which(clinical$Timepoint == 'EOT'),]$Timepoint <- '35';
clinical$Timepoint <- as.numeric(clinical$Timepoint);

clinical <- clinical[order(clinical$Patient, clinical$Timepoint),];

# collect list of all samples
all.samples <- as.character(clinical$Sample);

# get target regions
target_bed <- read.delim('CHARM-MMR_plus_EVOLVE_hg38.bed', header = FALSE);
#target_bed <- read.delim('/cluster/projects/pughlab/references/intervals/CHARM_EVOLVE_panel/hg38_liftover/CHARM-MMR_plus_EVOLVE_hg38_liftover.bed', header = FALSE);
colnames(target_bed) <- c('Chromosome','Start','End','ID','V5','Strand');

## For now, do not include mutations in MSI and Agena (genotyping) sites
target_genes <- subset(target_bed,!grepl("MSI|Agena",ID) | grepl("HBOC",ID))

### FORMAT DATA ####################################################################################
# because this is tightly targeted, fix these cases
mod.idx <- which(input.data$Hugo_Symbol == 'RUNDC3B' & grepl('ABCB1,intron_variant', input.data$all_effects));
input.data[mod.idx,]$Hugo_Symbol <- 'ABCB1';
input.data[mod.idx,]$Variant_Classification <- 'Intron';

mod.idx <- which(input.data$Hugo_Symbol == 'WRAP53' & grepl('TP53,intron_variant', input.data$all_effects));
input.data[mod.idx,]$Hugo_Symbol <- 'TP53';
input.data[mod.idx,]$Variant_Classification <- 'Intron';

#maf <- maf[-which(grepl('EVO-400-004', maf$Tumor_Sample_Barcode) & maf$Hugo_Symbol == 'BRCA2'),];
 
# indicate key fields
keep.fields <- c('Tumor_Sample_Barcode','Hugo_Symbol','Chromosome','Start_Position','End_Position','Variant_Classification','Reference_Allele','Tumor_Seq_Allele2','t_depth','t_ref_count','Variant_Type');

mutation.data <- input.data[which(input.data$Tumor_Sample_Barcode %in% all.samples), keep.fields];

# calculate tumour vaf
mutation.data$t_vaf <- 1-(mutation.data$t_ref_count/mutation.data$t_depth);

## apply some additional filters
# remove low coverage variants
mutation.data <- mutation.data[which(mutation.data$t_depth >= 100),];

# remove probable germline (high VAF)
mutation.data <- mutation.data[which(mutation.data$t_vaf < 0.8),];

# apply variant coding
mutation.data$Code <- variant.codes$Code[match(mutation.data$Variant_Classification, variant.codes$Classification)];
if (any(mutation.data$Code > 6)) {
	mutation.data[which(mutation.data$Code > 6),]$Code <- 9;
	}

# filter to target regions
target.gr <- makeGRangesFromDataFrame(target_genes, starts.in.df.are.0based = TRUE);
mutation.gr <- makeGRangesFromDataFrame(mutation.data,
	seqnames.field = 'Chromosome',
	start.field = 'Start_Position',
	end.field = 'End_Position',
	keep.extra.columns = TRUE
	);

mutation.data.filtered <- as.data.frame(subsetByOverlaps(mutation.gr, target.gr, ignore.strand = TRUE));

# reduce to 1 mutation per gene per sample [taking the higher priority code]
mutation.data.trimmed <- aggregate(
	Code ~ Tumor_Sample_Barcode + Hugo_Symbol + seqnames,
#	mutation.data.filtered[which(mutation.data.filtered$Variant_Type == 'SNP'),],
	mutation.data.filtered[which(mutation.data.filtered$Variant_Type != 'SNP'),],
	min
	);

mutation.data.trimmed <- merge(
	mutation.data.trimmed,
	mutation.data[,c('Tumor_Sample_Barcode','Hugo_Symbol','Code','t_vaf')],
	all.x = TRUE
	);
colnames(mutation.data.trimmed)[4] <- 'Chromosome';

mutation.data.trimmed <- mutation.data.trimmed[order(mutation.data.trimmed$t_vaf, decreasing = TRUE),];

# reshape data
plot.data <- reshape(
	mutation.data.trimmed[,c('Tumor_Sample_Barcode','Chromosome','Hugo_Symbol','Code')],
	direction = 'wide',
	timevar = 'Tumor_Sample_Barcode',
	idvar = c('Chromosome','Hugo_Symbol')
	);
colnames(plot.data) <- gsub('Code.','',colnames(plot.data));
rownames(plot.data) <- plot.data$Hugo_Symbol;

vaf.data <- reshape(
	mutation.data.trimmed[,c('Tumor_Sample_Barcode','Chromosome','Hugo_Symbol','t_vaf')],
	direction = 'wide',
	timevar = 'Tumor_Sample_Barcode',
	idvar = c('Chromosome','Hugo_Symbol')
	);
colnames(vaf.data) <- gsub('t_vaf.','',colnames(vaf.data));
rownames(vaf.data) <- vaf.data$Hugo_Symbol;

# get per-gene sample counts
missing.samples <- setdiff(all.samples,colnames(plot.data));
if (length(missing.samples) > 0) {
	plot.data[,missing.samples] <- NA;
	vaf.data[,missing.samples] <- NA;
	}
plot.data <- plot.data[,all.samples];

plot.data['PALB2',] <- NA;
vaf.data['PALB2',] <- NA;

plot.data$Count <- apply(plot.data[,all.samples],1,function(i) { length(i[!is.na(i)]) } );
plot.data <- plot.data[order(-plot.data$Count),];

# create heatmap for recurrent genes (ordered by recurrence)
#heatmap.data <- t(plot.data[,all.samples]);
#heatmap.data[!is.na(heatmap.data)] <- 1;
#heatmap.data <- heatmap.data[do.call(order, transform(heatmap.data)),];

#sample.order <- rownames(heatmap.data);
vaf.data <- vaf.data[rownames(plot.data),all.samples];
vaf.data[is.na(vaf.data)] <- 0;

# create function to determine spot size
modifier <- if (length(all.samples) < 12) { 1.5 } else { 1 }
spot.size.vaf <- function(x) { abs(x)* modifier; }
spot.colour.vaf <- function(x) {
	sapply(x, function(i) { if (i >= 0.5) { 'black' } else { 'grey60' } } )
	}

dot.key <- list(
	points = list(
		col = c('black','black','grey60','grey60','grey60'),
		cex = spot.size.vaf(c(0.8,0.5,0.2,0.1,0)),
		pch = 19
		),
	text = list(
		lab = c('80','50','20','10','0'),
		cex = 0.7
		),
	padding.text = 1.05,
	title = 'VAF (%)',
	cex.title = 0.8
	);

### SET-UP COVARIATES ##############################################################################
# set up covariates
covariate.data <- clinical;
covariate.data$Sample <- factor(covariate.data$Sample,levels = all.samples);
covariate.data <- covariate.data[order(covariate.data$Sample),];

type.colours <- rep('seagreen', length(all.samples));
type.colours[which(covariate.data$Group == 'Screening')] <- 'skyblue';
type.colours[which(covariate.data$Group == 'EOT')] <- 'pink';

covariates.obj <- covariates.grob(
	covariates = list(
		rect = list(col = 'transparent', fill = type.colours, lwd = 1)
		),
	ord = 1:nrow(covariate.data),
	side = 'top',
	grid.col = list(col = 'black', lwd = 1),
	col.lines = get.line.breaks(covariate.data$Patient)-0.5,
	grid.border = list(col = 'black', lwd = 1)
	);

# make the plot legend (mutation type/consequence)
functional.legend <- list(
	legend = list(
		colours = c('skyblue','seagreen','pink'),
		labels = c('baseline','on-trial','end-of-trial'),
		title = 'Clinical'
		),
	legend = list(
		colours = variant.colours,
		labels = names(variant.colours),
		title = 'Function'
		)
	);

clinical.legend <- legend.grob(
	legends = functional.legend,
	title.just = 'left',
	label.cex = 0.7,
	size = 1.2
	);

### PLOTTING #######################################################################################
# make the mutation heatmap with VAF dots
create.dotmap(
	x = vaf.data,
	bg.data = plot.data[rownames(vaf.data),colnames(vaf.data)],
	spot.size.function = spot.size.vaf,
	spot.colour.function = spot.colour.vaf,
	pch = 19,
	colour.scheme = variant.colours,
	at = seq(0,length(variant.colours),1),
	total.colours = length(variant.colours)+1,
	colourkey = FALSE,
	legend = list(
		right = list(fun = draw.key, args = list(key = dot.key)),
		inside = list(fun = covariates.obj, x = 0.5, y = 1.1),
		inside = list(fun = clinical.legend, x = 1.02, y = 1.15)
		),
	top.padding = 5,
	right.padding = 5,
	xaxis.lab = gsub('_ctDNA','',gsub('_allUNIQUE', '', colnames(vaf.data))),
	yaxis.lab = rownames(vaf.data),
	xaxis.cex = 0.7,
	yaxis.cex = 1,
	xaxis.tck = 0,
	yaxis.tck = 0,
	xaxis.fontface = 'plain',
	yaxis.fontface = 'plain',
	xaxis.rot = 90,
	bg.alpha = 1,
	lwd = 1,
	col.lwd = c(2,0.5),
	row.lwd = 2,
	col.colour = 'grey80',
	row.colour = 'grey80',
	height = if (nrow(plot.data) > 20) { 8 } else { 6 },
	width = if (length(all.samples) > 30) { 11 } else { 8 },
#	resolution = 200,
	filename = generate.filename('EVOLVE_ctDNA_allUNIQUE', '_indel_landscape','png')
	);

### SAVE SESSION INFO ##############################################################################
save.session.profile(generate.filename('SNVSummary','SessionProfile','txt'));
