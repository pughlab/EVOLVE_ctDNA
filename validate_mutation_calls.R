### validate_mutation_calls.R ######################################################################
# Use combined callset to evaluate sensitivity per-tool and for the ensemble method

### PREAMBLE #######################################################################################
library(BoutrosLab.plotting.general);
library(GenomicRanges);

source('/cluster/home/sprokope/git/analysis/helper_functions/session.functions.R');

setwd('/cluster/projects/ovgroup/projects/OV_Superset/EVOLVE/ctDNA/pipeline_suite/Ensemble_calls');

### READ DATA ######################################################################################
# load timeline info
timeline <- read.delim('../configs/2022-09-06_sample_info_with_batch.txt', stringsAsFactors = F);

# load in combined mutation calls
load('2022-07-13_EVOLVE_ctDNA_CombinedMutationData.RData');

# get exome mutation calls (considered true positives)
exome <- read.delim('/cluster/projects/pughlab/projects/OV_Superset/cBioportal/EVOLVE/mutation_data_extended.txt');

# get targeted panel regions
target_bed <- read.delim('../CHARM-MMR_plus_EVOLVE_hg38.bed', header = F)
colnames(target_bed) <- c('Chromosome','Start','End','ID','V5','Strand');

# move to output directory
setwd('wxs_comparison');

### FORMAT DATA ####################################################################################
# format target regions
target_bed$Start <- target_bed$Start - 100;
target_bed$End <- target_bed$End + 100;

target.gr <- makeGRangesFromDataFrame(
	target_bed[,1:3],
	ignore.strand = TRUE,
	starts.in.df.are.0based = TRUE
	);

# extract sample IDs
combined.data$Sample <- combined.data$Tumor_Sample_Barcode;
combined.data$Sample <- gsub('_allUNIQUE','',gsub('_DCS','', gsub('_SSCS','', gsub('_DCS_SSCS','',combined.data$Sample))));

combined.data <- combined.data[!grepl('Screening2',combined.data$Tumor_Sample_Barcode),];

# indicate ensemble method
combined.data$ENSEMBLE <- 0;
combined.data$ENSEMBLE.mod <- 0;

# current method takes ALL MuTect2 calls
combined.data[which(combined.data$FILTER == 'PASS'),]$ENSEMBLE.mod <- 1;

# original method takes MuTect2 + low VAF
combined.data[which(combined.data$Count >= 3),]$ENSEMBLE <- 1;
combined.data[which(combined.data$MuTect2 == 1 & combined.data$VAF < 0.1),]$ENSEMBLE <- 1;

passed.variants <- unique(combined.data[which(combined.data$ENSEMBLE == 1),c('Chromosome','Start_Position','End_Position','Reference_Allele','Tumor_Seq_Allele2','Variant_Type','Patient')]);
passed.variants$C4 <- 1;

combined.data <- merge(
	combined.data[,!grepl('C4',colnames(combined.data))],
	passed.variants,
	all.x = TRUE
	);

combined.data[which(combined.data$C4 == 1),]$ENSEMBLE <- 1;

# format exome mutation calls
exome$Patient <- substr(exome$Tumor_Sample_Barcode,0,11);
germline.keep <- exome[which(exome$Mutation_Status == 'germline'),];
exome <- exome[which(exome$Mutation_Status == 'somatic'),];

t_vaf <- exome$t_alt_count / exome$t_depth;
exome <- exome[which(t_vaf >= 0.05),];
to.remove <- setdiff(
	which(exome$FLAG.high_pop == TRUE),
	which(exome$VARIANT_IN_ONCOKB == 'True')
	);
exome <- exome[-to.remove,];

# remove patients without ctdna data
smps.to.remove <- c('EVO-009-025','EVO-400-001','EVO-400-002','EVO-400-006');

exome <- exome[-which(exome$Patient %in% smps.to.remove),];

# limit exome calls to target regions
exome.gr <- makeGRangesFromDataFrame(
	exome[,c('Chromosome','Start_Position','End_Position','Variant_Type','Patient','Hugo_Symbol')],
	seqnames.field = 'Chromosome',
	start.field = 'Start_Position',
	end.field = 'End_Position',
	ignore.strand = TRUE,
	keep.extra.columns = TRUE
	);

keep.exome <- data.frame(subsetByOverlaps(exome.gr, target.gr));
keep.exome <- keep.exome[,c(1:3,6:8)];
colnames(keep.exome)[1:3] <- c('Chromosome','Start_Position','End_Position');

# limit ctdna calls to target regions (just in case)
ctdna.gr <- makeGRangesFromDataFrame(
	combined.data,
	seqnames.field = 'Chromosome',
	start.field = 'Start_Position',
	end.field = 'End_Position',
	ignore.strand = TRUE,
	keep.extra.columns = TRUE
	);

keep.ctdna <- data.frame(subsetByOverlaps(ctdna.gr, target.gr));
keep.ctdna <- keep.ctdna[,c(1:3,6:10,23,12,14:21,24:25)];

# combine datasets
merged <- merge(
	unique(keep.exome),
	keep.ctdna,
	by.x = c('Chromosome','Start_Position','End_Position','Variant_Type','Hugo_Symbol','Patient'),
	by.y = c('seqnames','start','end','Variant_Type','Hugo_Symbol','Patient'),
	all.x = TRUE
	);

### STATISTICS #####################################################################################
# get sensitivity results
sensitivity.results <- list();
toolset <- c('MuTect','MuTect2','Strelka','VarScan','VarDict','Pindel');

for (tool in toolset) {
	merged[is.na(merged[,tool]),tool] <- 0;
	}

# pindel only calls INDELs
merged[which(merged$Variant_Type == 'SNP'),]$Pindel <- NA;

# mutect (v1) only calls SNPs
merged[which(merged$Variant_Type != 'SNP'),]$MuTect <- NA;


### PER-BAM ###
all.bams <- unique(combined.data$Tumor_Sample_Barcode);
all.bams <- all.bams[!grepl('Screening2',all.bams)];

results <- data.frame(
	BAM = sort(all.bams),
	ENSEMBLE = NA,
	ENSEMBLE.mod = NA
	);
results[,toolset] <- NA;

for (bam in all.bams) {

	patient <- substr(bam,0,11);
	tmp <- merged[which(merged$Patient == patient),];
	total <- nrow(unique(tmp[,1:6]));

	if (total == 0) { next; }
	if (total > 0) { results[which(results$BAM == bam),2:ncol(results)] <- 0; }	

	for (tool in toolset) {
		bam.count <- nrow(tmp[which(tmp[,tool] == 1 & tmp$Tumor_Sample_Barcode == bam),]);

		if ('MuTect' == tool) {
			snp.total <- nrow(unique(tmp[which(tmp$Variant_Type == 'SNP'),1:6]));
			results[which(results$BAM == bam),tool] <- bam.count / snp.total;
			}
		else if ('Pindel' == tool) {
			indel.total <- nrow(unique(tmp[which(tmp$Variant_Type != 'SNP'),1:6]));
			results[which(results$BAM == bam),tool] <- bam.count / indel.total;
			}
		else {
			results[which(results$BAM == bam),tool] <- bam.count / total;
			}
		}

	ensembl.count <- nrow(tmp[which(tmp$Tumor_Sample_Barcode == bam & tmp$ENSEMBLE == 1),]);
	results[which(results$BAM == bam),]$ENSEMBLE <- ensembl.count / total;

	ensembl.count <- nrow(tmp[which(tmp$Tumor_Sample_Barcode == bam & tmp$ENSEMBLE.mod == 1),]);
	results[which(results$BAM == bam),]$ENSEMBLE.mod <- ensembl.count / total;
	
	}

sensitivity.results$per.bam <- results;

### PER-SAMPLE ###
all.samples <- unique(combined.data$Sample);
all.samples <- all.samples[!grepl('Screening2',all.samples)];

results <- data.frame(
	SAMPLE = sort(all.samples),
	ENSEMBLE = NA,
	ENSEMBLE.mod = NA
	);
results[,toolset] <- NA;

for (smp in all.samples) {

	patient <- substr(smp,0,11);
	tmp <- merged[which(merged$Patient == patient),];
	total <- nrow(unique(tmp[,1:6]));

	if (total == 0) { next; }
	if (total > 0) { results[which(results$SAMPLE == smp),2:ncol(results)] <- 0; }	

	for (tool in toolset) {
		smp.count <- nrow(unique(tmp[which(tmp[,tool] == 1 & tmp$Sample == smp),1:6]));

		if ('MuTect' == tool) {
			snp.total <- nrow(unique(tmp[which(tmp$Variant_Type == 'SNP'),1:6]));
			results[which(results$SAMPLE == smp),tool] <- smp.count / snp.total;
			}
		else if ('Pindel' == tool) {
			indel.total <- nrow(unique(tmp[which(tmp$Variant_Type != 'SNP'),1:6]));
			results[which(results$SAMPLE == smp),tool] <- smp.count / indel.total;
			}
		else {
			results[which(results$SAMPLE == smp),tool] <- smp.count / total;
			}
		}

	ensembl.count <- nrow(unique(tmp[which(tmp$Sample == smp & tmp$ENSEMBLE == 1),1:6]));
	results[which(results$SAMPLE == smp),]$ENSEMBLE <- ensembl.count / total;

	ensembl.count <- nrow(unique(tmp[which(tmp$Sample == smp & tmp$ENSEMBLE.mod == 1),1:6]));
	results[which(results$SAMPLE == smp),]$ENSEMBLE.mod <- ensembl.count / total;
	
	}

sensitivity.results$per.sample <- results;

### PER-PATIENT ###
all.patients <- unique(combined.data$Patient);

results <- data.frame(
	PATIENT = sort(all.patients),
	ENSEMBLE = NA,
	ENSEMBLE.mod = NA
	);
results[,toolset] <- NA;

for (patient in all.patients) {

	tmp <- merged[which(merged$Patient == patient),];
	total <- nrow(unique(tmp[,1:6]));

	if (total == 0) { next; }
	if (total > 0) { results[which(results$PATIENT == patient),2:ncol(results)] <- 0; }	

	for (tool in toolset) {
		pat.count <- nrow(unique(tmp[which(tmp[,tool] == 1 & tmp$Patient == patient),1:6]));

		if ('MuTect' == tool) {
			snp.total <- nrow(unique(tmp[which(tmp$Variant_Type == 'SNP'),1:6]));
			results[which(results$PATIENT == patient),tool] <- pat.count / snp.total;
			}
		else if ('Pindel' == tool) {
			indel.total <- nrow(unique(tmp[which(tmp$Variant_Type != 'SNP'),1:6]));
			results[which(results$PATIENT == patient),tool] <- pat.count / indel.total;
			}
		else {
			results[which(results$PATIENT == patient),tool] <- pat.count / total;
			}
		}

	ensembl.count <- nrow(unique(tmp[which(tmp$Patient == patient & tmp$ENSEMBLE == 1),1:6]));
	results[which(results$PATIENT == patient),]$ENSEMBLE <- ensembl.count / total;

	ensembl.count <- nrow(unique(tmp[which(tmp$Patient == patient & tmp$ENSEMBLE.mod == 1),1:6]));
	results[which(results$PATIENT == patient),]$ENSEMBLE.mod <- ensembl.count / total;
	
	}

sensitivity.results$per.patient <- results;

### VISUALIZE DATA #################################################################################
# make some boxplots
plot.objects <- list();

plot.order <- names(sort(apply(sensitivity.results$per.bam[,-1],2,mean,na.rm = T)));
plot.layout.skip <- rep(TRUE, 15);

# per-BAM
plot.data <- reshape(sensitivity.results$per.bam, direction = 'long', varying = list(2:9), v.names = 'Score', timevar = 'Tool', times = colnames(sensitivity.results$per.bam)[2:9]);
plot.data$Tool <- factor(plot.data$Tool, levels = plot.order);

plot.objects[[1]] <- create.boxplot(
	Score ~ Tool | 'BAM',
	plot.data,
	add.stripplot = TRUE,
	points.alpha = 0.5,
	add.text = TRUE,
	text.labels = rep('\u2015',7),
	text.x = 1:8,
	text.y = sort(apply(sensitivity.results$per.bam[,-1],2,mean,na.rm = T)),
	text.cex = 1,
	text.col = 'red'
	);

plot.objects[[2]] <- create.boxplot(
	Score ~ Tool | 'allUNIQUE',
	plot.data[grepl('allUNIQUE', plot.data$BAM),]
	);
plot.objects[[3]] <- create.boxplot(
	Score ~ Tool | 'DCS',
	plot.data[grepl('DCS', plot.data$BAM) & !grepl('SSCS', plot.data$BAM),]
	);
plot.objects[[4]] <- create.boxplot(
	Score ~ Tool | 'DCS+SSCS',
	plot.data[grepl('DCS', plot.data$BAM) & grepl('SSCS', plot.data$BAM),]
	);
plot.objects[[5]] <- create.boxplot(
	Score ~ Tool | 'SSCS',
	plot.data[!grepl('DCS', plot.data$BAM) & grepl('SSCS', plot.data$BAM),]
	);

plot.layout.skip[1:5] <- FALSE;

# per-sample
plot.data <- reshape(sensitivity.results$per.sample, direction = 'long', varying = list(2:9), v.names = 'Score', timevar = 'Tool', times = colnames(sensitivity.results$per.sample)[2:9]);
plot.data$Tool <- factor(plot.data$Tool, levels = plot.order);

plot.objects[[6]] <- create.boxplot(
	Score ~ Tool | 'Sample',
	plot.data,
	add.stripplot = TRUE,
	points.alpha = 0.5,
	add.text = TRUE,
	text.labels = rep('\u2015',7),
	text.x = 1:8,
	text.y = sort(apply(sensitivity.results$per.sample[,-1],2,mean,na.rm = T)),
	text.cex = 1,
	text.col = 'red'
	);

plot.objects[[7]] <- create.boxplot(
	Score ~ Tool | 'baseline',
	plot.data[which(plot.data$SAMPLE %in% timeline[which(timeline$Group == 'baseline'),]$Sample),]
	);
plot.objects[[8]] <- create.boxplot(
	Score ~ Tool | 'cycle2',
	plot.data[which(plot.data$SAMPLE %in% timeline[which(timeline$Group == 'on.trial'),]$Sample),]
	);
plot.objects[[9]] <- create.boxplot(
	Score ~ Tool | 'EOT',
	plot.data[which(plot.data$SAMPLE %in% timeline[which(timeline$Group == 'EOT'),]$Sample),]
	);

plot.layout.skip[6:9] <- FALSE;


# per-patient
plot.data <- reshape(sensitivity.results$per.patient, direction = 'long', varying = list(2:9), v.names = 'Score', timevar = 'Tool', times = colnames(sensitivity.results$per.patient)[2:9]);
plot.data$Tool <- factor(plot.data$Tool, levels = plot.order);

plot.objects[[10]] <- create.boxplot(
	Score ~ Tool | 'Patient',
	plot.data,
	add.stripplot = TRUE,
	points.alpha = 0.5,
	add.text = TRUE,
	text.labels = rep('\u2015',7),
	text.x = 1:8,
	text.y = sort(apply(sensitivity.results$per.patient[,-1],2,mean,na.rm = T)),
	text.cex = 1,
	text.col = 'red'
	);

plot.layout.skip[11] <- FALSE;

# combine them!
create.multiplot(
	plot.objects,
	plot.layout = c(5,3),
	layout.skip = plot.layout.skip,
	filename = generate.filename('EVOLVE_ctDNA','_sensitivity_scores','png'),
	resolution = 200,
	height = 8,
	width = 11,
	ylab.label = 'Sensitivity',
	ylab.cex = 1.5,
	xaxis.cex = 1,
	yaxis.cex = 1,
	xaxis.fontface = 'plain',
	yaxis.fontface = 'plain',
	xaxis.tck = c(0.5,0),
	yaxis.tck = c(0.5,0),
	xaxis.rot = 90,
	axes.lwd = 1
	);

# save data
save(
	sensitivity.results,
	file = generate.filename('EVOLVE_ctDNA','_sensitivity_scores','RData')
	);

# extract key data (all exome mutations and any matching ctDNA calls)
validation <- merged[,c(1:6,10,9,11:13,16,17,15,14,19,20)];
colnames(validation)[which(colnames(validation) == 'VAF')] <- 'M2.vaf';

validation[is.na(validation$ENSEMBLE),]$ENSEMBLE <- 0;
validation[is.na(validation$ENSEMBLE.mod),]$ENSEMBLE.mod <- 0;

write.table(
	validation,
	file = generate.filename('EVOLVE_ctDNA','_validation_data','tsv'),
	row.names = FALSE,
	col.names = TRUE,
	sep = '\t'
	);

### EXAMINE OVERLAP ################################################################################
# visualize overlap (per timepoint + overall)
overlap.data <- data.frame(
	unique(validation[,c('Patient','Chromosome','Start_Position','End_Position','Variant_Type','Hugo_Symbol')]),
	exome = 1,
	baseline = 0,
	cycle2 = 0,
	eot = 0,
	ctdna = 0
	);

baseline.smps <- unique(timeline[which(timeline$Group == 'baseline'),]$Sample);
ontrial.smps <- unique(timeline[which(timeline$Group == 'on.trial'),]$Sample);
eot.smps <- unique(timeline[which(timeline$Group == 'EOT'),]$Sample);

for (i in 1:nrow(overlap.data)) {

	patient <- which(validation$Patient == overlap.data[i,]$Patient);
	chrom <- which(validation$Chromosome == overlap.data[i,]$Chromosome);
	pos1 <- which(validation$Start_Position == overlap.data[i,]$Start_Position);
	pos2 <- which(validation$End_Position == overlap.data[i,]$End_Position);
	tp <- which(validation$ENSEMBLE.mod == 1);

	idx <- intersect(tp, intersect(patient, intersect(chrom, intersect(pos1, pos2))));

	# is it in any ctDNA sample?
	if (length(idx) > 0) {
		overlap.data[i,]$ctdna <- 1;
		} else {
		next;
		}

	# is it in the baseline sample?
	if (any(validation[idx,]$Sample %in% baseline.smps)) { overlap.data[i,]$baseline <- 1; }

	# is it in the cycle2 sample?
	if (any(validation[idx,]$Sample %in% ontrial.smps)) { overlap.data[i,]$cycle2 <- 1; }

	# is it in the final (EOT) sample?
	if (any(validation[idx,]$Sample %in% eot.smps)) { overlap.data[i,]$eot <- 1; }
	}

# indicate cases with no data
overlap.data[which(overlap.data$Patient %in% c('EVO-009-026','EVO-400-008')),]$baseline <- NA;
overlap.data[which(overlap.data$Patient %in% c('EVO-009-001','EVO-009-008','EVO-009-012','EVO-009-016','EVO-009-021','EVO-009-022','EVO-400-004','EVO-400-007')),]$cycle2 <- NA;
overlap.data[which(overlap.data$Patient %in% c('EVO-009-026','EVO-400-008')),]$eot <- NA;

write.table(
	overlap.data,
	file = generate.filename('EVOLVE_ctDNA','_validation_data__overlaps','tsv'),
	row.names = FALSE,
	col.names = TRUE,
	sep = '\t'
	);

library(VennDiagram);
venn.diagram(
	x = list(
		exome = which(overlap.data$exome == 1),
		eot = which(overlap.data$eot == 1),
		baseline = which(overlap.data$baseline == 1),
		"on-trial" = which(overlap.data$cycle2 == 1)
		),
	cex = 1.8,
	cat.cex = 1.8,
	filename = generate.filename('EVOLVE_ctDNA','_validation_venn','png'),
	resolution = 200,
	units = 'in',
	height = 6,
	width = 6
	);


