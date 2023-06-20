### examine_fragment_sizes_mutations.R #############################################################
# Examine differences in fragment sizes in cfDNA at somatic mutation sites.

### PREPARE SESSION ################################################################################
# import libraries
library(BoutrosLab.plotting.general);

source('/cluster/home/sprokope/git/analysis/helper_functions/session.functions.R');

working.dir <- '/cluster/projects/ovgroup/projects/OV_Superset/EVOLVE/ctDNA/pipeline_suite/fragment_size';

setwd(working.dir);

### READ DATA ######################################################################################
# look at fragment size at known mutation sites
maf <- read.delim('../cBioportal/EVOLVE_ctDNA_20220908/mutation_data_extended.txt');

# also get purity estimates
purity <- read.delim('../Ensemble_calls/tumour_content/2022-09-06_EVOLVE_ctDNA__estimated_tumour_content.tsv');

# and results from above analysis (global fragment sizes)
global.fs <- read.delim('2022-09-07_EVOLVE_ctDNA_fragment_size_summary.tsv');
global.fs$Patient <- substr(global.fs$Sample, 0, 11);

# finally, get fragment sizes for each variant
input.files <- list.files(pattern = 'alt.reads', recursive = TRUE);
input.files <- input.files[!grepl('Screening2', input.files)]; # junk sample (replicate bam for no reason)
input.files <- input.files[!grepl('reversion', input.files)]; # ignore these ones

### FORMAT DATA ####################################################################################
# format mutation data
maf <- maf[which(maf$Mutation_Status == 'somatic'),];
maf <- maf[which(maf$Variant_Type == 'SNP'),];

maf$Patient <- substr(maf$Tumor_Sample_Barcode,0,11);

exome.data <- unique(maf[!grepl('ctDNA',maf$Tumor_Sample_Barcode),c('Patient','Hugo_Symbol','Chromosome','Start_Position','End_Position','Reference_Allele','Tumor_Seq_Allele2')]);
novel.data <- unique(maf[grepl('ctDNA',maf$Tumor_Sample_Barcode),c('Patient','Hugo_Symbol','Chromosome','Start_Position','End_Position','Reference_Allele','Tumor_Seq_Allele2')]);
exome.data$Exome <- 1;
novel.data$Exome <- 0;

snp.data <- rbind(exome.data, novel.data);
snp.data$MUTATION <- paste(
	snp.data$Chromosome,
	snp.data$Start_Position,
	snp.data$Reference_Allele,
	snp.data$Tumor_Seq_Allele2,
	sep = '_'
	);

snp.data <- aggregate(Exome ~ Patient + Hugo_Symbol + MUTATION, snp.data, max);

# compare fragment sizes (alt vs ref)
stats.data <- merge(
	purity[,c('Patient.ID','Sample','Group.ctdna')],
	snp.data[,c('Patient','Hugo_Symbol','MUTATION')],
	by.x = 'Patient.ID',
	by.y = 'Patient',
	all.x = TRUE
	);

# indicate fields to fill in
stats.data[,c('Median.REF','Median.ALT','Median.FC','p.value','short.ratio','Mean.REF','Mean.ALT')] <- NA;

# read in and process data
fragment.data <- list();

for (file in input.files) {

	smp <- unlist(strsplit(basename(file), '_allUNIQUE'))[1];
	variant <- sub('_alt.reads','',unlist(strsplit(basename(file), '_allUNIQUE__'))[2]);

	idx <- which(stats.data$Sample == smp & stats.data$MUTATION == variant);
	if (length(idx) == 0) { next; }

	# get fragment sizes for each mutation
	ref <- as.numeric(read.delim(sub('alt','ref',file), header = FALSE)$V1);
	alt <- tryCatch(
		expr = as.numeric(read.delim(file, header = FALSE)$V1),
		error = function(e) { NA }
		);

	# compare ref/alt at each site
	stats.data[idx,c('Median.REF','Median.ALT')] <- c(median(ref), median(alt));
	stats.data[idx,c('Mean.REF','Mean.ALT')] <- c(mean(ref), mean(alt));
	stats.data[idx,]$Median.FC <- median(alt) / median(ref);
	stats.data[idx,]$p.value <- tryCatch(
		expr = wilcox.test(ref, alt)$p.value,
		error = function(e) { NA }
		);

	# find ratio of short:normal reads (overall at site)
	all.reads <- c(ref,alt);
	stats.data[idx,]$short.ratio <- length(all.reads[which(all.reads >= 90 & all.reads < 151)]) /
		length(all.reads[which(all.reads >= 151 & all.reads < 231)]);

	# save data to per-sample ref/alt list
	if (smp %in% names(fragment.data)) {
		fragment.data[[smp]]$alt <- c(fragment.data[[smp]]$alt, alt);
		fragment.data[[smp]]$ref <- c(fragment.data[[smp]]$ref, ref);
		} else {
		fragment.data[[smp]]$alt <- alt;
		fragment.data[[smp]]$ref <- ref;
		}

	gc();
	}

na.counts <- apply(stats.data[,c('Median.REF','Median.ALT')], 1, function(i) { all(is.na(i)) } );
stats.data <- stats.data[!na.counts,];

write.table(
	stats.data,
	file = generate.filename('EVOLVE_ctDNA','_known_mutation_fragment_size_summary','tsv'),
	row.names = FALSE,
	col.names = TRUE,
	sep = '\t'
	);

# summarize data
results <- data.frame(
	Sample = names(mutation.data),
	Median.REF = NA,
	Median.ALT = NA,
	Median.FC = NA,
	p.value = NA,
	Mean.REF = NA,
	Mean.ALT = NA,
	short.ratio = NA
	);

for (i in 1:nrow(results)) {
	smp <- results[i,]$Sample;
	results[i,]$Median.REF <- median(mutation.data[[smp]]$ref);
	results[i,]$Median.ALT <- median(mutation.data[[smp]]$alt);
	results[i,]$Median.FC <- results[i,]$Median.ALT / results[i,]$Median.REF;
	results[i,]$p.value <- tryCatch(
		expr = wilcox.test(mutation.data[[smp]]$ref, mutation.data[[smp]]$alt)$p.value,
		error = function(e) { NA }
		);
	results[i,]$Mean.REF <- mean(mutation.data[[smp]]$ref);
	results[i,]$Mean.ALT <- mean(mutation.data[[smp]]$alt);
	all.reads <- c(mutation.data[[smp]]$ref,mutation.data[[smp]]$alt);
	results[i,]$short.ratio <- length(all.reads[which(all.reads >= 90 & all.reads < 151)]) /
		length(all.reads[which(all.reads >= 151 & all.reads < 231)]);
	}

combined <- merge(
		xaxis.rot = 90,
		xlab.label = NULL,
		xlab.top.label = label,
		xlab.top.cex = 1.2,
		xlab.top.y = 1,
		ylab.label = if (i == 'baseline') { 'Median Fragment Size' } else { NULL },
		ylab.cex = 1.5,
		ylab.axis.padding = 3,
		xaxis.tck = c(1,0),
		yaxis.tck = c(1,0),
		xaxis.cex = 1.2,
		yaxis.cex = 1.2,
		ylimits = if (max(plot.data[which(plot.data$Group == i),]$Median) > 200) {
			c(130,330) } else { c(135,195) },
		add.text = TRUE,
#		text.labels = if (stats.result < 0.001) { '***' } else if (stats.result < 0.01) {
#			'**' } else if (stats.result < 0.1) { '*' } else { 'ns' },
		text.labels = display.statistical.result(
			stats.result,
			statistic.type = 'p',
			symbol = ' = '
			),
		text.x = 1.5,
		text.y = if (max(plot.data[which(plot.data$Group == i),]$Median) > 200) { 318 } else { 190 },
		text.cex = 0.8, # if (stats.result >= 0.1) { 0.8 } else { 1.2 },
#		add.rectangle = TRUE,
#		xleft.rectangle = 1,
#		xright.rectangle = 2,
#		ytop.rectangle = if (max(plot.data[which(plot.data$Group == i),]$Median) > 200) { 311 } else { 187 },
#		ybottom.rectangle = if (max(plot.data[which(plot.data$Group == i),]$Median) > 200) { 310 } else { 186.7 },
#		col.rectangle = 'black',
		style = 'Nature'	
		);
	}

create.multipanelplot(
	plot.objects = plot.objects,
	layout.width = 3,
	layout.height = 1,
	plot.objects.widths = c(1.1,1,1),
	left.legend.padding = 0,
	right.legend.padding = 0,
	bottom.legend.padding = 0,
	top.legend.padding = 0,
	ylab.axis.padding = 3,
	height = 4,
	width = 8,
	style = 'Nature',
	filename = generate.filename('EVOLVE_ctDNA','_mutation_median_fragment_sizes', 'png')
	);

### SAVE SESSION INFO ##############################################################################
save.session.profile(generate.filename('FSbyMutationSummary','SessionProfile','txt'));
### examine_fragment_sizes_mutations.R #############################################################
# Examine differences in fragment sizes in cfDNA at somatic mutation sites.

### PREPARE SESSION ################################################################################
# import libraries
library(BoutrosLab.plotting.general);

source('/cluster/home/sprokope/git/analysis/helper_functions/session.functions.R');

working.dir <- '/cluster/projects/ovgroup/projects/OV_Superset/EVOLVE/ctDNA/pipeline_suite/fragment_size';

setwd(working.dir);

### READ DATA ######################################################################################
# look at fragment size at known mutation sites
maf <- read.delim('../cBioportal/EVOLVE_ctDNA_20220908/mutation_data_extended.txt');

# also get purity estimates
purity <- read.delim('../Ensemble_calls/tumour_content/2022-09-06_EVOLVE_ctDNA__estimated_tumour_content.tsv');

# and results from above analysis (global fragment sizes)
global.fs <- read.delim('2022-09-07_EVOLVE_ctDNA_fragment_size_summary.tsv');
global.fs$Patient <- substr(global.fs$Sample, 0, 11);

# finally, get fragment sizes for each variant
input.files <- list.files(pattern = 'alt.reads', recursive = TRUE);
input.files <- input.files[!grepl('Screening2', input.files)]; # junk sample (replicate bam for no reason)
input.files <- input.files[!grepl('reversion', input.files)]; # ignore these ones

### FORMAT DATA ####################################################################################
# format mutation data
maf <- maf[which(maf$Mutation_Status == 'somatic'),];
maf <- maf[which(maf$Variant_Type == 'SNP'),];

maf$Patient <- substr(maf$Tumor_Sample_Barcode,0,11);

exome.data <- unique(maf[!grepl('ctDNA',maf$Tumor_Sample_Barcode),c('Patient','Hugo_Symbol','Chromosome','Start_Position','End_Position','Reference_Allele','Tumor_Seq_Allele2')]);
novel.data <- unique(maf[grepl('ctDNA',maf$Tumor_Sample_Barcode),c('Patient','Hugo_Symbol','Chromosome','Start_Position','End_Position','Reference_Allele','Tumor_Seq_Allele2')]);
exome.data$Exome <- 1;
novel.data$Exome <- 0;

snp.data <- rbind(exome.data, novel.data);
snp.data$MUTATION <- paste(
	snp.data$Chromosome,
	snp.data$Start_Position,
	snp.data$Reference_Allele,
	snp.data$Tumor_Seq_Allele2,
	sep = '_'
	);

snp.data <- aggregate(Exome ~ Patient + Hugo_Symbol + MUTATION, snp.data, max);

# compare fragment sizes (alt vs ref)
stats.data <- merge(
	purity[,c('Patient.ID','Sample','Group.ctdna')],
	snp.data[,c('Patient','Hugo_Symbol','MUTATION')],
	by.x = 'Patient.ID',
	by.y = 'Patient',
	all.x = TRUE
	);

# indicate fields to fill in
stats.data[,c('Median.REF','Median.ALT','Median.FC','p.value','short.ratio','Mean.REF','Mean.ALT')] <- NA;

# read in and process data
fragment.data <- list();

for (file in input.files) {

	smp <- unlist(strsplit(basename(file), '_allUNIQUE'))[1];
	variant <- sub('_alt.reads','',unlist(strsplit(basename(file), '_allUNIQUE__'))[2]);

	idx <- which(stats.data$Sample == smp & stats.data$MUTATION == variant);
	if (length(idx) == 0) { next; }

	# get fragment sizes for each mutation
	ref <- as.numeric(read.delim(sub('alt','ref',file), header = FALSE)$V1);
	alt <- tryCatch(
		expr = as.numeric(read.delim(file, header = FALSE)$V1),
		error = function(e) { NA }
		);

	# compare ref/alt at each site
	stats.data[idx,c('Median.REF','Median.ALT')] <- c(median(ref), median(alt));
	stats.data[idx,c('Mean.REF','Mean.ALT')] <- c(mean(ref), mean(alt));
	stats.data[idx,]$Median.FC <- median(alt) / median(ref);
	stats.data[idx,]$p.value <- tryCatch(
		expr = wilcox.test(ref, alt)$p.value,
		error = function(e) { NA }
		);

	# find ratio of short:normal reads (overall at site)
	all.reads <- c(ref,alt);
	stats.data[idx,]$short.ratio <- length(all.reads[which(all.reads >= 90 & all.reads < 151)]) /
		length(all.reads[which(all.reads >= 151 & all.reads < 231)]);

	# save data to per-sample ref/alt list
	if (smp %in% names(fragment.data)) {
		fragment.data[[smp]]$alt <- c(fragment.data[[smp]]$alt, alt);
		fragment.data[[smp]]$ref <- c(fragment.data[[smp]]$ref, ref);
		} else {
		fragment.data[[smp]]$alt <- alt;
		fragment.data[[smp]]$ref <- ref;
		}

	gc();
	}

na.counts <- apply(stats.data[,c('Median.REF','Median.ALT')], 1, function(i) { all(is.na(i)) } );
stats.data <- stats.data[!na.counts,];

write.table(
	stats.data,
	file = generate.filename('EVOLVE_ctDNA','_known_mutation_fragment_size_summary','tsv'),
	row.names = FALSE,
	col.names = TRUE,
	sep = '\t'
	);

# summarize data
results <- data.frame(
	Sample = names(mutation.data),
	Median.REF = NA,
	Median.ALT = NA,
	Median.FC = NA,
	p.value = NA,
	Mean.REF = NA,
	Mean.ALT = NA,
	short.ratio = NA
	);

for (i in 1:nrow(results)) {
	smp <- results[i,]$Sample;
	results[i,]$Median.REF <- median(mutation.data[[smp]]$ref);
	results[i,]$Median.ALT <- median(mutation.data[[smp]]$alt);
	results[i,]$Median.FC <- results[i,]$Median.ALT / results[i,]$Median.REF;
	results[i,]$p.value <- tryCatch(
		expr = wilcox.test(mutation.data[[smp]]$ref, mutation.data[[smp]]$alt)$p.value,
		error = function(e) { NA }
		);
	results[i,]$Mean.REF <- mean(mutation.data[[smp]]$ref);
	results[i,]$Mean.ALT <- mean(mutation.data[[smp]]$alt);
	all.reads <- c(mutation.data[[smp]]$ref,mutation.data[[smp]]$alt);
	results[i,]$short.ratio <- length(all.reads[which(all.reads >= 90 & all.reads < 151)]) /
		length(all.reads[which(all.reads >= 151 & all.reads < 231)]);
	}

combined <- merge(
	merge(global.fs, results, all = TRUE),
	purity[,c('Sample','Purity','Max.VAF')],
	all = TRUE
	);

combined <- combined[,c('Patient','Sample','Median','Prop.short','Prop.long','Group','Median.REF','Median.ALT','Median.FC','p.value','short.ratio','Purity','Max.VAF')];

write.table(
	combined,
	file = generate.filename('EVOLVE_ctDNA','_fragment_size_and_purity_estimates','tsv'),
	row.names = FALSE,
	col.names = TRUE,
	sep = '\t'
	);

	filename = generate.filename('EVOLVE_ctDNA','_mutation_median_fragment_sizes', 'png')
	);

### SAVE SESSION INFO ##############################################################################
save.session.profile(generate.filename('FSbyMutationSummary','SessionProfile','txt'));
