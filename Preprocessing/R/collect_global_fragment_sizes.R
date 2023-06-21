### collect_global_fragment_sizes.R ################################################################
# Collect global fragment sizes from tumour and healthy control samples

### PREPARE SESSION ################################################################################
# import libraries
library(BoutrosLab.plotting.general);

source('/cluster/home/sprokope/git/analysis/helper_functions/session.functions.R');

working.dir <- '/cluster/projects/ovgroup/projects/OV_Superset/EVOLVE/ctDNA/pipeline_suite/fragment_size';

setwd(working.dir);

### READ DATA ######################################################################################
# read in timeline information
timeline <- read.delim('../configs/2022-09-06_sample_info_with_batch.txt', stringsAsFactors = FALSE);

# find data files
fragment.files <- list.files(pattern = 'insertsize.tsv', recursive = TRUE);

fragment.files <- fragment.files[!grepl('Screening2',fragment.files)];

tumour.data <- list();
for (i in fragment.files) {
	smp <- unlist(strsplit(basename(i), '_allUNIQUE'))[1];
	tumour.data[[smp]] <- as.numeric(read.delim(i, header = FALSE)$V1);
	gc();
	}

# and find normal files
normal.files <- list.files(path = 'PanelOfNormals', pattern = 'txt', full.names = TRUE);

normal.data <- list();
for (i in normal.files) {
	smp <- substr(basename(i),0,10);
	normal.data[[smp]] <- as.numeric(read.delim(i, header = FALSE)$V1);
	}

### FORMAT DATA ####################################################################################
# summarize data
results <- rbind(
	data.frame(
		Sample = names(tumour.data),
		Median = unlist(lapply(tumour.data, median)),
		Prop.short = unlist(lapply(tumour.data, function(i) {
			length(i[which(i < 150)]) / length(i) } )),
		Prop.norm = unlist(lapply(tumour.data, function(i) {
			length(i[which(i >= 150 & i <= 250)]) / length(i) } )),
		Prop.long = unlist(lapply(tumour.data, function(i) { 
			length(i[which(i < 320 & i > 250)]) / length(i) } )),
		Prop.plus = unlist(lapply(tumour.data, function(i) {
			length(i[which(i >= 320)]) / length(i) } )),
		Ratio.short = unlist(lapply(tumour.data, function(i) {
			length(i[which(i > 90 & i <= 150)]) / 
			length(i[which(i > 150 & i <= 230)]) } )),
		Ratio.long = unlist(lapply(tumour.data, function(i) {
			length(i[which(i > 230 & i <= 320)]) / 
			length(i[which(i > 150 & i <= 230)]) } ))
		),
	data.frame(
		Sample = names(normal.data),
		Median = unlist(lapply(normal.data, median)),
		Prop.short = unlist(lapply(normal.data, function(i) {
			length(i[which(i < 150)]) / length(i) } )),
		Prop.norm = unlist(lapply(normal.data, function(i) {
			length(i[which(i >= 150 & i <= 250)]) / length(i) } )),
		Prop.long = unlist(lapply(normal.data, function(i) { 
			length(i[which(i < 320 & i > 250)]) / length(i) } )),
		Prop.plus = unlist(lapply(normal.data, function(i) {
			length(i[which(i >= 320)]) / length(i) } )),
		Ratio.short = unlist(lapply(normal.data, function(i) {
			length(i[which(i > 90 & i <= 150)]) / 
			length(i[which(i > 150 & i <= 230)]) } )),
		Ratio.long = unlist(lapply(normal.data, function(i) {
			length(i[which(i > 230 & i <= 320)]) / 
			length(i[which(i > 150 & i <= 230)]) } ))
		)
	);

## get densities per sample for combined plots
# summarize normals
normal.summary <- data.frame(density(normal.data[[1]], from = 20, to = 350, n = 512)[c('x','y')]);
colnames(normal.summary) <- c('Breakpoints', names(normal.data)[1]);

for (i in 2:length(normal.data)) {
	tmp <- data.frame(density(normal.data[[i]], from = 20, to = 350, n = 512)[c('x','y')]);
	colnames(tmp)[2] <- names(normal.data)[i];

	normal.summary <- merge(normal.summary, tmp, by.x = 'Breakpoints', by.y = 'x');
	}

normal.summary$Median <- apply(normal.summary[,names(normal.data)],1,median);

# summarize tumours
tumour.summary <- data.frame(density(tumour.data[[1]], from = 20, to = 350, n = 512)[c('x','y')]);
colnames(tumour.summary) <- c('Breakpoints', names(tumour.data)[1]);

for (i in 2:length(tumour.data)) {
	tmp <- data.frame(density(tumour.data[[i]], from = 20, to = 350, n = 512)[c('x','y')]);
	colnames(tmp)[2] <- names(tumour.data)[i];

	tumour.summary <- merge(tumour.summary, tmp, by.x = 'Breakpoints', by.y = 'x');
	}

# format data to save
colnames(tumour.summary) <- gsub('_ctDNA', '', colnames(tumour.summary));
colnames(normal.summary) <- gusb('TGL49','HC', colnames(normal.summary));

fragment.summary <- results;
fragment.summary$Sample <- gsub('_ctDNA','', fragment.summary$Sample);
fragment.summary$Sample <- gsub('TGL49','HC', fragment.summary$Sample);

save(
	tumour.summary,
	normal.summary,
	fragment.summary,
	file = generate.filename('EVOLVE_ctDNA','fragmentation_summary','RData')
	);

### SAVE SESSION INFO ##############################################################################
save.session.profile(generate.filename('FragmentSizeSummary','SessionProfile','txt'));
