### plot_interval_coverage.R #######################################################################
# Create plot for mean interval coverage to assess CCNE1 amplification

### 
library(BoutrosLab.plotting.general);

setwd('');

### READ DATA
# get intervals
intervals <- read.delim('CHARM-MMR_plus_EVOLVE_hg38_filtered.bed', header = F);

# get coverage values
cov.data <- do.call(rbind, lapply(list.files(pattern = 'interval_meanCoverage.tsv', recursive = TRUE), function(i) { x <- read.delim(i, header = F)[,c(1:4,7)]; colnames(x) <- c('Chromosome','Start','End','ID','Coverage'); x$Sample <- gsub('_interval_meanCoverage.tsv','',basename(i)); return(x); } ));

cov.data$ID <- factor(cov.data$ID, levels = as.character(intervals$V4));

# is CCNE1 coverage different?
stats.data <- data.frame(Sample = unique(cov.data$Sample), CCNE1 = NA, Other = NA, p = NA);

for (i in 1:nrow(stats.data)) {
	smp <- stats.data[i,]$Sample
	tmp.data <- cov.data[which(cov.data$Sample == smp),]
	results <- t.test(
		tmp.data[grepl('CCNE1', tmp.data$ID),]$Coverage,
		tmp.data[!grepl('CCNE1', tmp.data$ID),]$Coverage
		);
	stats.data[i,2:4] <- c(results$estimate[[1]], results$estimate[[2]],results$p.value)
	}

# plot data
plot.data <- cov.data[grepl('allUNIQUE', cov.data$Sample),];
plot.data$Patient <- substr(plot.data$Sample, 0, 11);

create.scatterplot(
	Coverage ~ ID | Patient,
	plot.data,
	groups = plot.data$Sample,
	type = 'l',
	file = 'mean_coverage_per_interval.png',
	resolution = 200,
	width = 15,
	height = 10,
	layout = c(2,15),
	as.table = TRUE,
	x.spacing = 0.5,
	y.spacing = 0.1,
	xaxis.lab = rep('',nrow(intervals)),
	xaxis.tck = 0,
	xlab.label = 'target intervals',
	xlab.cex = 1,
	yaxis.tck = 0,
	ylab.label = 'Mean coverage',
	ylab.cex = 1,
	yaxis.lab = NULL,
	cex = 0.5,
	col = c('black','grey30','grey60','grey90'),
	add.rectangle = TRUE,
	xleft.rectangle = 1265,
	xright.rectangle = 1288,
	ybottom.rectangle = 0,
	ytop.rectangle = 10**4,
	col.rectangle = 'red',
	alpha.rectangle = 0.1
	);

timeline <- read.delim('../cBioportal/EVOLVE_ctDNA_20220622/data_clinical_timeline_specimen.txt');
timeline <- timeline[order(timeline$PATIENT_ID, timeline$START_DATE),];

plot.data <- cov.data[grepl('allUNIQUE', cov.data$Sample),];
plot.data$Patient <- substr(plot.data$Sample, 0, 11);
plot.data$Sample <- gsub('_allUNIQUE','',plot.data$Sample);

timeline <- timeline[which(timeline$SAMPLE_ID %in% plot.data$Sample),];

plot.data <- merge(
	plot.data,
	timeline[,c('PATIENT_ID','SAMPLE_ID','START_DATE')],
	by.x = c('Patient','Sample'),
	by.y = c('PATIENT_ID','SAMPLE_ID')
	);

plot.data <- droplevels(plot.data);
plot.data$Sample <- factor(plot.data$Sample, levels = as.character(timeline$SAMPLE_ID));
plot.data <- plot.data[order(plot.data$Sample),];

plot.data$Gene <- 0;
plot.data[grepl('CCNE1',plot.data$ID),]$Gene <- 1;

plot.data$Order <- NA;
for (patient in unique(plot.data$Patient)) {
	idx <- which(plot.data$Patient == patient);
	tmp <- droplevels(plot.data[idx,]);

	plot.data[idx,]$Order <- paste0(as.numeric(tmp$Sample),'.',tmp$Gene);
	}

create.boxplot(
	Coverage ~ Order | Patient,
	plot.data,
	x.relation = 'free',
	col = c('transparent','pink'),
	file = 'mean_coverage_per_interval_boxplot.png',
	resolution = 200,
	width = 10,
	height = 15,
	layout = c(5,6),
	as.table = TRUE,
	x.spacing = 0.5,
	y.spacing = 0.1,
	xaxis.lab = rep('',nrow(intervals)),
	xaxis.tck = 0,
	xlab.label = 'Timepoint',
	xlab.cex = 1,
	yaxis.tck = 0,
	ylab.label = 'Mean coverage',
	ylab.cex = 1,
	yaxis.lab = NULL
	);

