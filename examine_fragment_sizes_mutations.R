### examine_fragment_sizes_mutations.R #############################################################
# Examine differences in fragment sizes in cfDNA at known mutation sites.

### PREPARE SESSION ################################################################################
# import libraries
library(BoutrosLab.plotting.general);

source('/cluster/home/sprokope/git/analysis/helper_functions/session.functions.R');

working.dir <- '/cluster/projects/ovgroup/projects/OV_Superset/EVOLVE/ctDNA/pipeline_suite/fragment_size';

setwd(working.dir);

### READ DATA ######################################################################################
# look at fragment size at known mutation sites
maf <- read.delim('../Ensemble_calls/2022-07-26_EVOLVE_ctdna_plus_exome__mutation_data_cbioportal.txt');

# also get purity estimates
purity <- read.delim('../Ensemble_calls/tumour_content/2022-09-06_EVOLVE_ctDNA__estimated_tumour_content.tsv');

# and results from above analysis (global fragment sizes)
global.fs <- read.delim('2022-09-07_EVOLVE_ctDNA_fragment_size_summary.tsv');
global.fs$Patient <- substr(global.fs$Sample, 0, 11);

# finally, get fragment sizes for each variant
input.files <- list.files(pattern = 'alt.reads', recursive = TRUE);
input.files <- input.files[!grepl('Screening2', input.files)]; # junk sample (replicate bam for no reason)
input.files <- input.files[!grepl('EVO-400-004_ctDNA', input.files)]; # uses germline snps only

### FORMAT DATA ####################################################################################
# format mutation data
maf <- maf[-which(grepl('EVO-400-004', maf$Tumor_Sample_Barcode) & maf$Hugo_Symbol == 'BRCA2'),];
maf <- maf[which(maf$Mutation_Status == 'somatic'),];
maf <- maf[which(maf$Variant_Type == 'SNP'),];
maf$Patient <- substr(maf$Tumor_Sample_Barcode,0,11);

snp.data <- unique(maf[which(maf$Hugo_Symbol %in% c('TP53','BRCA1','BRCA2')),c('Patient','Hugo_Symbol','Chromosome','Start_Position','End_Position','Reference_Allele','Tumor_Seq_Allele2')]);
snp.data$MUTATION <- paste(
	snp.data$Chromosome,
	snp.data$Start_Position,
	snp.data$Reference_Allele,
	snp.data$Tumor_Seq_Allele2,
	sep = '_'
	);

known.sites <- snp.data[order(snp.data$Patient),c('Patient','MUTATION')];

# compare fragment sizes (alt vs ref)
stats.data <- merge(
	purity[,c('Patient.ID','Sample')],
	known.sites,
	by.x = 'Patient.ID',
	by.y = 'Patient',
	all.x = TRUE
	);

# indicate fields to fill in
stats.data[,c('Median.REF','Median.ALT','FC','p.value')] <- NA;

# read in and process data
mutation.data <- list();

for (file in input.files) {

	smp <- unlist(strsplit(basename(file), '_allUNIQUE'))[1];
	variant <- sub('_alt.reads','',unlist(strsplit(basename(file), '_allUNIQUE__'))[2]);

	idx <- which(stats.data$Sample == smp & stats.data$MUTATION == variant);

	# get fragment sizes for each mutation
	ref <- as.numeric(read.delim(sub('alt','ref',file), header = FALSE)$V1);
	alt <- as.numeric(read.delim(file, header = FALSE)$V1);

	# compare ref/alt at each site
	stats.data[idx,c('Median.REF','Median.ALT')] <- c(median(ref), median(alt));
	stats.data[idx,]$FC <- log2(median(alt) / median(ref));
	stats.data[idx,]$p.value <- wilcox.test(ref, alt)$p.value;

	# save data to per-sample ref/alt list
	if (smp %in% names(mutation.data)) {
		mutation.data[[smp]]$alt <- c(mutation.data[[smp]]$alt, alt);
		mutation.data[[smp]]$ref <- c(mutation.data[[smp]]$ref, ref);
		} else {
		mutation.data[[smp]]$alt <- alt;
		mutation.data[[smp]]$ref <- ref;
		}
	gc();
	}

write.table(
	stats.data,
	file = generate.filename('EVOLVE_ctDNA','_per_mutation_fragment_size_summary','tsv'),
	row.names = FALSE,
	col.names = TRUE,
	sep = '\t'
	);

# summarize data
results <- data.frame(
	Sample = names(mutation.data),
	Median.REF = NA,
	Median.ALT = NA,
	FC = NA,
	p.value = NA
	);

for (i in 1:nrow(results)) {
	smp <- results[i,]$Sample;
	results[i,]$Median.REF <- median(mutation.data[[smp]]$ref);
	results[i,]$Median.ALT <- median(mutation.data[[smp]]$alt);
	results[i,]$FC <- log2(results[i,]$Median.ALT / results[i,]$Median.REF);
	results[i,]$p.value <- wilcox.test(mutation.data[[smp]]$ref, mutation.data[[smp]]$alt)$p.value;
	}

combined <- merge(
	merge(global.fs, results, all = TRUE),
	purity[,c('Sample','Purity','Max.VAF')],
	all = TRUE
	);

combined <- combined[,c('Patient','Sample','Median','Prop.short','Prop.long','Group','Median.REF','Median.ALT','FC','p.value','Purity','Max.VAF')];

write.table(
	combined,
	file = generate.filename('EVOLVE_ctDNA','_fragment_size_and_purity_estimates','tsv'),
	row.names = FALSE,
	col.names = TRUE,
	sep = '\t'
	);

# this one case only had a single read; add more to allow plotting
mutation.data[['EVO-009-011_ctDNA_EOT']]$alt <- rep(mutation.data[['EVO-009-011_ctDNA_EOT']]$alt,2);

# do some quick tests
ref.data <- list();
alt.data <- list();

for (smp in names(mutation.data)) {
	ref.data[[smp]] <- mutation.data[[smp]]$ref;
	alt.data[[smp]] <- mutation.data[[smp]]$alt;
	}

baseline.idx <- grep('Screening|C1D1', names(ref.data));
on.trial.idx <- grep('C2D1', names(ref.data));
eot.idx <- setdiff(1:length(ref.data), c(baseline.idx, on.trial.idx));

wilcox.test(unlist(ref.data), unlist(alt.data))$p.value;
wilcox.test(unlist(ref.data[baseline.idx]), unlist(alt.data[baseline.idx]))$p.value;
wilcox.test(unlist(ref.data[on.trial.idx]), unlist(alt.data[on.trial.idx]))$p.value;
wilcox.test(unlist(ref.data[eot.idx]), unlist(alt.data[eot.idx]))$p.value;

# add group to results
results <- merge(results, unique(combined[,c('Sample','Group')]), all.x = TRUE);

### VISUALIZE DATA #################################################################################
# plot results
colour.scheme <- c('skyblue','seagreen','pink');
names(colour.scheme) <- c('baseline','on.trial','EOT');

line.key <- list(
	lines = list(
		size = 3,
		col = c('black','black','transparent', colour.scheme),
		lty = c(1,2,rep(1,4)),
		lwd = 1.5
		),
	text = list(
		lab = c('mutant','wildtype','','baseline','on-trial','end-of-treatment'),
		cex = 1
		)
	);

# plot each patient
patients <- unique(substr(names(mutation.data),0,11));

for (patient in patients) {

	if (all(is.na(combined[which(combined$Patient == patient),]$p.value))) { next; }
	my.colours <- c();
	max.x <- 350;

	plot.data <- list();
	if (any(combined[which(combined$Patient == patient),]$Group == 'baseline')) {
		smp <- combined[which(combined$Patient == patient & combined$Group == 'baseline'),]$Sample;
		plot.data[['baseline_alt']] <- mutation.data[[smp]]$alt;
		plot.data[['baseline_ref']] <- mutation.data[[smp]]$ref;
		my.colours <- c(my.colours, colour.scheme[1]);
		if (any(plot.data[['baseline_alt']] > 400) | any(plot.data[['baseline_ref']] > 400)) {
			max.x <- 500; }
		}
	if (any(combined[which(combined$Patient == patient),]$Group == 'on.trial')) {
		smp <- combined[which(combined$Patient == patient & combined$Group == 'on.trial'),]$Sample;
		plot.data[['cycle2_alt']] <- mutation.data[[smp]]$alt;
		plot.data[['cycle2_ref']] <- mutation.data[[smp]]$ref;
		my.colours <- c(my.colours, colour.scheme[2]);
		if (any(plot.data[['cycle2_alt']] > 400) | any(plot.data[['cycle2_ref']] > 400)) {
			max.x <- 500; }
		}
	if (any(combined[which(combined$Patient == patient),]$Group == 'eot')) {
		smp <- combined[which(combined$Patient == patient & combined$Group == 'eot'),]$Sample;
		plot.data[['eot_alt']] <- mutation.data[[smp]]$alt;
		plot.data[['eot_ref']] <- mutation.data[[smp]]$ref;
		my.colours <- c(my.colours, colour.scheme[3]);
		if (any(plot.data[['eot_alt']] > 400) | any(plot.data[['eot_ref']] > 400)) {
			max.x <- 500; }
		}

	tmp.plot <- create.densityplot(
		x = plot.data,
		col = rep(my.colours, each = 2),
		lty = rep(c(1,2), times = 3),
		lwd = 1.5,
		xlab.label = expression('Fragment Length'),
		ylab.label = expression('Density'),
		xaxis.cex = 1.2,
		yaxis.cex = 1.2,
		xlab.cex = 1.5,
		ylab.cex = 1.5,
		xaxis.fontface = 'plain',
		yaxis.fontface = 'plain',
		xaxis.tck = c(1,0),
		yaxis.tck = c(1,0),
		xlimits = c(0,max.x),
		xat = seq(0,500,100),
		abline.v = if (max.x == 350) { 167 } else { seq(0,550,167) },
		abline.lty = 2,
		abline.lwd = 0.5,
		abline.col = 'black',
		legend = list(
			inside = list(fun = draw.key, args = list(key = line.key),
				x = 0.6, y = 0.95),
			inside = list(fun = draw.key, args = list(key = list(text = list(lab = patient))),
				x = 0.02, y = 0.98)
			)
		);

	tmp.plot$par.settings$axis.line$lwd <- 1;

	write.plot(
		tmp.plot,
		filename = paste0(
			patient,'/',generate.filename('mutant_vs_wildtype','_fragment_lengths','png')),
		resolution = 200
		);	
	}

# format data for boxplots
plot.data <- reshape(
	combined[!is.na(combined$FC),c('Sample','Group','Median.REF','Median.ALT')],
	direction = 'long',
	varying = list(3:4),
	v.names = 'Median',
	timevar = 'base',
	times = c('wildtype','mutant')
	);

plot.data$base <- factor(plot.data$base, levels = c('wildtype','mutant'));

# create some boxplots
plot.objects <- list();

for (i in c('baseline','on.trial','EOT')) {
	
	label <- if (i == 'baseline') { i
		} else if (i == 'on.trial') { 'on-trial'
		} else { 'end-of-treatment' }

	stats.result <- wilcox.test(
		x = combined[which(combined$Group == i),]$Median.REF,
		y = combined[which(combined$Group == i),]$Median.ALT,
		paired = TRUE
		)$p.value;

	plot.objects[[i]] <- create.boxplot(
		Median ~ base,
		plot.data[which(plot.data$Group == i),],
		add.stripplot = TRUE,
		points.col = colour.scheme[i],
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
		text.labels = if (stats.result < 0.001) { '***' } else if (stats.result < 0.01) {
			'**' } else if (stats.result < 0.1) { '*' } else { 'ns' },
		text.x = 1.5,
		text.y = if (max(plot.data[which(plot.data$Group == i),]$Median) > 200) { 318 } else { 190 },
		text.cex = if (stats.result >= 0.1) { 0.8 } else { 1.2 },
		add.rectangle = TRUE,
		xleft.rectangle = 1,
		xright.rectangle = 2,
		ytop.rectangle = if (max(plot.data[which(plot.data$Group == i),]$Median) > 200) { 311 } else { 187 },
		ybottom.rectangle = if (max(plot.data[which(plot.data$Group == i),]$Median) > 200) { 310 } else { 186.7 },
		col.rectangle = 'black',
		#add.text = TRUE,
		#text.labels = display.statistical.result(
		#	stats.result,
		#	statistic.type = 'p',
		#	symbol = ' = '
		#	),
		#text.cex = 0.8,
		#text.x = 1.5,
		#text.y = if (max(plot.data[which(plot.data$Group == i),]$Median) > 200) { 300 } else { 192 },
		#text.fontface = 'italic',
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

# try a volcano plot?
create.scatterplot(
	-log10(p.value) ~ FC,
	combined,
	xlab.label = expression('log'['2']*' (ALT / REF)'),
	xlab.cex = 1.8,
	ylab.label = expression('-log'['10']*'(p-value)'),
	ylab.cex = 2,
	xlab.axis.padding = 4,
	ylab.axis.padding = 3,
	xaxis.tck = c(1,0),
	yaxis.tck = c(1,0),
	xaxis.cex = 1.5,
	yaxis.cex = 1.5,
	col = colour.scheme[match(combined$Group, names(colour.scheme))],
	cex = 1,
	alpha = 0.9,
	abline.h = -log10(0.05),
	abline.v = 0,
	abline.lty = 2,
	key = list(
		x = 0.6, y = 0.95,
		points = list(
			col = colour.scheme,
			cex = 1, pch = 19
			),
		text = list(
			lab = c('baseline','on-trial','end-of-treatment'),
			cex = 1.2
			)
		),
	height = 6,
	width = 6,
	style = 'Nature',
	filename = generate.filename('EVOLVE_ctDNA','_mutation_fragment_sizes_volcano', 'png')
	);

combined <- combined[!is.na(combined$p.value),];

for (i in c('baseline','on.trial','EOT')) {
	
	xlimits <- if (i == 'baseline') { c(-0.2, 0.15) 
		} else if (i == 'on.trial') { c(-0.3, 0.9)
		} else { c(-0.3, 0.1) }
	xat <- if (i == 'baseline') { seq(-0.2, 0.15, 0.1) 
		} else if (i == 'on.trial') { seq(-0.3, 0.9, 0.3)
		} else { seq(-0.3, 0.1, 0.1) }

	tmp.plot <- create.scatterplot(
		-log10(p.value) ~ FC,
		combined[which(combined$Group == i),],
		xlab.label = expression('log'['2']*' (alt / ref)'),
		xlab.cex = 1.8,
		ylab.label = expression('-log'['10']*'(p-value)'),
		ylab.cex = 2,
		ylab.axis.padding = 3,
		xaxis.tck = c(1,0),
		yaxis.tck = c(1,0),
		xaxis.cex = 1.5,
		yaxis.cex = 1.5,
		xlimits = xlimits,
		xat = xat,
		col = sapply(combined[which(combined$Group == i),]$p.value, function(i) {
			if (i >= 0.05) { 'turquoise3'} else { 'violetred3' } } ),
		cex = 1.2,
		abline.h = -log10(0.05),
		abline.v = 0,
		abline.lty = 2,
		style = 'Nature'
		);

	tmp.plot$par.settings$layout.heights$axis.xlab.padding <- 3;

	write.plot(
		tmp.plot,
		filename = generate.filename('EVOLVE_ctDNA',
			paste0('_mutation_fragment_sizes_volcano',i), 'png')
		);
	}

### SAVE SESSION INFO ##############################################################################
save.session.profile(generate.filename('FSbyMutationSummary','SessionProfile','txt'));
