### examine_fragment_sizes.R #######################################################################
# Attempt to use fragment sizes to estimate tumour fraction in cfDNA.

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
#fragment.files <- fragment.files[grepl('Screening|C1D1',fragment.files)];

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

save(
	tumour.data,
	normal.data,
	file = generate.filename('EVOLVE_ctDNA','fragment_sizes','RData')
	);

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

results <- merge(
	results,
	timeline[,c('Sample','Group','Batch')],
	by = 'Sample',
	all.x = TRUE
	);
results[grep('TGL', results$Sample),]$Group <- 'normal';
results$Group <- factor(results$Group, levels = c('normal', 'baseline','on.trial','EOT'));

write.table(
	results,
	file = generate.filename('EVOLVE_ctDNA','fragment_size_summary','tsv'),
	row.names = FALSE,
	col.names = TRUE,
	sep = '\t'
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

normal.summary$Normal <- apply(normal.summary[,names(normal.data)],1,median);

# summarize tumours
tumour.summary <- data.frame(density(tumour.data[[1]], from = 20, to = 350, n = 512)[c('x','y')]);
colnames(tumour.summary) <- c('Breakpoints', names(tumour.data)[1]);

for (i in 2:length(tumour.data)) {
	tmp <- data.frame(density(tumour.data[[i]], from = 20, to = 350, n = 512)[c('x','y')]);
	colnames(tmp)[2] <- names(tumour.data)[i];

	tumour.summary <- merge(tumour.summary, tmp, by.x = 'Breakpoints', by.y = 'x');
	}

# get some comparison metrics
normal.idx <- which(results$Group == 'normal');
baseline.idx <- which(results$Group == 'baseline');
on.trial.idx <- which(results$Group == 'on.trial');
eot.idx <- which(results$Group == 'EOT');

wilcox.test(results[baseline.idx,]$Prop.short, results[normal.idx,]$Prop.short)$p.value;
wilcox.test(results[on.trial.idx,]$Prop.short, results[normal.idx,]$Prop.short)$p.value;
wilcox.test(results[eot.idx,]$Prop.short, results[normal.idx,]$Prop.short)$p.value;

wilcox.test(results[baseline.idx,]$Prop.long, results[normal.idx,]$Prop.long)$p.value;
wilcox.test(results[on.trial.idx,]$Prop.long, results[normal.idx,]$Prop.long)$p.value;
wilcox.test(results[eot.idx,]$Prop.long, results[normal.idx,]$Prop.long)$p.value;

wilcox.test(results[baseline.idx,]$Ratio.short, results[normal.idx,]$Ratio.short)$p.value;
wilcox.test(results[on.trial.idx,]$Ratio.short, results[normal.idx,]$Ratio.short)$p.value;
wilcox.test(results[eot.idx,]$Ratio.short, results[normal.idx,]$Ratio.short)$p.value;

wilcox.test(results[baseline.idx,]$Ratio.long, results[normal.idx,]$Ratio.long)$p.value;
wilcox.test(results[on.trial.idx,]$Ratio.long, results[normal.idx,]$Ratio.long)$p.value;
wilcox.test(results[eot.idx,]$Ratio.long, results[normal.idx,]$Ratio.long)$p.value;

### VISUALIZE DATA #################################################################################
# start with boxplots to compare proportion of short/long between tumour/normal
create.boxplot(
	Prop.short ~ Group,
	results,
	add.stripplot = TRUE,
	points.cex = 1,
	points.col = c('grey70','skyblue','seagreen','pink')[match(results$Group, c('normal','baseline','on.trial','EOT'))],
	xaxis.rot = 45,
	xlab.label = NULL,
	ylab.label = expression('Proportion fragments <150bp'),
	ylab.cex = 1.5,
	ylab.axis.padding = 3,
	xaxis.cex = 1.2,
	yaxis.cex = 1.2,
	xaxis.lab = gsub('\\.','-',levels(results$Group)),
	ylimits = c(0.08,0.36),
	yat = seq(0.1,0.35,0.05),
	add.text = TRUE,
	text.labels = c('**','ns','***'),
	text.x = c(1.5, 2.5, 3.5),
	text.y = c(0.355,0.346,0.335),
	text.cex = c(1.2,0.8,1.2),
	add.rectangle = TRUE,
	xleft.rectangle = c(1, 1, 1),
	xright.rectangle = c(2, 3, 4),
	ytop.rectangle = c(0.351, 0.341, 0.331),
	ybottom.rectangle = c(0.350, 0.340, 0.330),
	col.rectangle = 'black',
	style = 'Nature',
	filename = generate.filename('EVOLVE_ctDNA','_proportion150_boxplot', 'png')
	);

create.boxplot(
	Prop.long ~ Group,
	results,
	add.stripplot = TRUE,
	points.cex = 1,
	points.col = c('grey70','skyblue','seagreen','pink')[match(results$Group, c('normal','baseline','on.trial','EOT'))],
	xaxis.rot = 45,
	xlab.label = NULL,
	ylab.label = expression('Proportion fragments 250-320bp'),
	ylab.cex = 1.5,
	ylab.axis.padding = 3,
	xaxis.cex = 1.2,
	yaxis.cex = 1.2,
	xaxis.lab = gsub('\\.','-',levels(results$Group)),
	ylimits = c(-0.005,0.14),
	yat = seq(0,0.12,0.02),
	add.text = TRUE,
	text.labels = c('**','***','*'),
	text.x = c(1.5, 2.5, 3.5),
	text.y = c(0.135,0.130,0.125),
	text.cex = 1.2,
	add.rectangle = TRUE,
	xleft.rectangle = c(1, 1, 1),
	xright.rectangle = c(2, 3, 4),
	ytop.rectangle = c(0.133, 0.128, 0.123),
	ybottom.rectangle = c(0.1325, 0.1275, 0.1225),
	col.rectangle = 'black',
	style = 'Nature',
	filename = generate.filename('EVOLVE_ctDNA','_proportion_250-320_boxplot', 'png')
	);

create.boxplot(
	Ratio.short ~ Group,
	results,
	add.stripplot = TRUE,
	points.cex = 1,
	points.col = c('grey70','skyblue','seagreen','pink')[match(results$Group, c('normal','baseline','on.trial','EOT'))],
	xaxis.rot = 45,
	xlab.label = NULL,
	ylab.label = expression('Ratio (short:normal)'),
	ylab.cex = 1.5,
	ylab.axis.padding = 3,
	xaxis.cex = 1.2,
	yaxis.cex = 1.2,
	xaxis.lab = gsub('\\.','-',levels(results$Group)),
	ylimits = c(-0.005,0.6),
	yat = seq(0,0.6,0.1),
	add.text = TRUE,
	text.labels = c('***','ns','***'),
	text.x = c(1.5, 2.5, 3.5),
	text.y = c(0.58,0.565,0.545),
	text.cex = c(1.2,0.8,1.2),
	add.rectangle = TRUE,
	xleft.rectangle = c(1, 1, 1),
	xright.rectangle = c(2, 3, 4),
	ytop.rectangle = c(0.573, 0.555, 0.537),
	ybottom.rectangle = c(0.572, 0.554, 0.536),
	col.rectangle = 'black',
	style = 'Nature',
	filename = generate.filename('EVOLVE_ctDNA','_ratio_short_boxplot', 'png')
	);

create.boxplot(
	Ratio.long ~ Group,
	results,
	add.stripplot = TRUE,
	points.cex = 1,
	points.col = c('grey70','skyblue','seagreen','pink')[match(results$Group, c('normal','baseline','on.trial','EOT'))],
	xaxis.rot = 45,
	xlab.label = NULL,
	ylab.label = expression('Ratio (long:normal)'),
	ylab.cex = 1.5,
	ylab.axis.padding = 3,
	xaxis.cex = 1.2,
	yaxis.cex = 1.2,
	xaxis.lab = gsub('\\.','-',levels(results$Group)),
	ylimits = c(-0.005,0.3),
	yat = seq(0,0.3,0.05),
	add.text = TRUE,
	text.labels = c('***','***','*'),
	text.x = c(1.5, 2.5, 3.5),
	text.y = c(0.29,0.280,0.27),
	text.cex = 1.2,
	add.rectangle = TRUE,
	xleft.rectangle = c(1, 1, 1),
	xright.rectangle = c(2, 3, 4),
	ytop.rectangle = c(0.285, 0.275, 0.265),
	ybottom.rectangle = c(0.284, 0.274, 0.264),
	col.rectangle = 'black',
	style = 'Nature',
	filename = generate.filename('EVOLVE_ctDNA','_ratio_long_boxplot', 'png')
	);

# group by type
main.plot.objects <- list();
insert.plot.objects <- list();

plot.data <- merge(
	results[,c('Sample','Group')],
	reshape(tumour.summary,direction = 'long', varying = list(2:ncol(tumour.summary)), times = names(tumour.data), v.names = 'Density'),
	by.x = 'Sample',
	by.y = 'time'
	);

tmp <- data.frame(
	Sample = rep('NORMAL', nrow(normal.summary)),
	Group = rep('NORMAL', nrow(normal.summary)),
	Breakpoints = normal.summary$Breakpoints,
	Density = normal.summary$Normal,
	id = 1:nrow(normal.summary)
	);

plot.data <- droplevels(rbind(plot.data, tmp));

for (group in c('baseline','on.trial','EOT')) {

	legend.lab <- if (group == 'baseline') {
		'baseline' } else if (group == 'EOT') {
		'end-of-treatment' } else { 'on-trial' }

	main.plot.objects[[group]] <- create.scatterplot(
		Density ~ Breakpoints,
		plot.data[which(plot.data$Group %in% c(group,'NORMAL')),],
		groups = droplevels(plot.data[which(plot.data$Group %in% c(group,'NORMAL')),]$Sample),
		type = 'l',
		col = c(
			gray.colors(length(which(results$Group == group)), start = 0.4, end = 0.9, rev = TRUE),
			'red'
			),
		lwd = 1,
		xlimits = c(0,600),
		xat = seq(0,300,100),
		ylimits = c(0,0.041),
		yat = seq(0,0.04,0.01),
		xlab.label = if (group == 'on.trial') { 'Fragment Length' } else { NULL }, 
		ylab.label = if (group == 'baseline') { 'Density' } else { NULL },
		xaxis.cex = 1.2,
		yaxis.cex = 1.2,
		xlab.cex = 1.5,
		ylab.cex = 1.5,
		xaxis.fontface = 'plain',
		yaxis.fontface = 'plain',
		xaxis.tck = c(1,0),
		yaxis.tck = c(1,0),
		abline.v = 167,
		abline.lty = 2,
		abline.col = 'black',
		abline.lwd = 0.5,
		style = 'Nature',
		key = list(corner = c(0,1), background = 'white', text = list(lab = legend.lab, cex = 1.2))
		);

	insert.plot.objects[[group]] <- create.scatterplot(
		Density ~ Breakpoints,
		plot.data[which(plot.data$Group %in% c(group,'NORMAL')),],
		groups = droplevels(plot.data[which(plot.data$Group %in% c(group,'NORMAL')),]$Sample),
		type = 'l',
		col = c(
			gray.colors(length(which(results$Group == group)), start = 0.5, end = 0.9, rev = TRUE),
			'red'
			),
		lwd = 1,
		xlimits = c(50,150),
		xat = seq(50,150,50),
		ylimits = c(0,0.01),
		yat = seq(0,0.01,0.01),
		xlab.label = NULL,
		ylab.label = NULL,
		xaxis.cex = 1,
		yaxis.cex = 1,
		xaxis.fontface = 'plain',
		yaxis.fontface = 'plain',
		xaxis.tck = c(0.5,0),
		yaxis.tck = c(0.5,0),
		style = 'Nature'
		);
	}

combined.main.plot <- create.multipanelplot(
	plot.objects = main.plot.objects,
	layout.width = 3,
	layout.height = 1,
	plot.objects.widths = c(1.1,1,1),
	left.legend.padding = 0,
	right.legend.padding = 0,
	bottom.legend.padding = 0,
	top.legend.padding = 0,
	ylab.axis.padding = 3,
	xlab.axis.padding = 3,
	legend = list(
		inside = list(
			fun = draw.key,
			x = 0.15, y = 0.06,
			args = list(key = list(
				lines = list(col = c('grey50','red'), size = 2, lwd = 1.5),
				text = list(lab = c('EVOLVE','median normal'), cex = 1)
				))
			)
		)
	);

write.plot(
	combined.main.plot,
	additional.trellis.objects = insert.plot.objects,
	additional.trellis.locations = list(
		xleft = c(0.17,0.49,0.805),
		xright = c(0.37,0.69,1),
		ytop = c(0.9,0.9,0.9),
		ybottom = c(0.3,0.3,0.3)
		),
	height = 5,
	width = 12,
	resolution = 800,
	filename = generate.filename('EVOLVE_ctDNA','fragment_size_distribution_groups','png')
	);
	
main.plot <- create.scatterplot(
	Density ~ Breakpoints,
	plot.data,
	groups = plot.data$Sample,
	type = 'l',
	col = c(gray.colors(length(tumour.data), start = 0.4, end = 0.9, rev = TRUE), 'red'),
	lwd = 1,
	xlimits = c(50,350),
	xat = seq(50,350,50),
	ylimits = c(0,0.041),
	yat = seq(0,0.04,0.01),
	xlab.label = 'Fragment Length',
	ylab.label = 'Density',
	xaxis.cex = 1.2,
	yaxis.cex = 1.2,
	xlab.cex = 1.5,
	ylab.cex = 1.5,
	xaxis.fontface = 'plain',
	yaxis.fontface = 'plain',
	xaxis.tck = c(1,0),
	yaxis.tck = c(1,0),
	abline.v = 167,
	abline.lty = 2,
	abline.col = 'black',
	abline.lwd = 0.5,
	style = 'Nature',
	right.padding = 50,
	legend = list(
		inside = list(
			fun = draw.key,
			x = 1.35, y = 0.01,
			args = list(key = list(
				lines = list(col = c('grey50','red'), size = 3, lwd = 2),
				text = list(lab = c('EVOLVE','median normal'), cex = 1.2)
				))
			)
		)
	);

insert.short <- create.scatterplot(
	Density ~ Breakpoints,
	plot.data,
	groups = plot.data$Sample,
	type = 'l',
	col = c(gray.colors(length(tumour.data), start = 0.5, end = 0.9, rev = TRUE), 'red'),
	lwd = 1,
	xlimits = c(50,150),
	xat = seq(50,150,50),
	ylimits = c(0,0.015),
	yat = seq(0,0.015,0.005),
	xlab.label = NULL,
	ylab.label = NULL,
	xaxis.cex = 1,
	yaxis.cex = 1,
	xaxis.fontface = 'plain',
	yaxis.fontface = 'plain',
	xaxis.tck = c(0.5,0),
	yaxis.tck = c(0.5,0),
	style = 'Nature'
	);

insert.long <- create.scatterplot(
	Density ~ Breakpoints,
	plot.data,
	groups = plot.data$Sample,
	type = 'l',
	col = c(gray.colors(length(tumour.data), start = 0.5, end = 0.9, rev = TRUE), 'red'),
	lwd = 1,
	xlimits = c(250,320),
	xat = c(250,300,320),
	ylimits = c(0,0.003),
	yat = seq(0,0.005,0.001),
	xlab.label = NULL,
	ylab.label = NULL,
	xaxis.cex = 1,
	yaxis.cex = 1,
	xaxis.fontface = 'plain',
	yaxis.fontface = 'plain',
	xaxis.tck = c(0.5,0),
	yaxis.tck = c(0.5,0),
	style = 'Nature'
	);

write.plot(
	main.plot,
	additional.trellis.objects = list(insert.short, insert.long),
	additional.trellis.locations = list(
		xleft = c(0.35, 0.5),
		xright = c(0.98, 0.98),
		ytop = c(0.98, 0.5),
		ybottom = c(0.45, 0.17)
		),
	height = 5,
	width = 9,
	resolution = 1600,
	filename = generate.filename('EVOLVE_ctDNA','fragment_size_distribution_all','png')
	);

### SAVE SESSION INFO ##############################################################################
save.session.profile(generate.filename('FragmentSizeSummary','SessionProfile','txt'));
