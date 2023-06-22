### examine_fragment_sizes.R #######################################################################
# Visualize differences in global fragment sizes between tumour/normal samples, and 
# mutation-specific fragment sizes within tumour cohorts.

### PREPARE SESSION ################################################################################
# import libraries
library(BoutrosLab.plotting.general);
library(survival);

# load in helper functions (used for survival analyses)
source(paste0(getwd(),'/helper_functions/session.functions.R'));
source(paste0(getwd(),'/helper_functions/fit.coxmodel.R'));
source(paste0(getwd(),'/helper_functions/ph.fails.R'));
source(paste0(getwd(),'/helper_functions/calculate.number.at.risk.R'));
source(paste0(getwd(),'/helper_functions/create.km.plot.R'));

# get clinical covariates
load(paste0(getwd(),'/data/EVOLVE_ctDNA__clinical_timeline.RData'));

# get global fragment data 
load(paste0(getwd(),'/data/EVOLVE_ctDNA__fragmentation_summary.RData'));

# get mutation-specific fragment data 
load(paste0(getwd(),'/data/EVOLVE_ctDNA__mutation_frgment_size.RData'));

### FORMAT DATA ####################################################################################
# add annotations to global summary data
results <- merge(
	fragment.summary,
	timeline[,c('Sample','Group','Batch')],
	by = 'Sample',
	all.x = TRUE
	);

results[grep('HC', results$Sample),]$Group <- 'normal';
results$Group <- factor(results$Group, levels = c('normal', 'baseline','on.trial','EOT'));

# format survival data
survival.data <- merge(
	merge(
		timeline[which(timeline$Group == 'baseline'),c('Patient','Sample','Group')],
		clinical,
		by.x = 'Patient',
		by.y = 'Patient.ID'
		),
	ref.v.alt.results[,c('Sample','p.value.all')],
	by = 'Sample'
	);

# organize groups
groups <- rep(0, nrow(survival.data));
groups[which(survival.data$p.value.all < 0.05)] <- 1;

# collect survival stats
survtime <- survival.data$DFS.Months;
survstat <- as.numeric(survival.data$DFS);
survobj <- Surv(survtime, survstat);

# get some comparison metrics (global)
normal.idx <- which(results$Group == 'normal');
baseline.idx <- which(results$Group == 'baseline');
on.trial.idx <- which(results$Group == 'on.trial');
eot.idx <- which(results$Group == 'EOT');

p.short.baseline.v.normal <- wilcox.test(results[baseline.idx,]$Prop.short, results[normal.idx,]$Prop.short)$p.value;
p.short.ontrial.v.normal <- wilcox.test(results[on.trial.idx,]$Prop.short, results[normal.idx,]$Prop.short)$p.value;
p.short.eot.v.normal <- wilcox.test(results[eot.idx,]$Prop.short, results[normal.idx,]$Prop.short)$p.value;

p.long.baseline.v.normal <- wilcox.test(results[baseline.idx,]$Prop.long, results[normal.idx,]$Prop.long)$p.value;
p.long.ontrial.v.normal <- wilcox.test(results[on.trial.idx,]$Prop.long, results[normal.idx,]$Prop.long)$p.value;
p.long.eot.v.normal <- wilcox.test(results[eot.idx,]$Prop.long, results[normal.idx,]$Prop.long)$p.value;

### VISUALIZE DATA #################################################################################
# start with boxplots to compare proportion of short/long between tumour/normal
simplify.pvalue <- function(x) {
	if (x < 0.001) { '***' } else if (x < 0.01) { '**' } else if (x < 0.1) { '*' } else { 'ns' }
	}

figure3b.plot <- create.boxplot(
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
	text.labels = sapply(
		c(p.short.baseline.v.normal,p.short.ontrial.v.normal,p.short.eot.v.normal),
		simplify.pvalue
		),
	text.x = c(1.5, 2.5, 3.5),
	text.y = c(0.355,0.346,0.335),
	text.cex = c(1.2,0.8,1.2),
	add.rectangle = TRUE,
	xleft.rectangle = c(1, 1, 1),
	xright.rectangle = c(2, 3, 4),
	ytop.rectangle = c(0.351, 0.341, 0.331),
	ybottom.rectangle = c(0.350, 0.340, 0.330),
	col.rectangle = 'black',
	style = 'Nature'
	);

write.plot(
	figure3b.plot,
	resolution = 200,
	filename = generate.filename('EVOLVE_ctDNA','proportion150__Figure3B', 'png')
	);

supp.figure4b.plot <- create.boxplot(
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
	text.labels = sapply(
		c(p.long.baseline.v.normal,p.long.ontrial.v.normal,p.long.eot.v.normal),
		simplify.pvalue
		),
	text.x = c(1.5, 2.5, 3.5),
	text.y = c(0.135,0.130,0.125),
	text.cex = 1.2,
	add.rectangle = TRUE,
	xleft.rectangle = c(1, 1, 1),
	xright.rectangle = c(2, 3, 4),
	ytop.rectangle = c(0.133, 0.128, 0.123),
	ybottom.rectangle = c(0.1325, 0.1275, 0.1225),
	col.rectangle = 'black',
	style = 'Nature'
	);

write.plot(
	supp.figure4b.plot,
	resolution = 200,
	filename = generate.filename('EVOLVE_ctDNA','proportionLong__SuppFigure4B', 'png')
	);

# plot density functions, organized by timepoint
main.plot.objects <- list();
insert.plot.objects <- list();

plot.data <- merge(
	results[,c('Sample','Group')],
	reshape(
		tumour.summary,
		direction = 'long',
		varying = list(2:ncol(tumour.summary)),
		times = names(tumour.summary)[2:ncol(tumour.summary)],
		v.names = 'Density'
		),
	by.x = 'Sample',
	by.y = 'time'
	);

tmp <- data.frame(
	Sample = rep('NORMAL', nrow(normal.summary)),
	Group = rep('NORMAL', nrow(normal.summary)),
	Breakpoints = normal.summary$Breakpoints,
	Density = normal.summary$Median,
	id = 1:nrow(normal.summary)
	);

plot.data <- droplevels(rbind(plot.data, tmp));
plot.data$Sample <- as.factor(plot.data$Sample);

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

supp.figure4a.plot <- create.multipanelplot(
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
	supp.figure4a.plot,
	additional.trellis.objects = insert.plot.objects,
	additional.trellis.locations = list(
		xleft = c(0.17,0.49,0.805),
		xright = c(0.37,0.69,1),
		ytop = c(0.9,0.9,0.9),
		ybottom = c(0.3,0.3,0.3)
		),
	height = 5,
	width = 12,
	resolution = 200,
	filename = generate.filename('EVOLVE_ctDNA','fragment_size_distribution__SuppFigure4A', 'png')
	);

# plot combined density functions
main.plot <- create.scatterplot(
	Density ~ Breakpoints,
	plot.data,
	groups = plot.data$Sample,
	type = 'l',
	col = c(gray.colors(ncol(tumour.summary)-1, start = 0.4, end = 0.9, rev = TRUE), 'red'),
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
	col = c(gray.colors(ncol(tumour.summary)-1, start = 0.5, end = 0.9, rev = TRUE), 'red'),
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
	col = c(gray.colors(ncol(tumour.summary)-1, start = 0.5, end = 0.9, rev = TRUE), 'red'),
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
	resolution = 200,
	filename = generate.filename('EVOLVE_ctDNA','fragment_size_distribution__Figure3A', 'png')
	);

# plot mutation-specific fragment summary
plot.objects <- list();

for (group in c('baseline','on.trial','EOT')) {

	smp.set <- unique(mutation.data[which(mutation.data$Group == group),]$Sample);

	plot.data <- rbind(
		data.frame(
			BASE = rep('REF', length(unlist(ref.data[smp.set]))),
			Size = unlist(ref.data[smp.set])
			),
		data.frame(
			BASE = rep('ALT', length(unlist(alt.data[smp.set]))),
			Size = unlist(alt.data[smp.set])
			)
		);

	plot.data$BASE <- factor(plot.data$BASE, levels = c('REF','ALT'));

	stat.text <- display.statistical.result(
		wilcox.test(unlist(ref.data[smp.set]), unlist(alt.data[smp.set]))$p.value,
		statistic.type = 'p',
		symbol = ' = '
		);

	plot.objects[[group]] <- create.violinplot(
		Size ~ BASE,
		plot.data,
		main = if (group == 'EOT') { 'end-of-treatment' } else if (group == 'baseline') {
			'   baseline' } else { gsub('\\.','-',group) },
		main.x = 0.65,
		main.y = 2,
		main.cex = 1.2,
		xaxis.cex = 1.2,
		yaxis.cex = 1.2,
		ylab.cex = 1.5,
		xlab.label = NULL,
		ylab.label = if (group == 'baseline') { 'Fragment Size' } else { NULL },
		ylimits = c(0,400),
		yat = c(seq(0,400,100),167,334),
		col = 'grey90',
		add.rectangle = TRUE,
		xleft.rectangle = c(0,0),
		xright.rectangle = c(3,3),
		ybottom.rectangle = c(166.5,333.5),
		ytop.rectangle = c(167.5,334.5),
		col.rectangle = 'grey70',
		legend = list(
			inside = list(fun = draw.key, args = list(key = list(text = list(lab = stat.text, cex = 0.9))), x = 0.15, y = 1.12)
			),
		style = 'Nature'
		);
	}

figure3c.plot <- create.multipanelplot(
	plot.objects = plot.objects,
	layout.width = 3,
	layout.height = 1,
	plot.objects.widths = c(1.1,1,1),
	left.legend.padding = 0,
	right.legend.padding = 0,
	bottom.legend.padding = 0,
	top.legend.padding = 1,
	ylab.axis.padding = 3,
	style = 'Nature'
	);

write.plot(
	figure3c.plot,
	size.units = 'in',
	height = 4,
	width = 8,
	resolution = 200,
	filename = generate.filename('EVOLVE_ctDNA','fragment_size_distribution__Figure3C', 'png')
	);

# make a km plot showing DFS by fragment size (difference or no difference in REF/ALT size at mutation sites)
labels <- c('no difference in FS','difference (p < 0.05)');
alt.labels <- c('no difference', '     difference');

create.km.plot(
	survival.object = survobj,
	patient.groups = groups,
	xlab.label = 'Time (months)',
	ylab.label = 'DFS Proportion',
	ylab.axis.padding = 5,
	ylab.cex = 2,
	xlab.cex = 2,
	yaxis.cex = 1.5,
	xaxis.cex = 1.5,
	line.colours = c('turquoise3','violetred3'),
	risk.labels = alt.labels,
	risk.label.pos = -3.6,
	left.padding = 4,
	key.groups.labels = labels,
	key.groups.cex = 1.5,
	key.stats.corner = c(1,1),
	key.stats.x.pos = 1,
	key.stats.y.pos = 1,
	key.groups.corner = c(1,1),
	key.groups.x.pos = 1,
	key.groups.y.pos = 0.85,
	statistical.method = 'cox',
	explicit.HR.label = FALSE,
	filename = generate.filename('EVOLVE_ctDNA', 'dfs_by_fragment_size_pvalue__Figure3D', 'png'),
	style = 'Nature'
	);

### SAVE SESSION INFO ##############################################################################
save.session.profile(generate.filename('Figure3','SessionProfile','txt'));
