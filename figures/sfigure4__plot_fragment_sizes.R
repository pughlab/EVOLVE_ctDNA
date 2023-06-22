### sfigure4__plot_fragment_sizes.R ################################################################
# Visualize differences in global fragment sizes between tumour/normal samples, and 
# mutation-specific fragment sizes within tumour cohorts.

### PREPARE SESSION ################################################################################
# import libraries
library(BoutrosLab.plotting.general);

# load in helper functions
source(paste0(getwd(),'/helper_functions/session.functions.R'));

# get clinical covariates
load(paste0(getwd(),'/../data/EVOLVE_ctDNA__clinical_timeline.RData'));

# get fragment data 
load(paste0(getwd(),'/../data/EVOLVE_ctDNA__fragmentation_summary.RData'));

### FORMAT DATA ####################################################################################
# add annotations to global summary data
results <- merge(
	global.fragment.summary,
	timeline[,c('Sample','Group','Batch')],
	by = 'Sample',
	all.x = TRUE
	);

results[grep('HC', results$Sample),]$Group <- 'normal';
results$Group <- factor(results$Group, levels = c('normal', 'baseline','on.trial','EOT'));

# get some comparison metrics (global)
normal.idx <- which(results$Group == 'normal');
baseline.idx <- which(results$Group == 'baseline');
on.trial.idx <- which(results$Group == 'on.trial');
eot.idx <- which(results$Group == 'EOT');

p.long.baseline.v.normal <- wilcox.test(results[baseline.idx,]$Prop.long, results[normal.idx,]$Prop.long)$p.value;
p.long.ontrial.v.normal <- wilcox.test(results[on.trial.idx,]$Prop.long, results[normal.idx,]$Prop.long)$p.value;
p.long.eot.v.normal <- wilcox.test(results[eot.idx,]$Prop.long, results[normal.idx,]$Prop.long)$p.value;

### VISUALIZE DATA #################################################################################
# start with boxplots to compare proportion of short/long between tumour/normal
simplify.pvalue <- function(x) {
	if (x < 0.001) { '***' } else if (x < 0.01) { '**' } else if (x < 0.1) { '*' } else { 'ns' }
	}

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
		global.tumour.summary,
		direction = 'long',
		varying = list(2:ncol(global.tumour.summary)),
		times = names(global.tumour.summary)[2:ncol(global.tumour.summary)],
		v.names = 'Density'
		),
	by.x = 'Sample',
	by.y = 'time'
	);

tmp <- data.frame(
	Sample = rep('NORMAL', nrow(global.normal.summary)),
	Group = rep('NORMAL', nrow(global.normal.summary)),
	Breakpoints = global.normal.summary$Breakpoints,
	Density = global.normal.summary$Median,
	id = 1:nrow(global.normal.summary)
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

### SAVE SESSION INFO ##############################################################################
save.session.profile(generate.filename('SuppFig4parts','SessionProfile','txt'));
