### plot_coverage.R ################################################################################
# create summary plots

### PREAMBLE #######################################################################################
library(BoutrosLab.plotting.general);

source('/cluster/home/sprokope/git/analysis/helper_functions/session.functions.R');

working.dir <- '/cluster/projects/ovgroup/projects/OV_Superset/EVOLVE/ctDNA/pipeline_suite/BAMQC/Coverage';

setwd(working.dir);

### READ DATA ######################################################################################
clinical <- read.delim('../../configs/2022-09-06_sample_info_with_batch.txt');

metric.data <- read.delim('2022-08-18_EVOLVE__raw_Coverage_summary.tsv');
cov.data <- read.delim('2022-08-18_EVOLVE__raw_Coverage_statistics.tsv', row.names = 1);

### FORMAT DATA ####################################################################################
# format sample names to match clinical
metric.data$X <- gsub('_raw','',metric.data$X);

# add clinical to metric data
master.matrix <- merge(
	clinical,
	metric.data,
	by.x = 'Sample',
	by.y = 'X'
	);

master.matrix$Batch <- factor(master.matrix$Batch, levels = c(1,2), labels = c('2021','2022'));
master.matrix$Group <- factor(master.matrix$Group, levels = c('baseline','on.trial','EOT'));

# format coverage data
smp.idx <- grep('EVO',colnames(cov.data));

total.coverage <- apply(cov.data[,smp.idx],2,sum);
tmp <- cov.data[,smp.idx]/total.coverage;

cov.list <- list();
plot.depths <- c(0, 1:10, seq(15,50,5), seq(60,200,10), seq(250,500,50));

for (smp in colnames(cov.data)[smp.idx]) {

	cov.list[[smp]] <- data.frame(
		Sample = rep(smp, length(plot.depths)),
		Depth = plot.depths,
		Proportion = NA
		);

	for (i in 1:length(plot.depths)) {

		depth <- plot.depths[i];
		cov.list[[smp]][i,]$Proportion <- sum(tmp[(depth+1):nrow(tmp), smp]);

		}
	}

plot.data <- do.call(rbind, cov.list);
rownames(plot.data) <- 1:nrow(plot.data);
plot.data$Sample <- gsub('\\.','-',gsub('_raw','',plot.data$Sample));

plot.data <- merge(
	clinical,
	plot.data,
	by = 'Sample'
	);

plot.data$Batch <- factor(plot.data$Batch, levels = c(1,2));

### PLOT DATA ######################################################################################
# plot total coverage
plot.objects <- list();

plot.objects[[1]] <- create.boxplot(
	total/(10**9) ~ Group,
	droplevels(master.matrix[which(master.matrix$Batch == '2021'),]),
	add.stripplot = TRUE,
	points.cex = 1,
	points.alpha = 0.8,
	ylimits = c(0,8),
	yat = seq(0,8,2),
	xlab.cex = 1.2,
	ylab.cex = 1.2,
	xlab.label = '2021',
	ylab.label = expression('total coverage (x 10'^9*')'),
	xaxis.cex = 1,
	yaxis.cex = 1,
	xaxis.tck = c(1,0),
	yaxis.tck = c(1,0),
	style = 'Nature'
	);

plot.objects[[2]] <- create.boxplot(
	total/(10**9) ~ Group,
	droplevels(master.matrix[which(master.matrix$Batch == '2022'),]),
	add.stripplot = TRUE,
	points.cex = 1,
	points.alpha = 0.8,
	ylimits = c(0,8),
	yat = seq(0,8,2),
	yaxis.lab = rep('', 5),
	xlab.label = '2022',
	ylab.label = NULL,
	xlab.cex = 1.2,
	xaxis.cex = 1,
	yaxis.cex = 1,
	xaxis.tck = c(1,0),
	yaxis.tck = 0,
	style = 'Nature'
	);

create.multipanelplot(
	plot.objects = plot.objects,
	layout.width = 2,
	layout.height = 1,
	plot.objects.widths = c(1,0.4),
	left.legend.padding = 0,
	right.legend.padding = 0,
	bottom.legend.padding = 0,
	top.legend.padding = 0,
	ylab.axis.padding = 1,
	xlab.axis.padding = 2,
	file = generate.filename('EVOLVE_','raw_total_coverage','png'),
	width = 6,
	height = 6,
	resolution = 800
	);

main.plot <- create.scatterplot(
	Proportion ~ Depth,
	plot.data,
	groups = plot.data$Sample,
	col = c('grey30','grey70')[match(unique(plot.data[,c('Sample','Batch')])$Batch, c(1,2))],
	type = 'l',
	xlimits = c(0,500),
	ylimits = c(0.7,1),
	xat = seq(0,500,100),
	yat = seq(0.7,1,0.05),
	xlab.cex = 1.5,
	ylab.cex = 1.5,
	xlab.label = 'depth of coverage',
	ylab.label = 'proportion',
	xaxis.cex = 1,
	yaxis.cex = 1,
	xaxis.tck = c(1,0),
	yaxis.tck = c(1,0),
	xaxis.col = c('black','black','red','black','black','black'),
	add.rectangle = TRUE,
	xleft.rectangle = 200,
	xright.rectangle = 201,
	ytop.rectangle = 0.99,
	ybottom.rectangle = 0.79,
	col.rectangle = 'red',
	right.padding = 25,
	style = 'Nature',
	key = list(
		lines = list(col = c('grey30','grey70'), lwd = 2, lty = 1),
		text = list(lab = c('2021','2022'), cex = 1),
		x = 0.05,
		y = 0.15
		)
	);

insert.plot <- create.boxplot(
	Proportion ~ Batch,
	plot.data[which(plot.data$Depth == 200),],
	add.stripplot = TRUE,
	points.col = c('skyblue','seagreen','pink')[match(plot.data[which(plot.data$Depth == 200),]$Group, c('baseline','on.trial','EOT'))],
	points.alpha = 0.5,
	points.cex = 0.8,
	col = 'transparent',
	ylimits = c(0.8,1),
	yat = seq(0.8,1,0.1),
	ylab.cex = 0.8,
	xlab.label = NULL,
	ylab.label = 'proportion',
	xaxis.lab = c('2021','2022'),
	xaxis.cex = 0.8,
	yaxis.cex = 0.8,
	xaxis.tck = c(1,0),
	yaxis.tck = c(1,0),
	style = 'Nature',
	key = list(
		points = list(col = c('skyblue','seagreen','pink'), cex = 1, alpha = 0.5, pch = 19),
		text = list(lab = c('baseline','on-trial','end-of-treatment'), cex = 0.8),
		x = 0,
		y = -0.2
		)
	);

write.plot(
	main.plot,
	additional.trellis.objects = list(insert.plot),
	additional.trellis.locations = list(xleft = 0.6, xright = 1, ytop = 0.9, ybottom = 0.3),
	file = generate.filename('EVOLVE_','raw_cumulative_coverage','png'),
	height = 6,
	width = 6,
	resolution = 800
	);


