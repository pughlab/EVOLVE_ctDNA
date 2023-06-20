### compare_tumour_content_estimates.R ##############################################################
# Compare tumour content estimation methods by 
#	- TP53 mutation with bam-readcount
#	- maximum somatic VAF
#	- proportion of short (<150bp) fragments

### PREPARE SESSION ################################################################################
# import libraries
library(BoutrosLab.plotting.general);

source('/cluster/home/sprokope/git/analysis/helper_functions/session.functions.R');

input.dir <- '/cluster/projects/ovgroup/projects/OV_Superset/EVOLVE/ctDNA/pipeline_suite';
output.dir <- 'Ensemble_calls/tumour_content';

setwd(input.dir);

### READ DATA ######################################################################################
# get clinical data (for ctDNA)
timeline <- read.delim('configs/2022-09-06_sample_info_with_batch.txt');
timeline <- unique(timeline[order(timeline$Patient, timeline$Timepoint),]);

# read in estimated tumour content
purity <- read.delim('fragment_size/2022-09-07_EVOLVE_ctDNA__fragment_size_and_purity_estimates.tsv', stringsAsFactors = FALSE);

# move to output directory
setwd(output.dir);

### FORMAT DATA ####################################################################################
# format purity data
purity.formatted <- merge(
	timeline,
	purity[,c('Patient','Sample','Prop.short','Prop.long','Purity','Max.VAF')],
	by.x = c('Patient','Sample'),
	by.y = c('Patient','Sample'),
	all.x = TRUE
	);

# compare estimates
cor(purity.formatted[,c('Purity','Max.VAF','Prop.short','Prop.long')], use = 'pairwise');

cor(purity.formatted[which(purity.formatted$Timepoint == 0),c('Purity','Max.VAF','Prop.short','Prop.long')], use = 'pairwise');
cor(purity.formatted[which(purity.formatted$Timepoint == 2),c('Purity','Max.VAF','Prop.short','Prop.long')], use = 'pairwise');
cor(purity.formatted[which(purity.formatted$Timepoint == 3),c('Purity','Max.VAF','Prop.short','Prop.long')], use = 'pairwise');

### PLOT DATA ######################################################################################
# next, let's compare estimates
plot.data <- purity.formatted[,c('Sample','Purity','Max.VAF','Prop.short','Prop.long')];
plot.data[is.na(plot.data$Purity),]$Purity <- -0.07;
plot.data[is.na(plot.data$Max.VAF),]$Max.VAF <- -0.07;

combinations.to.plot <- data.frame(t(combn(c('Purity','Max.VAF','Prop.short','Prop.long'),2)));
combinations.to.plot$Slope <- NA;
combinations.to.plot$Intercept <- NA;

plot.objects <- list();

for (i in 1:nrow(combinations.to.plot)) {
	metric.y <- as.character(combinations.to.plot[i,1]);
	metric.x <- as.character(combinations.to.plot[i,2]);
	best.fit.line <- lm(get(metric.y) ~ get(metric.x), plot.data[which(plot.data[,metric.x] >= 0 & plot.data[,metric.y] >= 0),])$coef;
	combinations.to.plot[i,]$Slope <- as.numeric(best.fit.line[2]);
	combinations.to.plot[i,]$Intercept <- as.numeric(best.fit.line[1]);
	}

# fill in manually because the panel function takes the most recent values of i
line.functions <- list();
line.functions[[1]] <- function(x) { x * 0.9898036 + 0.002473406 }
line.functions[[2]] <- function(x) { x * 0.8364555 + -0.079452183 }
line.functions[[3]] <- function(x) { x * 1.6137545 + 0.013863082 }
line.functions[[4]] <- function(x) { x * 0.6183754 + -0.030518511 }
line.functions[[5]] <- function(x) { x * 1.8235598 + 0.016606749 }
line.functions[[6]] <- function(x) { x * -0.6181525 + 0.205372539 }

# plot each pair
for (i in 1:nrow(combinations.to.plot)) {

	metric.y <- as.character(combinations.to.plot[i,1]);
	metric.x <- as.character(combinations.to.plot[i,2]);

	label.y <- if (metric.y == 'Purity') { expression(italic('TP53')*' VAF')
		} else if (metric.y == 'Max.VAF') { 'Maximum VAF'
		} else if (metric.y == 'Prop.short') { 'Proportion fragments <150bp'
		} else { 'Proportion fragments 250-320bp' }

	label.x <- if (metric.x == 'Purity') { expression(italic('TP53')*' VAF')
		} else if (metric.x == 'Max.VAF') { 'Maximum VAF'
		} else if (metric.x == 'Prop.short') { 'Proportion fragments <150bp'
		} else { 'Proportion fragments 250-320bp' }

	# create plot
	plot.objects[[i]] <- create.scatterplot(
		get(metric.y) ~ get(metric.x),
		plot.data,
		alpha = 0.8,
		ylimits = c(-0.1,1),
		yat = c(-0.07, seq(0,1,0.2)),
		yaxis.lab = c('NA', '0','0.2','0.4','0.6','0.8','1.0'),
		xlimits = c(-0.1,1),
		xat = c(-0.07, seq(0,1,0.2)),
		xaxis.lab = c('NA', '0','0.2','0.4','0.6','0.8','1.0'),
		yaxis.tck = c(1,0),
		xaxis.tck = c(1,0),
		ylab.label = label.y,
		xlab.label = label.x,
		ylab.cex = 1.2,
		xlab.cex = 1.2,
		yaxis.cex = 1,
		xaxis.cex = 1,
		abline.h = 0,
		abline.v = 0,
		abline.lty = 2,
		abline.lwd = 1,
		add.rectangle = TRUE,
		xleft.rectangle = c(-0.1,0),
		xright.rectangle = c(0,1),
		ytop.rectangle = c(1,0),
		ybottom.rectangle = c(0,-0.1),
		col.rectangle = c('grey50','grey50'),
		alpha.rectangle = 0.5,
		style = 'Nature',
		key = get.corr.key(
			x = plot.data[,metric.x],
			y = plot.data[,metric.y],
			label.items = c('pearson'),
			x.pos = 0.16,
			y.pos = 0.95,
			key.cex = 1.5
			),
		add.curves = TRUE,
		curves.from = 0,
		curves.to = max(0.3, max(plot.data[,metric.x])),
		curves.lwd = 1,
		curves.exprs = list(line.functions[[i]])
		);
	}

# fill in manually because the panel function takes the most recent values of i
#line.functions <- list();
#line.functions[[1]] <- function(x) { x * 0.9898036 + 0.002473406 }
#line.functions[[2]] <- function(x) { x * 0.8364555 + -0.079452183 }
#line.functions[[3]] <- function(x) { x * 1.6137545 + 0.013863082 }
#line.functions[[4]] <- function(x) { x * 0.6183754 + -0.030518511 }
#line.functions[[5]] <- function(x) { x * 1.8235598 + 0.016606749 }
#line.functions[[6]] <- function(x) { x * -0.6181525 + 0.205372539 }

# combine them
combined.plot <- create.multipanelplot(
	plot.objects = plot.objects,
	layout.width = 3,
	layout.height = 2,
	plot.objects.widths = c(1,1,1),
	plot.objects.heights = c(1,1),
	x.spacing = 1,
	y.spacing = 1,
	left.legend.padding = 0,
	right.legend.padding = 0,
	bottom.legend.padding = 0,
	top.legend.padding = 0
	);

write.plot(
	combined.plot,
	height = 8,
	width = 12,
	resolution = 800,
	filename = generate.filename('EVOLVE_ctDNA', '_compare_estimated_tumour_contents','png')
	);

### SAVE SESSION INFO ##############################################################################
save.session.profile(generate.filename('ComparePuritySummary','SessionProfile','txt'));
