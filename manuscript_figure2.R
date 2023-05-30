### compare_tumour_content_max_vaf.R ###############################################################
# Compare tumour content estimation (by TP53 mutation with bam-readcount) and by maximum VAF

### PREPARE SESSION ################################################################################
# import libraries
library(BoutrosLab.plotting.general);

source('/cluster/home/sprokope/git/analysis/helper_functions/session.functions.R');

input.dir <- '/cluster/projects/ovgroup/projects/OV_Superset/EVOLVE/ctDNA/pipeline_suite/Ensemble_calls/tumour_content';
output.dir <- '/cluster/projects/ovgroup/projects/OV_Superset/EVOLVE/ctDNA/pipeline_suite/paper_figures';

setwd(input.dir);

### READ DATA ######################################################################################
# create covariates
load('../../2022-09-06_EVOLVE_ctDNA_clinicalCovariates.RData');

# read in estimated tumour content
results <- read.delim('2022-09-06_EVOLVE_ctDNA__estimated_tumour_content.tsv', stringsAsFactors = FALSE);

setwd(output.dir);

### FORMAT DATA ####################################################################################
# format data for scatterplot (Figure 2A)
plot.data <- results[order(results$Patient.ID, results$Timepoint),];
plot.data$Group.ctdna <- factor(
	plot.data$Group.ctdna,
	levels = c('baseline','on.trial','EOT'),
	labels = c('baseline','on-trial','end-of-treatment')
	);

# format data for barplot (Figure 2B)
results.wide <- reshape(
	results[,c('Patient.ID','Group.ctdna','Final')],
	direction = 'wide',
	idvar = 'Patient.ID',
	timevar = 'Group.ctdna'
	);
colnames(results.wide) <- gsub('Final\\.','',colnames(results.wide));

# calculate clearance
results.wide$delta <- apply(
	results.wide[,c('baseline','on.trial')],
	1,
	function(i) {
		i <- na.omit(i);
		if (length(i) < 2) { return(NA) } else { return( 100* ( (i[2]-i[1]) / i[1]) ) }
		}
	);

results.wide$max <- apply(
	results.wide[,c('baseline','on.trial','EOT')],
	1,
	function(i) {
		if (all(is.na(i))) { return(NA) } else { max(i, na.rm = TRUE) }
		}	
	);

# sort data
results.wide <- results.wide[order(results.wide$delta, results.wide$max, decreasing = FALSE, na.last = FALSE),];
results.wide$Order <- 1:nrow(results.wide);

# add annotations to clearance data
clearance.data <- merge(clinical[,c(1,9:11,14)], results.wide, by = 'Patient.ID');
clearance.data <- clearance.data[order(clearance.data$Order),];

### COVARIATES AND LEGENDS #########################################################################
# colour for clearance bars
point.colours <- rep('cornflowerblue', nrow(clearance.data));
point.colours[is.na(clearance.data$delta)] <- NA;
point.colours[which(clearance.data$delta > 0)] <- 'lightcoral';

covariate.colours$Response$colours[3] <- "#81A1C1";
covariate.colours$BRCA$colours[3] <- 'grey90';

# create covariates for clearance barplot
smp.covariates <- list(
	rect = list(
		col = 'transparent',
		fill = covariate.colours$Cohort$colours[match(clearance.data$Cohort.cat, covariate.colours$Cohort$levels)],
		lwd = 1
		),
	rect = list(
		col = 'transparent',
		fill = covariate.colours$BRCA$colours[match(clearance.data$BRCA.cat, covariate.colours$BRCA$levels)],
		lwd = 1
		),
	rect = list(
		col = 'transparent',
		fill = covariate.colours$Age$colours[match(clearance.data$Age.cat, covariate.colours$Age$levels)],
		lwd = 1
		),
	rect = list(
		col = 'transparent',
		fill = covariate.colours$Response$colours[match(clearance.data$Response.cat, covariate.colours$Response$levels)],
		lwd = 1
		)
	);

smp.covariate.grob <- covariates.grob(
	covariates = smp.covariates,
	ord = 1:nrow(clearance.data),
	side = 'right',
	size = 1,
	grid.col = list(col = 'white', lwd = 2),
	col.lines = 1:length(smp.covariates)
	);

smp.legends <- list(
	legend = list(
		colours = covariate.colours$Cohort$colours,
		border = 'transparent',
		labels = covariate.colours$Cohort$labels,
		title = 'Cohort'
		),
	legend = list(
		colours = covariate.colours$BRCA$colours,
		border = 'transparent',
		labels = covariate.colours$BRCA$labels,
		title = 'Germline BRCA Status'
		),
	legend = list(
		colours = covariate.colours$Age$colours,
		border = 'transparent',
		labels = covariate.colours$Age$labels,
		title = 'Age'
		),
	legend = list(
		colours = covariate.colours$Response$colours,
		border = 'transparent',
		labels = covariate.colours$Response$labels,
		title = 'Response'
		),
	legend = list(
		colours = c('lightcoral','cornflowerblue'),
		labels = c('increase','decrease'),
		border = 'transparent',
		title = 'ctDNA Change'
		)
	);

smp.legend.grob <- legend.grob(
	legends = smp.legends,
	label.cex = 1,
	title.cex = 1,
	title.fontface = 'plain',
	title.just = 'left',
	size = 2,
	layout = c(1,length(smp.legends))
	);

# create a text-key for samples with no variants
nv.key <- list(
	text = list(lab = 'no data', cex = 1)
	);

### PLOT DATA ######################################################################################
# plot clearance separately
create.barplot(
	Order ~ delta,
	clearance.data,
	col = point.colours[!is.na(point.colours)],
	ylab.label = NULL,
	xlab.label = expression(Delta * 'ctDNA'),
	xlab.cex = 1.5,
	xaxis.cex = 1.2,
	yaxis.cex = 1.2,
	xlimits = c(-120,300),
	xat = c(-100,0,100,300),
	yaxis.lab = rep('    ', nrow(clearance.data)),
	xaxis.fontface = 'plain',
	yaxis.fontface = 'plain',
	xaxis.tck = c(1,0),
	yaxis.tck = 0,
	plot.horizontal = TRUE,
	legend = list(
		left = list(fun = smp.covariate.grob)
		),
	left.padding = 8,
	style = 'Nature',
	height = 9,
	width = 4,
	filename = generate.filename('EVOLVE_ctDNA', '_estimated_tumour_clearance','png')
	);

# try a scatterplot per patient?
create.scatterplot(
	log10(Final) ~ Group.ctdna | Patient.ID,
	plot.data,
	type = 'b',
	xaxis.rot = 90,
	x.spacing = 1,
	y.spacing = 1,
	as.table = TRUE,
	xaxis.tck = c(1,0),
	yaxis.tck = c(1,0),
	xlab.label = NULL,
	ylab.label = expression('Estimated Tumour Fraction'),
	ylab.axis.padding = 2,
	yat = log10(c(0.001,0.01,0.1,1)),
	yaxis.lab = c(0.001,0.01,0.1,1),
	xaxis.cex = 1.2,
	yaxis.cex = 1.2,
	ylab.cex = 1.5,
	xaxis.fontface = 'plain',
	yaxis.fontface = 'plain',
	strip.fontface = 'plain',
	strip.cex = 0.9,
	layout = c(6,5),
	legend = list(
		inside = list(fun = draw.key, args = list(key = nv.key), x = 0.187, y = 0.5),
		inside = list(fun = draw.key, args = list(key = nv.key), x = 0.187, y = 0.09),
		inside = list(fun = draw.key, args = list(key = nv.key), x = 0.53, y = 0.09)
		),
	height = 9,
	width = 9,
	filename = generate.filename('EVOLVE_ctDNA', '_estimated_tumour_fraction','png')
	);

# combine them
per.patient.plots <- list();
plot.data$Group <- as.numeric(plot.data$Group.ctdna);

for (patient in unique(plot.data$Patient.ID)) {

	per.patient.plots[[patient]] <- create.scatterplot(
		log10(Final) ~ Group | Patient.ID,
		plot.data[which(plot.data$Patient.ID == patient),],
		type = 'b',
		col = covariate.colours$Response$colours[match(
			clinical[which(clinical$Patient.ID == patient),]$Response.cat,
			covariate.colours$Response$levels)],
		xlimits = c(0.5,3.5),
		xat = c(1,2,3),
		xaxis.lab = levels(plot.data$Group.ctdna),
		xaxis.tck = c(1,0),
		yaxis.tck = c(1,0),
		xlab.label = NULL,
		ylab.label = NULL,
		ylimits = log10(c(0.001,1)),
		yat = log10(c(0.001,0.01,0.1,1)),
		yaxis.lab = c(0.001,0.01,0.1,1),
		strip.fontface = 'plain',
		strip.cex = 0.9
		);
	}

# combine them
figure2a <- create.multiplot(
	plot.objects = c(
		per.patient.plots[c(19,21,16,14,26,28)],
		per.patient.plots[c(22,3,6,12,15,18)],
		per.patient.plots[c(23,24,25,27,29,30)],
		per.patient.plots[c(9,10,11,13,17,20)],
		per.patient.plots[c(1,2,4,5,7,8)]
		),
	plot.layout = c(6,5),
	panel.heights = rep(1,5),
	panel.widths = rep(1,6),
	ylab.label = expression('Estimated Tumour Fraction'),
	ylab.padding = 3,
	axes.lwd = 1,
	xaxis.cex = 1.2,
	yaxis.cex = 1.2,
	xaxis.rot = 90,
	xaxis.fontface = 'plain',
	yaxis.fontface = 'plain',
	xaxis.tck = c(1,0),
	yaxis.tck = c(1,0),
	print.new.legend = TRUE,
	legend = list(
		inside = list(fun = draw.key, args = list(key = nv.key), x = 0.533, y = 0.09),
		inside = list(fun = draw.key, args = list(key = nv.key), x = 0.708, y = 0.09),
		inside = list(fun = draw.key, args = list(key = nv.key), x = 0.878, y = 0.09)
		)
	);

covariate.key <- list(
	text = list(
		lab = c('Cohort','gBRCA','Age','Response'),
		cex = 1.1,
		adj = 1
		),
	padding.text = 1.8
	);

figure2b <- create.barplot(
	Order ~ delta,
	clearance.data,
	col = point.colours[!is.na(point.colours)],
	ylab.label = NULL,
	xlab.label = expression(Delta * 'ctDNA'),
	xlab.cex = 1.5,
	xaxis.cex = 1.2,
	yaxis.cex = 1.2,
	xlimits = c(-120,300),
	xat = c(-100,0,100,300),
	yaxis.lab = rep('', nrow(clearance.data)),
	xaxis.fontface = 'plain',
	yaxis.fontface = 'plain',
	xaxis.tck = c(1,0),
	yaxis.tck = 0,
	plot.horizontal = TRUE,
	legend = list(left = list(fun = smp.covariate.grob),
		inside = list(fun = draw.key(vp = viewport(angle = 90), key = covariate.key), x = -0.638, y = -0.02)),
	style = 'Nature'
	);

# combine them
create.multipanelplot(
	plot.objects = list(figure2a, figure2b),
	layout.width = 2,
	layout.height = 1,
	plot.objects.widths = c(3,1),
	x.spacing = 9,
	right.padding = 1,
	top.padding = 0,
	bottom.padding = 7,
	legend = list( inside = list(fun = smp.legend.grob, x = 0.69, y = 0.56) ),
	left.legend.padding = 0,
	right.legend.padding = 0,
	bottom.legend.padding = 0,
	top.legend.padding = 0,
	xlab.axis.padding = -15,
	ylab.axis.padding = c(3,0),
	height = 8,
	width = 14,
	resolution = 1600,
	filename = generate.filename('EVOLVE_ctDNA', '_figure2','png')
	);

### SAVE SESSION INFO ##############################################################################
save.session.profile(generate.filename('Figure2','SessionProfile','txt'));
