### plot_sample_timeline.R #########################################################################
# Plot sample timelines per patient

### FUNCTIONS ######################################################################################
# function to generate a standardized filename
generate.filename <- function(project.stem, file.core, extension, include.date = TRUE) {

	# build up the filename
	file.name <- paste(project.stem, file.core, sep = '_');
	file.name <- paste(file.name, extension, sep = '.');

	if (include.date) {
		file.name <- paste(Sys.Date(), file.name, sep = '_');
		}

	return(file.name);
	}

# function to write session profile to file
save.session.profile <- function(file.name) {

	# open the file
	sink(file = file.name, split = FALSE);

	# write memory usage to file
	cat('### MEMORY USAGE ###############################################################');
	print(proc.time());

	# write sessionInfo to file
	cat("\n### SESSION INFO ###############################################################");
	print(sessionInfo());

	# close the file
	sink();

	}

### PREPARE SESSION ################################################################################
# import libraries
library(BoutrosLab.plotting.general);

load('/Users/sprokopec/git/EVOLVE_ctDNA/data/EVOLVE_ctDNA__clinical_timeline.RData');

### FORMAT DATA ####################################################################################
# get list of patients
all.patients <- intersect(clinical$Patient.ID, timeline$Patient);

# sort data by final timepoint available (time to progression)
plot.data <- aggregate(Timepoint ~ Patient, timeline, max);

# but do so within each cohort
tmp <- merge(clinical, plot.data, by.x = 'Patient.ID', by.y = 'Patient', all.x = TRUE);
tmp <- tmp[order(tmp$Cohort.cat, -tmp$Timepoint),];
patient.order <- tmp$Patient.ID;

# ensure covariates are ordered properly
clinical$Patient.ID <- factor(clinical$Patient.ID, levels = rev(patient.order));
clinical <- clinical[order(clinical$Patient),];

timeline$Patient <- factor(timeline$Patient, levels = rev(patient.order));
timeline$Group <- as.numeric(timeline$Patient);

plot.data$Patient <- factor(plot.data$Patient, levels = rev(patient.order));
plot.data <- plot.data[order(plot.data$Patient),];
plot.data$Group <- as.numeric(plot.data$Patient);

### PLOT DATA ######################################################################################
# create covariates
smp.covariates <- list(
	rect = list(
		col = 'transparent',
		fill = covariate.colours$Cohort$colours[match(clinical$Cohort.cat, covariate.colours$Cohort$levels)],
		lwd = 1
		),
	rect = list(
		col = 'transparent',
		fill = covariate.colours$Tissue$colours[match(clinical$Oncotree.Code, covariate.colours$Tissue$levels)],
		lwd = 1
		),
	rect = list(
		col = 'transparent',
		fill = covariate.colours$BRCA$colours[match(clinical$BRCA.cat, covariate.colours$BRCA$levels)],
		lwd = 1
		),
	rect = list(
		col = 'transparent',
		fill = covariate.colours$Age$colours[match(clinical$Age.cat, covariate.colours$Age$levels)],
		lwd = 1
		),
	rect = list(
		col = 'transparent',
		fill = covariate.colours$Response$colours[match(clinical$Response, covariate.colours$Response$levels)],
		lwd = 1
		),
	rect = list(
		col = 'transparent',
		fill = covariate.colours$DFS$colours[match(clinical$DFS, covariate.colours$DFS$levels)],
		lwd = 1
		),
	rect = list(
		col = 'transparent',
		fill = covariate.colours$OS$colours[match(clinical$OS, covariate.colours$OS$levels)],
		lwd = 1
		)
	);

smp.covariate.grob <- covariates.grob(
	covariates = smp.covariates,
	ord = 1:nrow(clinical),
	side = 'right',
	size = 0.6,
	grid.col = list(col = 'white', lwd = 1),
	col.lines = 1:length(smp.covariates)
	);

covariate.key <- list(
	text = list(
		lab = c('Cohort','Tissue','gBRCA','Age','Response','DFS','OS'),
		cex = 0.8,
		adj = 1
		),
	padding.text = 0.65
	);

dot.key <- list(
	points = list(col = 'black', fill = c('black','white'), pch = c(23,23,4), cex = 0.8),
	text = list(lab = c('success','fail','not run'), cex = 0.8),
	title = 'Sequencing Status',
	cex.title = 0.8
	);

# create timeline plot
timeline.plot <- create.scatterplot(
	Group ~ Timepoint,
	timeline,
	pch = c(23,23,4)[match(timeline$Status, c(1,0,-1))],
	cex = 0.8,
	col = c('black','transparent','black')[match(timeline$Status, c(1,0,-1))],
	add.rectangle = TRUE,
	xleft.rectangle = rep(0,nrow(plot.data)),
	xright.rectangle = plot.data$Timepoint,
	ytop.rectangle = plot.data$Group + 0.15,
	ybottom.rectangle = plot.data$Group - 0.15,
	col.rectangle = 'red',
	alpha.rectangle = 0.5,
	xlimits = c(-1,30),
	xat = seq(0,30,5),
	xaxis.lab = seq(0,30,5),
	right.padding = 1,
	left.padding = 18,
	ylimits = c(0.5, max(timeline$Group)+0.5),
	yat = seq(1, max(timeline$Group),1),
	yaxis.lab = paste0(levels(timeline$Patient), '                    '),
	yaxis.cex = 0.95,
	xaxis.cex = 1,
	xaxis.tck = c(0.5,0),
	yaxis.tck = c(0,0),
	xlab.label = 'Number of Cycles',
	xlab.cex = 1.5,
	ylab.label = NULL,
	legend = list(
		inside = list(fun = smp.covariate.grob, x = -0.345, y = 0.5),
		inside = list(fun = draw.key, args = list(dot.key), x = 0.63, y = 0.8),
		inside = list(fun = smp.legend.grob, x = -1.39, y = 0.99),
		inside = list(fun = draw.key(vp = viewport(angle = 90), key = covariate.key), x = -0.34, y = 0.01)
		),
	style = 'Nature'
	);

write.plot(
	timeline.plot,
	size.units = 'in',
	width = 6.5,
	height = 6.5,
	resolution = 200,
	filename = generate.filename('EVOLVE_ctDNA', 'timeline_plot__Figure1','png')
	);

### SAVE SESSION INFO ##############################################################################
save.session.profile(generate.filename('Figure1','SessionProfile','txt'));
