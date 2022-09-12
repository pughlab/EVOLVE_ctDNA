### plot_sample_timeline.R #########################################################################
# Plot sample timelines per patient

### PREPARE SESSION ################################################################################
# import libraries
library(BoutrosLab.plotting.general);

source('/cluster/home/sprokope/git/analysis/helper_functions/session.functions.R');

setwd('/cluster/projects/ovgroup/projects/OV_Superset/EVOLVE/ctDNA/pipeline_suite/');

### READ DATA ######################################################################################
# get clinical data
clinical <- read.delim('configs/EVOLVE_clinical_data_updated2022.tsv');
clinical <- clinical[which(clinical$SAMPLE_TYPE == 'Archival'),];
clinical <- unique(clinical[,c('Patient.ID','AGE','Germline_BRCA_Status','Best.Response','COHORT','Disease.Free.Status','Oncotree.Code','Overall.Survival.Status')]);

# get timeline data
timeline <- read.delim('configs/2022-09-06_sample_info_with_batch.txt');
timeline$Status <- 1;

# add in samples that failed library prep
failed.smps <- data.frame(
	Patient = unlist(strsplit('EVO-009-008 EVO-009-012 EVO-009-021 EVO-400-002 EVO-400-004 EVO-400-007',' ')), 
	Sample = rep('', 6),
	Tumor.Type = rep('C2D1', 6),
	Timepoint = rep(2, 6),
	Group = rep('on.trial',6),
	Batch = rep(2,6),
	Status = rep(0,6)
	);

# add in samples that were never run (ctDNA)
missing.smps <- data.frame(
	Patient = c('EVO-009-025','EVO-400-001','EVO-400-006'),
	Sample = rep('', 3),
	Tumor.Type = rep(NA, 3),
	Timepoint = rep(0, 3),
	Group = rep(NA,3),
	Batch = rep(NA,3),
	Status = rep(-1,3)
	);

# combine them
timeline <- rbind(timeline, failed.smps, missing.smps);

### FORMAT DATA ####################################################################################
# get list of patients
all.patients <- intersect(clinical$Patient.ID, timeline$Patient);

clinical <- clinical[which(clinical$Patient.ID %in% all.patients),];

# format data for covariates
clinical$Age.cat <- '40 - 50';
clinical[which(clinical$AGE >= 50),]$Age.cat <- '50 - 60';
clinical[which(clinical$AGE >= 60),]$Age.cat <- '60 - 70';
clinical[which(clinical$AGE >= 70),]$Age.cat <- '70 - 80';

clinical$Response.cat <- 'SD';
clinical[which(clinical$Best.Response == 'Progressive Disease (Objective)'),]$Response.cat <- 'PD';
clinical[which(clinical$Best.Response == 'Partial Response'),]$Response.cat <- 'PR';
clinical[which(clinical$Best.Response == 'Inevaluable'),]$Response.cat <- 'IE';

clinical$Cohort.cat <- '1';
clinical[which(clinical$COHORT == '2. Platinum Resistant'),]$Cohort.cat <- '2';
clinical[which(clinical$COHORT == '3. Exploratory'),]$Cohort.cat <- '3';

clinical$DFS <- '1';
clinical[which(clinical$Disease.Free.Status == 'DiseaseFree'),]$DFS <- '0';

clinical$OS <- '1';
clinical[which(clinical$Overall.Survival.Status == 'LIVING'),]$OS <- '0';

clinical$BRCA.cat <- factor(clinical$Germline_BRCA_Status, levels = c('BRCA1','BRCA2','WT'));

# format timeline for sample ordering
timeline <- unique(timeline[order(timeline$Patient, timeline$Timepoint),]);
timeline$Sample <- gsub('_ctDNA','',timeline$Sample);

plot.data <- aggregate(Timepoint ~ Patient, timeline, max);

# plot by final timepoint
tmp <- merge(clinical, plot.data, by.x = 'Patient.ID', by.y = 'Patient', all.x = TRUE);
tmp <- tmp[order(tmp$COHORT, -tmp$Timepoint),]; # tmp$Oncotree.Code, tmp$DFS, tmp$OS, tmp$Best.Response, tmp$AGE),];
patient.order <- tmp$Patient.ID;

# ensure covariates are ordered properly
clinical$Patient.ID <- factor(clinical$Patient.ID, levels = rev(patient.order));
clinical <- clinical[order(clinical$Patient),];

timeline$Patient <- factor(timeline$Patient, levels = rev(patient.order));
timeline$Group <- as.numeric(timeline$Patient);

plot.data$Patient <- factor(plot.data$Patient, levels = rev(patient.order));
plot.data <- plot.data[order(plot.data$Patient),];
plot.data$Group <- as.numeric(plot.data$Patient);

# set colours
covariate.colours <- list(
	Age = list(
		colours = c('navajowhite','darkgoldenrod1','darkorange','darkorange3'),
		levels = c('40 - 50','50 - 60','60 - 70','70 - 80'),
		labels = c('40 - 50','50 - 60','60 - 70','70 - 80')
		),
	BRCA = list(
		colours = c('black','grey30','white'),
		levels = c('BRCA1','BRCA2','WT'),
		labels = c('BRCA1','BRCA2','WT')
		),
	Tissue = list(
		colours = c('#008B8B','#5BD2DC'),
		levels = c('HGSOC','PSEC'),
		labels = c('HGSOC','PSEC')
		),
	Response = list(
		colours = c('#843A1C','#5BCB8E','#D6DEFF','grey80'),
		levels = c('PD','PR','SD','IE'),
		labels = c('Progressive Disease','Partial Response','Stable Disease','Inevaluable')
		),
	Cohort = list(
		colours = c('#40CC61','#6093F7','#D99FDC'),
		levels = c('1','2','3'),
		labels = c('Platinum Sensitive','Platinum Resistant','Exploratory')
		),
	DFS = list(
		colours = c('#99C19A','#725F7A'),
		levels = c('0','1'),
		labels = c('Censored','Progressed')
		),
	OS = list(
		colours = c('#BCEE68','#AB82FF'),
		levels = c('0','1'),
		labels = c('Living','Deceased')
		)
	);		

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
	size = 0.5,
	grid.col = list(col = 'white', lwd = 1),
	col.lines = 1:length(smp.covariates)
	);

smp.legends <- list(
	legend = list(
		colours = covariate.colours$Cohort$colours,
		labels = covariate.colours$Cohort$labels,
		title = 'Cohort'
		),
	legend = list(
		colours = covariate.colours$Tissue$colours,
		labels = covariate.colours$Tissue$labels,
		title = 'Tissue'
		),
	legend = list(
		colours = covariate.colours$BRCA$colours,
		labels = covariate.colours$BRCA$labels,
		title = 'Germline BRCA Status'
		),
	legend = list(
		colours = covariate.colours$Age$colours,
		labels = covariate.colours$Age$labels,
		title = 'Age'
		),
	legend = list(
		colours = covariate.colours$Response$colours,
		labels = covariate.colours$Response$labels,
		title = 'Response'
		),
	legend = list(
		colours = covariate.colours$DFS$colours,
		labels = covariate.colours$DFS$labels,
		title = 'DFS'
		),
	legend = list(
		colours = covariate.colours$OS$colours,
		labels = covariate.colours$OS$labels,
		title = 'OS'
		)
	);

smp.legend.grob <- legend.grob(
	legends = smp.legends,
	label.cex = 0.8,
	title.cex = 0.8,
	title.fontface = 'plain',
	title.just = 'left',
	size = 1,
	layout = c(1, length(smp.legends))
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
	left.padding = 20,
	ylimits = c(0.5, max(timeline$Group)+0.5),
	yat = seq(1, max(timeline$Group),1),
	yaxis.lab = paste0(levels(timeline$Patient), '              '),
	yaxis.cex = 0.95,
	xaxis.cex = 1,
	xaxis.tck = c(0.5,0),
	yaxis.tck = c(0,0),
	xlab.label = 'Number of Cycles',
	xlab.cex = 1.5,
	ylab.label = NULL,
	legend = list(
		inside = list(fun = smp.covariate.grob, x = -0.295, y = 0.5),
		inside = list(fun = draw.key, args = list(dot.key), x = 0.63, y = 0.8),
		inside = list(fun = smp.legend.grob, x = -1.45, y = 0.96)
		),
	style = 'Nature'
	);

write.plot(
	timeline.plot,
	size.units = 'in',
	width = 6.5,
	height = 6.5,
	filename = generate.filename('EVOLVE_ctDNA', 'timeline_plot','png')
	);

# save data
save(
	all.patients,covariate.colours,smp.legend.grob,clinical,timeline,
	file = generate.filename('EVOLVE_ctDNA','clinicalCovariates','RData')
	);

# and for manuscript
setwd('/cluster/projects/ovgroup/projects/OV_Superset/EVOLVE/ctDNA/pipeline_suite/paper_figures');

#write.plot(
#	timeline.plot,
#	size.units = 'in',
#	width = 6.5,
#	height = 6.5,
#	filename = generate.filename('EVOLVE_ctDNA', '_figure1','png')
#	);

### SAVE SESSION INFO ##############################################################################
save.session.profile(generate.filename('TimelinePlot','SessionProfile','txt'));
