### survival_statistics.R ##########################################################################
# Perform survival statistics on clinical features

### PREPARE SESSION ################################################################################
# import libraries
library(BoutrosLab.plotting.general);
library(survival);

source('/cluster/home/sprokope/git/analysis/helper_functions/session.functions.R');
source('/cluster/home/sprokope/git/analysis/helper_functions/survival.functions.R');

working.dir <- '/cluster/projects/ovgroup/projects/OV_Superset/EVOLVE/ctDNA/pipeline_suite/Ensemble_calls/tumour_content';

setwd(working.dir);

### READ DATA ######################################################################################
# find data files
clearance.data <- read.delim('2022-09-06_EVOLVE_ctDNA_tumourClearance_data.tsv');
phenodata <- read.delim('../../configs/EVOLVE_clinical_data.tsv');
phenodata <- unique(phenodata[,c(1,8,11,13,16,17,20,21)]);
colnames(phenodata) <- c('Patient','Age','Best.Response','Cohort','DFS.Months','DFS.Status','OS.Months','OS.Status');

### FORMAT DATA ####################################################################################
# format clinical info
phenodata$Best.Response <- factor(
	phenodata$Best.Response,
	levels = c('Stable Disease','Partial Response','Progressive Disease (Objective)'),
	labels = c('SD','PR','PD')
	);

phenodata$DFS.Status <- factor(
	phenodata$DFS.Status,
	levels = c('DiseaseFree','Recurred/Progressed'),
	labels = c(0,1)
	);

phenodata$OS.Status <- factor(
	phenodata$OS.Status,
	levels = c('LIVING','DECEASED'),
	labels = c(0,1)
	);

# merge mutation and clinical info
master.matrix <- merge(
	phenodata,
	clearance.data[,c('Patient.ID','delta','N.cycles')],
	by.x = 'Patient',
	by.y = 'Patient.ID'
	);

### SURVIVAL #######################################################################################
# organize groups
surv.data <- master.matrix[!is.na(master.matrix$delta),];

groups <- rep(0, nrow(surv.data));
groups[which(surv.data$delta > -10)] <- 1;

# collect survival stats
survtime <- surv.data$OS.Months;
survstat <- as.numeric(surv.data$OS.Status);
survobj <- Surv(survtime, survstat);

output <- fit.coxmodel(
	groups = groups,
	survobj = survobj,
#	other.data = surv.data[,c('Cohort','Age')],
	return.cox.model = TRUE
	);

# check if logrank should be used instead
if ( ph.fails(output) | (nlevels(groups) > 2) ) {

	output <- logrank.analysis(
		groups = groups,
		survival.object = survobj
		);
	}

# make a km plot
create.km.plot(
	survival.object = survobj,
	patient.groups = groups,
	xlab.label = 'Time (months)',
	ylab.label = 'OS Proportion',
	ylab.axis.padding = 3.5,
	ylab.cex = 2,
	xlab.cex = 2,
	yaxis.cex = 1.5,
	xaxis.cex = 1.5,
	risk.labels = c('-ve','+ve'),
	key.groups.labels = c('-ve','+ve'),
	key.groups.cex = 1.5,
	key.stats.corner = c(1,1),
	key.stats.x.pos = 1,
	key.stats.y.pos = 1,
	key.groups.corner = c(1,1),
	key.groups.x.pos = 1,
	key.groups.y.pos = 0.85,
	filename = generate.filename('EVOLVE_ctDNA', 'os_by_clearance', 'png'),
	style = 'Nature'
	);

# collect survival stats
survtime <- surv.data$DFS.Months;
survstat <- as.numeric(surv.data$DFS.Status);
survobj <- Surv(survtime, survstat);

output <- fit.coxmodel(
	groups = groups,
	survobj = survobj,
#	other.data = surv.data[,c('Cohort','Age')],
	return.cox.model = TRUE
	);

# check if logrank should be used instead
if ( ph.fails(output) | (nlevels(groups) > 2) ) {

	output <- logrank.analysis(
		groups = groups,
		survival.object = survobj
		);
	}

# make a km plot
create.km.plot(
	survival.object = survobj,
	patient.groups = groups,
	xlab.label = 'Time (months)',
	ylab.label = 'DFS Proportion',
	ylab.axis.padding = 3.5,
	ylab.cex = 2,
	xlab.cex = 2,
	yaxis.cex = 1.5,
	xaxis.cex = 1.5,
	line.colours = c('turquoise3','violetred3'),
	risk.labels = c('-ve','+ve'),
	key.groups.labels = c('-ve','+ve'),
	key.groups.cex = 1.5,
	key.stats.corner = c(1,1),
	key.stats.x.pos = 1,
	key.stats.y.pos = 1,
	key.groups.corner = c(1,1),
	key.groups.x.pos = 1,
	key.groups.y.pos = 0.85,
	statistical.method = 'cox',
#	predefined.p = 2*pnorm(-abs(summary(output)$coef[1,4])),
#	statistical.result.hr = summary(output)$conf.int[1,1],
	filename = generate.filename('EVOLVE_ctDNA', 'dfs_by_clearance', 'png'),
	style = 'Nature'
	);

### N CYCLES #######################################################################################
# organize groups
groups <- rep(0, nrow(master.matrix));
groups[which(master.matrix$N.cycles > 3)] <- 1;

# collect survival stats
survtime <- master.matrix$OS.Months;
survstat <- as.numeric(master.matrix$OS.Status);
survobj <- Surv(survtime, survstat);

output <- fit.coxmodel(
	groups = groups,
	survobj = survobj,
#	other.data = tmp.data[,key.clinical], # if length(key.clinical) > 1
#	other.data = matrix(as.numeric(master.matrix$Age)),
	return.cox.model = TRUE
	);

# check if logrank should be used instead
if ( ph.fails(output) | (nlevels(groups) > 2) ) {

	output <- logrank.analysis(
		groups = groups,
		survival.object = survobj
		);
	}

# make a km plot
create.km.plot(
	survival.object = survobj,
	patient.groups = groups,
	xlab.label = 'Time (months)',
	ylab.label = 'OS Proportion',
	ylab.axis.padding = 3.5,
	ylab.cex = 2,
	xlab.cex = 2,
	yaxis.cex = 1.5,
	xaxis.cex = 1.5,
	risk.labels = c('low','high'),
	key.groups.labels = c(paste0('\u2264','3 cycles'),'>3 cycles'),
	key.groups.cex = 1.5,
	key.stats.corner = c(1,1),
	key.stats.x.pos = 1,
	key.stats.y.pos = 1,
	key.groups.corner = c(1,1),
	key.groups.x.pos = 1,
	key.groups.y.pos = 0.85,
	filename = generate.filename('EVOLVE_ctDNA', 'os_by_nCycles', 'png'),
	style = 'Nature'
	);

# collect survival stats
survtime <- master.matrix$DFS.Months;
survstat <- as.numeric(master.matrix$DFS.Status);
survobj <- Surv(survtime, survstat);

output <- fit.coxmodel(
	groups = groups,
	survobj = survobj,
#	other.data = tmp.data[,key.clinical], # if length(key.clinical) > 1
#	other.data = matrix(master.matrix$Age),
	return.cox.model = TRUE
	);

# check if logrank should be used instead
if ( ph.fails(output) | (nlevels(groups) > 2) ) {

	output <- logrank.analysis(
		groups = groups,
		survival.object = survobj
		);
	}

# make a km plot
create.km.plot(
	survival.object = survobj,
	patient.groups = groups,
	xlab.label = 'Time (months)',
	ylab.label = 'DFS Proportion',
	ylab.axis.padding = 3.5,
	ylab.cex = 2,
	xlab.cex = 2,
	yaxis.cex = 1.5,
	xaxis.cex = 1.5,
	line.colours = c('turquoise3','violetred3'),
	risk.labels = c('low','high'),
	key.groups.labels = c(paste0('\u2264','3 cycles'),'>3 cycles'),
	key.groups.cex = 1.5,
	key.stats.corner = c(1,1),
	key.stats.x.pos = 1,
	key.stats.y.pos = 1,
	key.groups.corner = c(1,1),
	key.groups.x.pos = 1,
	key.groups.y.pos = 0.85,
	statistical.method = 'logrank',
#	predefined.p = 2*pnorm(-abs(summary(output)$coef[1,4])),
#	statistical.result.hr = summary(output)$conf.int[1,1],
	filename = generate.filename('EVOLVE_ctDNA', 'dfs_by_nCycles', 'png'),
	style = 'Nature'
	);

### SAVE SESSION INFO ##############################################################################
save.session.profile(generate.filename('SurvivalAnalysis','SessionProfile','txt'));
