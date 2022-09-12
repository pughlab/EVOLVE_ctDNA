### survival_statistics.R ##########################################################################
# Perform survival statistics on clinical features

### PREPARE SESSION ################################################################################
# import libraries
library(BoutrosLab.plotting.general);
library(survival);

source('/cluster/home/sprokope/git/analysis/helper_functions/session.functions.R');
source('/cluster/home/sprokope/git/analysis/helper_functions/survival.functions.R');

### READ DATA ######################################################################################
# find data files
phenodata <- read.delim('EVOLVE_clinical_data.tsv');
phenodata <- unique(phenodata[,c(1,8,11,16,17,20,21)]);
colnames(phenodata) <- c('Patient','Age','Best.Response','DFS.Months','DFS.Status','OS.Months','OS.Status');

### FORMAT DATA ####################################################################################
# format clinical info
phenodata$Best.Response <- factor(
	phenodata$Best.Response,
	levels = c('Stable Disease','Partial Response','Progressive Disease (Objective)'),
	labels = c('SD','PR','PD')
	);

phenodata$Age <- factor(
	phenodata$Age,
	levels = c('< 40', '40 - 50','50 - 60','60 - 70','>= 70'), 
	labels = c('< 40', '40 - 50','50 - 60','60 - 70','>= 70')
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
master.matrix <- phenodata;
master.matrix$Reversion <- 0;
master.matrix[which(master.matrix$Patient.ID %in% c('EVO-009-001','EVO-009-003','EVO-009-013','EVO-009-023','EVO-009-006','EVO-009-007','EVO-009-018','EVO-009-021')),]$Reversion <- 1;

### SURVIVAL #######################################################################################
# organize groups
groups <- master.matrix$Reversion;

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
	risk.labels = c('0','1'),
	key.groups.labels = c('0','1'),
	key.groups.cex = 1.5,
	key.stats.corner = c(1,1),
	key.stats.x.pos = 1,
	key.stats.y.pos = 1,
	key.groups.corner = c(1,1),
	key.groups.x.pos = 1,
	key.groups.y.pos = 0.85,
	filename = generate.filename('EVOLVE_ctDNA', 'os_by_brca_reversion', 'png'),
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
	risk.labels = c('0','1'),
	key.groups.labels = c('0','1'),
	key.groups.cex = 1.5,
	key.stats.corner = c(1,1),
	key.stats.x.pos = 1,
	key.stats.y.pos = 1,
	key.groups.corner = c(1,1),
	key.groups.x.pos = 1,
	key.groups.y.pos = 0.85,
	filename = generate.filename('EVOLVE_ctDNA', 'dfs_by_brca_reversion', 'png'),
	style = 'Nature'
	);

### SAVE SESSION INFO ##############################################################################
save.session.profile(generate.filename('SurvivalAnalysis','SessionProfile','txt'));
