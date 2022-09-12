### survival_statistics.R ##########################################################################
# Perform survival statistics on clinical features

### PREPARE SESSION ################################################################################
# import libraries
library(BoutrosLab.plotting.general);
library(survival);

source('/cluster/home/sprokope/git/analysis/helper_functions/session.functions.R');
source('/cluster/home/sprokope/git/analysis/helper_functions/survival.functions.R');

working.dir <- '/cluster/projects/ovgroup/projects/OV_Superset/EVOLVE/ctDNA/pipeline_suite/fragment_size';

setwd(working.dir);

### READ DATA ######################################################################################
# find data files
fragment.sizes <- read.delim('2022-09-07_EVOLVE_ctDNA__fragment_size_and_purity_estimates.tsv');
phenodata <- read.delim('../configs/EVOLVE_clinical_data.tsv');
phenodata <- unique(phenodata[,c(1,8,11,16,17,20,21)]);
colnames(phenodata) <- c('Patient','Age','Best.Response','DFS.Months','DFS.Status','OS.Months','OS.Status');

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

# format fragment info	
fragment.sizes <- fragment.sizes[grepl('Screening|C1D1', fragment.sizes$Sample),];
fragment.sizes$FC <- log2(fragment.sizes$Median.ALT) / log2(fragment.sizes$Median.REF);

# merge mutation and clinical info
master.matrix <- merge(
	phenodata,
	fragment.sizes[,c('Patient','Prop.short','Prop.long','p.value','FC')],
	by = 'Patient'
	);

### SURVIVAL #######################################################################################
# loop over each variable
for (metric in c('Prop.short','Prop.long','p.value','FC')) {

	# organize groups
	surv.data <- master.matrix[!is.na(master.matrix[,metric]),];

	groups <- rep(0, nrow(surv.data));

	if (grepl('Prop', metric)) {
		groups[which(surv.data[,metric] > median(surv.data[,metric]))] <- 1;
		labels <- c('low','high');
		} else if (metric == 'p.value') {
		groups[which(surv.data[,metric] < 0.05)] <- 1;
		labels <- c('ns','signif');
		} else {
		groups[which(surv.data[,metric] < 1)] <- 1;
		labels <- c('+ve','-ve');
		}

	# collect survival stats
	survtime <- surv.data$OS.Months;
	survstat <- as.numeric(surv.data$OS.Status);
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
		risk.labels = labels,
		key.groups.labels = labels,
		key.groups.cex = 1.5,
		key.stats.corner = c(1,1),
		key.stats.x.pos = 1,
		key.stats.y.pos = 1,
		key.groups.corner = c(1,1),
		key.groups.x.pos = 1,
		key.groups.y.pos = 0.85,
		filename = generate.filename('EVOLVE_ctDNA__os_by', metric, 'png'),
		style = 'Nature'
		);

	# collect survival stats
	survtime <- surv.data$DFS.Months;
	survstat <- as.numeric(surv.data$DFS.Status);
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
		risk.labels = labels,
		key.groups.labels = labels,
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
		filename = generate.filename('EVOLVE_ctDNA_dfs_by',metric, 'png'),
		style = 'Nature'
		);
	}


# organize groups
surv.data <- master.matrix[!is.na(master.matrix$p.value),];

groups <- rep(0, nrow(surv.data));
groups[which(surv.data$p.value < 0.05)] <- 1;

plot.groups <- rep(0, nrow(surv.data));
plot.groups[which(surv.data$p.value < 0.05)] <- 1;
plot.groups[which(surv.data$p.value < 0.05 & surv.data$FC < 1)] <- 2;
labels <- c('ns','+ve','-ve');

# collect survival stats
survtime <- surv.data$DFS.Months;
survstat <- as.numeric(surv.data$DFS.Status);
survobj <- Surv(survtime, survstat);

output <- fit.coxmodel(
	groups = groups,
	survobj = survobj,
#	other.data = tmp.data[,key.clinical], # if length(key.clinical) > 1
#	other.data = matrix(master.matrix$Age),
	return.cox.model = TRUE
	);

# make a km plot
create.km.plot(
	survival.object = survobj,
	patient.groups = plot.groups,
	xlab.label = 'Time (months)',
	ylab.label = 'DFS Proportion',
	ylab.axis.padding = 3.5,
	ylab.cex = 2,
	xlab.cex = 2,
	yaxis.cex = 1.5,
	xaxis.cex = 1.5,
	line.colours = c('turquoise3','#CB74F5','violetred3'),
	risk.labels = labels,
	key.groups.labels = c('(mt/wt) = 0','(mt/wt) > 0', '(mt/wt) < 0'),
	key.groups.cex = 1.5,
	key.stats.corner = c(1,1),
	key.stats.x.pos = 1,
	key.stats.y.pos = 1,
	key.groups.corner = c(1,1),
	key.groups.x.pos = 1,
	key.groups.y.pos = 0.85,
	statistical.method = 'cox',
	explicit.HR.label = FALSE,
	predefined.p = 2*pnorm(-abs(summary(output)$coef[1,4])),
	predefined.hr = summary(output)$conf.int[1,1],
	predefined.hr.ci = summary(output)$conf.int[1,3:4],
	filename = generate.filename('EVOLVE_ctDNA_dfs_by', paste0(metric,'_alt'), 'png'),
	style = 'Nature'
	);

### SAVE SESSION INFO ##############################################################################
save.session.profile(generate.filename('SurvivalAnalysis','SessionProfile','txt'));
