### figure5__plot_TP53_fragment_sizes.R ############################################################
# Examine differences in fragment sizes in cfDNA at TP53 mutation sites to differentiate between
# true germline/somatic/CHIP mutations.

### PREPARE SESSION ################################################################################
# import libraries
library(BoutrosLab.plotting.general);
library(ggplot2);

source(paste0(getwd(),'/helper_functions/session.functions.R'));

# get clinical covariates
load(paste0(getwd(),'/../data/EVOLVE_ctDNA__clinical_timeline.RData'));

# get fragment data from EVOLVE ctDNA samples
load(paste0(getwd(),'/../data/EVOLVE_ctDNA__TP53_fragment_sizes.RData'));

### MAIN ###########################################################################################	
plot.objects <- list();

# summarize all germline variant sizes from healthy controls (all germline variants in CHARM HBC dataset)
ref <- data.frame(V1 = germline.variants$ref);
alt <- data.frame(V1 = germline.variants$alt);
outfile <- 'hbc__all_germline_mutations.pdf';
mutation <- 'healthy (Germline)';

# is the CDF of ALT above (shorter) than ref?
p <- ks.test(alt$V1, ref$V1, alternative = 'g')$p.value
p <- display.statistical.result(p, statistic.type = 'p', symbol = ' = ');

wt <- length(ref$V1);
mut <- length(alt$V1);
median1 <- median(ref$V1, na.rm = TRUE);
median2 <- median(alt$V1, na.rm = TRUE);

# Plot cumulative distributions
plot.objects[[1]] <- ggplot() + 
	stat_ecdf(data = ref, aes(V1, color = "black"), geom = "step") +
	stat_ecdf(data = alt, aes(V1, color = "red"), geom = "step") +
	xlab("Fragment Length") + 
	ylab("Cumulative Distribution") +
	ggtitle(mutation) + 
	scale_colour_manual(name = 'Fragments', values =c("black" = "black", "red" = "red"), labels = c("Wildtype", "ALT")) +
	theme(plot.title = element_text(hjust = 0.5, face = "bold.italic", size = 15), 
		axis.line = element_line(colour = "black"),
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(),
		panel.border = element_rect(colour = "black", fill=NA, size=0.5),
		panel.background = element_blank(),
		legend.position = "none",
		legend.key = element_rect(fill = "white"),
		axis.text = element_text(size = 10),
		axis.title = element_text(size = 12),
		axis.text.x = element_text(angle = 0, hjust = 1)) + 
	scale_x_continuous(limits=c(75, 250, expand = c(0,0))) +
	geom_vline(xintercept = median1, linetype = "dashed", color = "black", size=0.5) +
	geom_vline(xintercept = median2, linetype = "dashed", color = "red", size=0.5) +
	annotate("text", x = 210, y = c(0.40, 0.27, 0.14), 
		label = c(p, paste0("WT = ", wt), paste0("Mutant = ", mut)),
		col = c("black", "black", "red"),
		size = 4);

# summarize known/novel/chip variant sizes
for (i in c('known','chip','novel')) {

	if (i == 'known') {
		ref <- data.frame(V1 = known.tp53.variants$ref);
		alt <- data.frame(V1 = known.tp53.variants$alt);
		outfile <- 'known_somatic_tp53_mutations.pdf';
		title <- 'TP53 (Somatic)';
		} else if (i == 'chip') {
		ref <- data.frame(V1 = suspected.chip.variants$ref);
		alt <- data.frame(V1 = suspected.chip.variants$alt);
		outfile <- 'suspected_CHIP_tp53_mutations.pdf';
		title <- 'TP53 (Suspected CHIP)';
		} else {
		ref <- data.frame(V1 = novel.tp53.variants$ref);
		alt <- data.frame(V1 = novel.tp53.variants$alt);
		outfile <- 'suspected_somatic_tp53_mutations.pdf';
		title <- 'TP53 (Suspected Somatic)';
		}

	wt <- length(ref$V1);
	mut <- length(alt$V1);
	median1 <- median(ref$V1, na.rm = TRUE);
	median2 <- median(alt$V1, na.rm = TRUE);
	
	p <- ks.test(alt$V1, ref$V1, alternative = 'g')$p.value;
	p <- display.statistical.result(p, statistic.type = 'p', symbol = ' = ', lower.cutoff = 1e-8, digits = 1);

	# Plot cumulative distributions
	plot <- ggplot() + 
		stat_ecdf(data = ref, aes(V1, color = "black"), geom = "step") +
		stat_ecdf(data = alt, aes(V1, color = "red"), geom = "step") +
		xlab("Fragment Length") + 
		ylab("Cumulative Distribution") +
		ggtitle(title) + 
		scale_colour_manual(name = 'Fragments', values =c("black" = "black", "red" = "red"), labels = c("Wildtype", "ALT")) +
		theme(plot.title = element_text(hjust = 0.5, face = "bold.italic", size = 15), 
			axis.line = element_line(colour = "black"),
			panel.grid.major = element_blank(),
			panel.grid.minor = element_blank(),
			panel.border = element_rect(colour = "black", fill=NA, size=0.5),
			panel.background = element_blank(),
			legend.position = "none",
			legend.key = element_rect(fill = "white"),
			axis.text = element_text(size = 10),
			axis.title = element_text(size = 12),
			axis.text.x = element_text(angle = 0, hjust = 1)) + 
		scale_x_continuous(limits=c(75, 250, expand = c(0,0))) +
		geom_vline(xintercept = median1, linetype = "dashed", color = "black", size=0.5) +
		geom_vline(xintercept = median2, linetype = "dashed", color = "red", size=0.5) +
		annotate("text", x = 215, y = c(0.40, 0.27, 0.14), 
			label = c(p, paste0("WT = ", wt), paste0("Mutant = ", mut)),
			col = c("black", "black", "red"),
			size = 4);
	}

# combine them!
create.multipanelplot(
	plot.objects = plot.objects,
	plot.objects.heights = c(1,1),
	plot.objects.widths = c(1,1),
	layout.height = 2,
	layout.width = 2,
	x.spacing = 1,
	y.spacing = 0,
	left.legend.padding = 0,
	right.legend.padding = 0,
	top.legend.padding = 0,
	bottom.legend.padding = 0,
	height = 5,
	width = 7,
	resolution = 200,
	filename = generate.filename('EVOLVE_ctDNA', 'TP53_classification__Figure5','png')
	);

### SAVE SESSION INFO ##############################################################################
save.session.profile(generate.filename('Figure5','SessionProfile','txt'));
