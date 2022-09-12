library(GenomicRanges)

charm.bed <- '/cluster/projects/ovgroup/projects/OV_Superset/EVOLVE/ctDNA/pipeline_suite/configs/intervals/CHARM-MMR_plus_EVOLVE_hg38_liftover_sorted_padding100bp.bed';

charm_bed <- read.delim(charm.bed, header = FALSE, comment.char = '#');
colnames(charm_bed)[1:3] <- c('Chromosome','Start','End');
charm.gr <- makeGRangesFromDataFrame(charm_bed, starts.in.df.are.0based = TRUE);

exome.bed <- '/cluster/projects/ovgroup/projects/OV_Superset/Combined/intervals/S04380219_Covered_3col_headless__hg38_trimmed_sorted_padding100bp.bed';

exome_bed <- read.delim(exome.bed, header = FALSE, comment.char = '#');
colnames(exome_bed)[1:3] <- c('Chromosome','Start','End');
exome.gr <- makeGRangesFromDataFrame(exome_bed, starts.in.df.are.0based = TRUE);

target.regions <- intersect(charm.gr, exome.gr);


### SV #############################################################################################
exome <- read.delim('2021-11-16_EVOLVE_wxs_sv_data_annotated_filtered.tsv', row.names = NULL);
ctdna <- read.delim('2021-11-18_EVOLVE_ctDNA_sv_data_for_cbioportal_annotated_filtered.tsv');
#ctdna <- rbind(
#	read.delim('2021-11-16_EVOLVE_ctDNA_sv_data_allUNIQUE_annotated_filtered.tsv'),
#	read.delim('2021-11-16_EVOLVE_ctDNA_sv_data_DCS_SSCS_annotated_filtered.tsv'),
#	read.delim('2021-11-16_EVOLVE_ctDNA_sv_data_SSCS_annotated_filtered.tsv')
#	);

exome$Patient <- substr(exome$Sample_ID, 0, 11);
ctdna$Patient <- substr(ctdna$Sample_ID, 0, 11);

# create genomic ranges object for each breakpoint
first_bp <- data.frame(
	Chromosome = paste0('chr',exome$Site1_Chromosome),
	Start = exome$Site1_Position,
	End = exome$Site1_Position
	);

second_bp <- data.frame(
	Chromosome = paste0('chr',exome$Site2_Chromosome),
	Start = exome$Site2_Position,
	End = exome$Site2_Position
	);

bp1.gr <- makeGRangesFromDataFrame(first_bp, starts.in.df.are.0based = FALSE);
bp2.gr <- makeGRangesFromDataFrame(second_bp, starts.in.df.are.0based = FALSE);

# find overlaps
overlaps.p1 <- as.data.frame(findOverlaps(bp1.gr, target.regions));
overlaps.p2 <- as.data.frame(findOverlaps(bp2.gr, target.regions));

overlap.data <- merge(
	overlaps.p1,
	overlaps.p2,
	by = 'queryHits',
	suffixes = c('.1','.2'),
	all = TRUE
	);

# only keep entries for which both breakpoints are within target regions
to.remove <- which(is.na(overlap.data$subjectHits.1) | is.na(overlap.data$subjectHits.2));
keep.exome.idx <- unique(overlap.data[-to.remove,]$queryHits);

exome.panel <- exome[keep.exome.idx,c('Patient','Sample_ID','Fusion','Class','Connection_Type','Site2_Effect_On_Frame','Site1_Exon','Site2_Exon')];

# create genomic ranges object for each breakpoint
first_bp <- data.frame(
	Chromosome = paste0('chr',ctdna$Site1_Chromosome),
	Start = ctdna$Site1_Position,
	End = ctdna$Site1_Position
	);

second_bp <- data.frame(
	Chromosome = paste0('chr',ctdna$Site2_Chromosome),
	Start = ctdna$Site2_Position,
	End = ctdna$Site2_Position
	);

bp1.gr <- makeGRangesFromDataFrame(first_bp, starts.in.df.are.0based = FALSE);
bp2.gr <- makeGRangesFromDataFrame(second_bp, starts.in.df.are.0based = FALSE);

# find overlaps
overlaps.p1 <- as.data.frame(findOverlaps(bp1.gr, target.regions));
overlaps.p2 <- as.data.frame(findOverlaps(bp2.gr, target.regions));

overlap.data <- merge(
	overlaps.p1,
	overlaps.p2,
	by = 'queryHits',
	suffixes = c('.1','.2'),
	all = TRUE
	);

# only keep entries for which both breakpoints are within target regions
to.remove <- which(is.na(overlap.data$subjectHits.1) | is.na(overlap.data$subjectHits.2));
keep.ctdna.idx <- unique(overlap.data[-to.remove,]$queryHits);

ctdna.panel <- ctdna[keep.ctdna.idx,c('Patient','Sample_ID','Fusion','Class','Connection_Type','Site2_Effect_On_Frame','Site1_Exon','Site2_Exon')];

sv.overlap.template <- data.frame(
	Patient = sort(as.character(unique(ctdna.panel$Patient))),
	Group = NA,
	N.exome = 0,
	N.ctdna = 0,
	Overlap = 0
	);

ctdna.subset <- ctdna.panel;
results.data <- list();

for (group in c('allUNIQUE','DCS-SSCS','SSCS')) {

	sv.overlap <- sv.overlap.template;
	sv.overlap$Group <- group;
	panel.data <- ctdna.subset[grepl(group, ctdna.subset$Sample_ID),];

	for (i in 1:nrow(sv.overlap)) {

		patient <- as.character(sv.overlap[i,]$Patient);
		tmp1 <- unique(exome.panel[which(exome.panel$Patient == patient),c('Fusion','Class','Connection_Type')]);
		tmp2 <- unique(panel.data[which(panel.data$Patient == patient),c('Fusion','Class','Connection_Type')]);

		sv.overlap[i,3] <- nrow(tmp1);
		sv.overlap[i,4] <- nrow(tmp2);

		if (nrow(tmp1) == 0) { next; }
		if (nrow(tmp2) == 0) { next; }

		tmp1$Call <- 1;
		tmp2$Call <- 1;

		tmp3 <- merge(
			tmp1,
			tmp2,
			by = c('Fusion','Class','Connection_Type'),
			all = TRUE,
			suffixes = c('.wxs','.panel')
			);

		sv.overlap[i,5] <- nrow(tmp3[which(tmp3$Call.wxs == 1 & tmp3$Call.panel == 1),]);

		}

	sv.overlap$TP <- sv.overlap$Overlap;
	sv.overlap$FP <- sv.overlap$N.ctdna - sv.overlap$Overlap;
	sv.overlap$FN <- sv.overlap$N.exome - sv.overlap$Overlap;

	sv.overlap$F1 <- apply(
		sv.overlap[,c('TP','FP','FN')],1,
		function(i) { 2 * i[1] / (2 * i[1] + i[2] + i[3] ) }
		);

	results.data[[group]] <- sv.overlap;
	ctdna.subset <- ctdna.subset[!grepl(group, ctdna.subset$Tumor_Sample_Barcode),];
	}

results.output <- do.call(rbind, results.data);

write.table(
	results.output,
	file = 'sv_overlap_statistics.tsv',
	row.names = FALSE,
	col.names = TRUE,
	sep = '\t'
	);

### SOMATIC ########################################################################################
exome <- read.delim('2021-04-29_EVOLVE_wxs_somatic_oncoKB_annotated_filtered.tsv');
ctdna <- read.delim('2021-08-11_EVOLVE_ctDNA_somatic_oncokb_annotated_filtered.tsv');

exome$Patient <- substr(exome$Tumor_Sample_Barcode, 0, 11);
ctdna$Patient <- substr(ctdna$Tumor_Sample_Barcode, 0, 11);

exome$t_vaf <- exome$t_alt_count / exome$t_depth;
ctdna$t_vaf <- ctdna$t_alt_count / ctdna$t_depth;

exome <- exome[,c('Patient','Tumor_Sample_Barcode','Chromosome','Start_Position','End_Position','Allele','t_vaf')];
ctdna <- ctdna[,c('Patient','Tumor_Sample_Barcode','Chromosome','Start_Position','End_Position','Allele','t_vaf')];

exome.gr <- makeGRangesFromDataFrame(exome, start.field = 'Start_Position', end.field = 'End_Position');
target.exome.overlaps <- as.data.frame(findOverlaps(exome.gr, target.regions));
exome.panel <- exome[unique(target.exome.overlaps$queryHits),];

ctdna.gr <- makeGRangesFromDataFrame(ctdna, start.field = 'Start_Position', end.field = 'End_Position');
target.ctdna.overlaps <- as.data.frame(findOverlaps(ctdna.gr, target.regions));
ctdna.panel <- ctdna[unique(target.ctdna.overlaps$queryHits),];

somatic.overlap.template <- data.frame(
	Patient = unique(ctdna.panel$Patient),
	Group = NA,
	N.exome = 0,
	N.ctdna = 0,
	Overlap = 0
	);

ctdna.subset <- ctdna.panel;
results.data <- list();

for (group in c('allUNIQUE','DCS_SSCS','DCS','SSCS')) {

	somatic.overlap <- somatic.overlap.template;
	somatic.overlap$Group <- group;
	panel.data <- ctdna.subset[grepl(group, ctdna.subset$Tumor_Sample_Barcode),];

	for (i in 1:nrow(somatic.overlap)) {

		patient <- as.character(somatic.overlap[i,]$Patient);
		tmp1 <- unique(exome.panel[which(exome.panel$Patient == patient),c('Chromosome','Start_Position','End_Position','Allele')]);
		tmp2 <- unique(panel.data[which(panel.data$Patient == patient),c('Chromosome','Start_Position','End_Position','Allele')]);

		somatic.overlap[i,3] <- nrow(tmp1);
		somatic.overlap[i,4] <- nrow(tmp2);

		if (nrow(tmp1) == 0) { next; }
		if (nrow(tmp2) == 0) { next; }

		tmp1$Call <- 1;
		tmp2$Call <- 1;

		tmp3 <- merge(
			tmp1,
			tmp2,
			by = c('Chromosome','Start_Position','End_Position','Allele'),
			all = TRUE,
			suffixes = c('.wxs','.panel')
			);

		somatic.overlap[i,5] <- nrow(tmp3[which(tmp3$Call.wxs == 1 & tmp3$Call.panel == 1),]);

		}

	somatic.overlap$TP <- somatic.overlap$Overlap;
	somatic.overlap$FP <- somatic.overlap$N.ctdna - somatic.overlap$Overlap;
	somatic.overlap$FN <- somatic.overlap$N.exome - somatic.overlap$Overlap;

	somatic.overlap$F1 <- apply(
		somatic.overlap[,c('TP','FP','FN')],1,
		function(i) { 2 * i[1] / (2 * i[1] + i[2] + i[3] ) }
		);

	results.data[[group]] <- somatic.overlap;
	ctdna.subset <- ctdna.subset[!grepl(group, ctdna.subset$Tumor_Sample_Barcode),];
	}

results.output <- do.call(rbind, results.data);

write.table(
	results.output,
	file = 'somatic_overlap_statistics.tsv',
	row.names = FALSE,
	col.names = TRUE,
	sep = '\t'
	);

### GERMLINE #######################################################################################
exome <- read.delim('2021-04-22_EVOLVE_wxs_germline_oncokb_annotated.tsv')
ctdna <- read.delim('2021-06-04_EVOLVE_ctDNA_germline_oncokb_annotated.tsv')

exome$Patient <- substr(exome$Tumor_Sample_Barcode, 0, 11);
ctdna$Patient <- substr(ctdna$Tumor_Sample_Barcode, 0, 11);

exome <- exome[,c('Patient','Tumor_Sample_Barcode','Chromosome','Start_Position','End_Position','Allele')];
ctdna <- ctdna[,c('Patient','Tumor_Sample_Barcode','Chromosome','Start_Position','End_Position','Allele')];

exome.gr <- makeGRangesFromDataFrame(exome, start.field = 'Start_Position', end.field = 'End_Position');
target.exome.overlaps <- as.data.frame(findOverlaps(exome.gr, target.regions));
exome.panel <- exome[unique(target.exome.overlaps$queryHits),];

ctdna.gr <- makeGRangesFromDataFrame(ctdna, start.field = 'Start_Position', end.field = 'End_Position');
target.ctdna.overlaps <- as.data.frame(findOverlaps(ctdna.gr, target.regions));
ctdna.panel <- ctdna[unique(target.ctdna.overlaps$queryHits),];

germline.overlap.template <- data.frame(
	Patient = unique(ctdna.panel$Patient),
	Group = NA,
	N.exome = 0,
	N.ctdna = 0,
	Overlap = 0
	);

ctdna.subset <- ctdna.panel;
results.data <- list();

for (group in c('allUNIQUE','DCS_SSCS','DCS','SSCS')) {

	germline.overlap <- germline.overlap.template;
	germline.overlap$Group <- group;
	panel.data <- ctdna.subset[grepl(group, ctdna.subset$Tumor_Sample_Barcode),];

	for (i in 1:nrow(germline.overlap)) {

		patient <- as.character(germline.overlap[i,]$Patient);
		tmp1 <- unique(exome.panel[which(exome.panel$Patient == patient),c('Chromosome','Start_Position','End_Position','Allele')]);
		tmp2 <- unique(panel.data[which(panel.data$Patient == patient),c('Chromosome','Start_Position','End_Position','Allele')]);

		germline.overlap[i,3] <- nrow(tmp1);
		germline.overlap[i,4] <- nrow(tmp2);

		if (nrow(tmp1) == 0) { next; }
		if (nrow(tmp2) == 0) { next; }

		tmp1$Call <- 1;
		tmp2$Call <- 1;

		tmp3 <- merge(
			tmp1,
			tmp2,
			by = c('Chromosome','Start_Position','End_Position','Allele'),
			all = TRUE,
			suffixes = c('.wxs','.panel')
			);

		germline.overlap[i,5] <- nrow(tmp3[which(tmp3$Call.wxs == 1 & tmp3$Call.panel == 1),]);

		}

	germline.overlap$TP <- germline.overlap$Overlap;
	germline.overlap$FP <- germline.overlap$N.ctdna - germline.overlap$Overlap;
	germline.overlap$FN <- germline.overlap$N.exome - germline.overlap$Overlap;

	germline.overlap$F1 <- apply(
		germline.overlap[,c('TP','FP','FN')],1,
		function(i) { 2 * i[1] / (2 * i[1] + i[2] + i[3] ) }
		);

	results.data[[group]] <- germline.overlap;
	ctdna.subset <- ctdna.subset[!grepl(group, ctdna.subset$Tumor_Sample_Barcode),];
	}

results.output <- do.call(rbind, results.data);

write.table(
	results.output,
	file = 'germline_overlap_statistics.tsv',
	row.names = FALSE,
	col.names = TRUE,
	sep = '\t'
	);
