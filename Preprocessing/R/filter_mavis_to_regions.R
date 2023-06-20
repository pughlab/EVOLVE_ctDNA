### combined_and_filter_mavis.R ####################################################################
# Combine calls from allUNIQUE, DCS, SSCS and DCS_SSCS bams and filter to target regions.

### PREAMBLE #######################################################################################
library(GenomicRanges);

### READ DATA ######################################################################################
# get mavis output
sv.data <- read.delim('2021-08-19_EVOLVE_ctDNA_mavis_output.tsv');

# get target regions
target_bed <- read.delim('../../CHARM-MMR_plus_EVOLVE_hg38.bed', header = F)
colnames(target_bed) <- c('Chromosome','Start','End','ID','V5','Strand');

# remove extra regions
target_genes <- subset(target_bed,!grepl("MSI|Agena",ID));

### FORMAT DATA ####################################################################################
# create genomic ranges object for target regions
target.gr <- makeGRangesFromDataFrame(target_genes, starts.in.df.are.0based = TRUE);

# create genomic ranges object for each breakpoint
first_bp <- data.frame(
	Chromosome = paste0('chr',sv.data$break1_chromosome),
	Start_Position = sv.data$break1_position_start,
	End_Position = sv.data$break1_position_end
	);

second_bp <- data.frame(
	Chromosome = paste0('chr',sv.data$break2_chromosome),
	Start_Position = sv.data$break2_position_start,
	End_Position = sv.data$break2_position_end
	);

bp1.gr <- makeGRangesFromDataFrame(
	first_bp,
	seqnames.field = 'Chromosome',
	start.field = 'Start_Position',
	end.field = 'End_Position'
	);

bp2.gr <- makeGRangesFromDataFrame(
	second_bp,
	seqnames.field = 'Chromosome',
	start.field = 'Start_Position',
	end.field = 'End_Position'
	);

# find overlaps
overlaps.p1 <- as.data.frame(findOverlaps(bp1.gr, target.gr));
overlaps.p2 <- as.data.frame(findOverlaps(bp2.gr, target.gr))

overlap.data <- merge(
	overlaps.p1,
	overlaps.p2,
	by = 'queryHits',
	suffixes = c('.1','.2'),
	all = TRUE
	);

to.remove <- which(is.na(overlap.data$subjectHits.1) | is.na(overlap.data$subjectHits.2));

# indicate indices to keep
keep.idx <- unique(overlap.data[-to.remove,]$queryHits);

# filter initial input
sv.data.filtered <- sv.data[keep.idx,];

### SAVE DATA ######################################################################################
# save filtered data
write.table(
	sv.data.filtered,
	file = generate.filename(arguments$project, 'mavis_output_filtered','tsv'),
	row.names = FALSE,
	col.names = TRUE,
	sep = '\t'
	);


