### extract_fragment_sizes.R #######################################################################

### PREPARE SESSION ################################################################################
# import libraries
library(argparse);

# import command line arguments
parser <- ArgumentParser();

parser$add_argument('-i', '--input', type = 'character', help = 'path to input file (BAM)');
parser$add_argument('-o', '--output', type = 'character', help = 'path to output directory');
parser$add_argument('-s', '--sample', type = 'character', help = 'project name');
parser$add_argument('-r', '--ref_type', type = 'character', help = 'reference genome (hg38 or hg19)',
	default = 'hg38');
parser$add_argument('-c', '--chr', type = 'logical', help = 'should "chr" be appended to chromosome names',
	default = TRUE);
parser$add_argument('-t', '--targets', type = 'character', help = 'path to target regions', default = NULL);

arguments <- parser$parse_args();

library(BSgenome);
library(biovizBase);

# Rscript /cluster/home/sprokope/git/analysis/EVOLVE_ctDNA/extract_fragment_sizes.R -i /cluster/projects/ovgroup/projects/OV_Superset/EVOLVE/ctDNA/ConsensusCruncher/output/consensus/TGL59_0024_Pl_T_PE_358_TS_211110_NB551056_0232_AHMGCHBGXK_1_TTGCGAAG-AACACGCT.sorted/dcs_sc/TGL59_0024_Pl_T_PE_358_TS_211110_NB551056_0232_AHMGCHBGXK_1_TTGCGAAG-AACACGCT.sorted.all.unique.dcs.sorted.bam -o /cluster/projects/ovgroup/projects/OV_Superset/EVOLVE/ctDNA/pipeline_suite/fragment_size/test/EVO-009-024/EVO-009-024_ctDNA_C2D1_allUNIQUE_per_bin_sizes.tsv -s EVO-009-024_ctDNA_C2D1_allUNIQUE -r hg38 -t /cluster/projects/ovgroup/projects/OV_Superset/EVOLVE/ctDNA/pipeline_suite/configs/intervals/CHARM-MMR_plus_EVOLVE_hg38_liftover_sorted.bed

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

# function to expand cigar string
expand.cigar <- function(cigar) {

	string <- cigar;
	split.list <- unlist(strsplit(string, paste(LETTERS, collapse = '|')));
	output <- '';

	for (i in split.list) {

		string <- sub(paste0('^',i),'',string);
		first.char <- substr(string,0,1);
		string <- sub(paste0('^',first.char),'',string);

		output <- paste0(output, paste(rep(first.char, as.numeric(i)), collapse = ''));

        	}

	return(output);
	}

# function to read in BAM and extract fragment information
createFragSet <- function (sampleID, filePath, refGenome = "BSgenome.Hsapiens.UCSC.hg38", chr.select = paste0("chr", 1:22), mapQthreshold = 30, targets = NULL) {

	indexed.bam <- gsub("$", ".bai", filePath)
	if (!file.exists(indexed.bam)) { indexBam(filePath); }

	param <- Rsamtools::ScanBamParam(
		mapqFilter = mapQthreshold,
		which = makeGRangesFromDataFrame(targets, keep.extra.columns = TRUE),
		flag = Rsamtools::scanBamFlag(
			isPaired = TRUE,
			isProperPair = TRUE,
			isDuplicate = FALSE,
			isSecondaryAlignment = FALSE,
			isUnmappedQuery = FALSE
			)
		);

	galp <- GenomicAlignments::readGAlignmentPairs(filePath, param = param);
	frags <- granges(
		GenomeInfoDb::keepSeqlevels(galp, chr.select, pruning.mode = "coarse"),
		on.discordant.seqnames = "drop"
		);

	# determine if reads are REF or ALT
	sam <- data.frame(galp);

	alt.pos <- targets$Start[1];
	alt.size <- (targets$End[1] - targets$Start[1]) + 1;

	sam$Group <- 'REF';

	for (i in 1:nrow(sam)) {

		# is the target in the first or second read?
		read <- if (sam[i,]$end.first > alt.pos) { 'first' } else { 'last' };

		# first, see if there are any mismatches
		expected <- paste0(sam[i,paste0('width.', read)],'M');
		if (sam[i,paste0('cigar.',read)] == expected) { next; }

		# how far is it from the beginning of the read?
		pos1 <- alt.pos - sam[i,paste0('start.',read)] + 1;
		pos2 <- alt.pos - sam[i,paste0('start.',read)] + alt.size;

		# what should the cigar be?
		type <- unlist(strsplit(
			substr(expand.cigar(sam[i,paste0('cigar.',read)]),pos1,pos2),''));

		# do we see this?
		if (all(type == targets$Type)) { sam[i,]$Group <- 'ALT'; } 
		}

	frags$Group <- sam$Group;

	fs <- list(
		sampleID = sampleID,
		refGenome = refGenome,
		frags_Granges = frags,
		galp_len = length(galp),
		params = param
		);

	return(fs);
	}

### testing
#arguments$input <- '/cluster/projects/ovgroup/projects/OV_Superset/EVOLVE/ctDNA/ConsensusCruncher/output_batch1/consensus/Screening_indel_realined_TGL59_0001_Pl_n_PE_416_TS_210423_NB551051_0207_AHVVHVBGXH_1_TGAGCTAG-GAACGGTT_R1.fastq.gz.sorted/dcs_sc/Screening_indel_realined_TGL59_0001_Pl_n_PE_416_TS_210423_NB551051_0207_AHVVHVBGXH_1_TGAGCTAG-GAACGGTT_R1.fastq.gz.sorted.all.unique.dcs.sorted.bam';
#arguments$sample <- 'EVO-009-001_ctDNA_Screening';
#arguments$ref_type <- 'hg38';
#arguments$chr <- TRUE;
#arguments$targets <- '/cluster/projects/ovgroup/projects/OV_Superset/EVOLVE/ctDNA/pipeline_suite/configs/intervals/CHARM-MMR_plus_EVOLVE_hg38_liftover_sorted.bed';
#arguments$output <- '/cluster/projects/ovgroup/projects/OV_Superset/EVOLVE/ctDNA/pipeline_suite/fragment_size/test/EVO-009-001/EVO-009-001_ctDNA_Screening_allUNIQUE_per_bin_sizes.tsv'
###

### MAIN ###########################################################################################
# initiate genome settings
if (arguments$chr) {
	chromosomes <- paste0('chr',1:22);
	} else {
	chromosomes <- 1:22;
	}
if (arguments$ref_type == 'hg38') {
	genome <- 'BSgenome.Hsapiens.UCSC.hg38';
	} else {
	genome <- 'BSgenome.Hsapiens.UCSC.hg19';
	}
# read in target positions (+ any annotations; ie, Sample or Gene names)
target.gr <- NULL;
if (!is.null(arguments$targets)) {
	intervals <- read.delim(arguments$targets, header = FALSE, stringsAsFactors = FALSE);
	colnames(intervals) <- c('Chromosome','Start','End','Name','Type');
	chromosomes <- intersect(chromosomes, unique(intervals$Chromosome));
	}

results <- intervals[,1:4];
results[,c('short','normal','long','ratio.short','ratio.long','nfrags','coverage')] <- NA;

# loop over each variant
for (indel in 1:nrow(intervals)) {

	# format fragment data
	fset <- createFragSet(
		sampleID = arguments$sample,
		filePath = arguments$input,
		refGenome = genome,
		chr.select = chromosomes,
		targets = intervals[indel,],
		mapQthreshold = 30
		);

	# overlap fragment data with target intervals
	frags.gr <- fset$frags_Granges;
	ref.fragments <- frags.gr[which(frags.gr$Group == 'REF')];
	alt.fragments <- frags.gr[which(frags.gr$Group == 'ALT')];

	write(width(ref.fragments), file = paste0(arguments$sample, '__', intervals$Name[indel],'_ref.reads'), ncolumns = 1);
	write(width(alt.fragments), file = paste0(arguments$sample, '__', intervals$Name[indel],'_alt.reads'), ncolumns = 1);

	# count fragment sizes per bin
	counts <- table(c(width(ref.fragments), width(alt.fragments)));

	short.idx <- as.character(90:150);
	norm.idx <- as.character(151:230);
	long.idx <- as.character(231:320);

	results[indel,]$short <- sum(counts[intersect(names(counts),short.idx)]);	# 90-150 bp
	results[indel,]$normal <- sum(counts[intersect(names(counts),norm.idx)]);	# 151-230 bp
	results[indel,]$long <- sum(counts[intersect(names(counts),long.idx)]);   # 231-320 bp

	ratio.short <- results[indel,]$short/results[indel,]$normal;
	ratio.short[is.nan(ratio.short) | is.infinite(ratio.short)] <- NA;
	results[indel,]$ratio.short <- ratio.short;

	ratio.long <- results[indel,]$long/results[indel,]$normal;
	ratio.long[is.nan(ratio.long) | is.infinite(ratio.long)] <- NA;
	results[indel,]$ratio.long <- ratio.long;

	results[indel,]$nfrags <- apply(results[indel,c('short', 'normal', 'long')],1,sum);
	results[indel,]$coverage <- results[indel,]$nfrags/sum(results[indel,]$nfrags, na.rm = TRUE);
	}

# format for output
write.table(
	results,
	file = paste0(arguments$sample, '__fragment_ratios.tsv'),
	row.names = FALSE,
	col.names = TRUE,
	sep = '\t'
	);

### SAVE SESSION INFO ##############################################################################
save.session.profile(generate.filename('ExtractFragmentSizes','SessionProfile','txt'));
