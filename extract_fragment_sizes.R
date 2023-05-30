### extract_fragment_sizes.R #######################################################################

### PREPARE SESSION ################################################################################
# import libraries
library(argparse);

# import command line arguments
parser <- ArgumentParser();

parser$add_argument('-i', '--input', type = 'character', help = 'path to input file (BAM)');
parser$add_argument('-o', '--output_dir', type = 'character', help = 'path to output directory');
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

# function to read in BAM and extract fragment information
createFragSet <- function (sampleID, filePath, refGenome = "BSgenome.Hsapiens.UCSC.hg38", chr.select = paste0("chr", 1:22), mapQthreshold = 30) {

	indexed.bam = gsub("$", ".bai", filePath)
	if (!file.exists(indexed.bam)) { indexBam(filePath); }

	param <- Rsamtools::ScanBamParam(
		mapqFilter = mapQthreshold,
		flag = Rsamtools::scanBamFlag(
			isPaired = TRUE,
			isProperPair = TRUE,
			isDuplicate = FALSE,
			isSecondaryAlignment = FALSE,
			isUnmappedQuery = FALSE
			)
		);

	galp <- GenomicAlignments::readGAlignmentPairs(filePath, param = param);
	galp_len <- length(galp);
	frags <- granges(
		GenomeInfoDb::keepSeqlevels(galp, chr.select, pruning.mode = "coarse"),
		on.discordant.seqnames = "drop"
		);

	w.all <- width(frags);
	hsapiens <- getBSgenome(refGenome);
	gcs <- GCcontent(hsapiens, unstrand(frags));
	frags$gc <- gcs;

	fs <- list(
		sampleID = sampleID,
		refGenome = refGenome,
		frags_Granges = frags,
		galp_len = galp_len,
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
	intervals <- read.delim(arguments$targets, header = FALSE, stringsAsFactors = FALSE)[,1:4];
	colnames(intervals) <- c('Chromosome','Start','End','Name');
	chromosomes <- intersect(chromosomes, unique(intervals$Chromosome));
	target.gr <- makeGRangesFromDataFrame(intervals, keep.extra.columns = TRUE);
	}

# format fragment data
fset <- createFragSet(
	sampleID = arguments$sample,
	filePath = arguments$input,
	refGenome = genome,
	chr.select = chromosomes,
	mapQthreshold = 30
	);

frags.gr <- fset$frags_Granges;

# move to output directory
setwd(arguments$output_dir);

# save fragment sizes
save(
	frags.gr,
	file = generate.filename(arguments$sample, 'fragment_sizes','RData')
	);

# get counts per interval
if (!is.null(arguments$targets)) {

	# overlap fragment data with target intervals
	fragments <- frags.gr[queryHits(findOverlaps(frags.gr, target.gr))];

	# count fragment sizes per bin
	w <- width(fragments);
	frag.list <- split(fragments, w);
	counts <- sapply(frag.list, function(x) countOverlaps(target.gr, x));

	short.idx <- as.character(90:150); # define short reads (90-150 bp)
	norm.idx <- as.character(151:230); # define 'normal' reads (151-230 bp)
	long.idx <- as.character(231:320); # define long reads (231-320 bp)

	intervals$short <- rowSums(counts[,intersect(colnames(counts),short.idx)]);	# 90-150 bp
	intervals$normal <- rowSums(counts[,intersect(colnames(counts),norm.idx)]);	# 151-230 bp
	intervals$long <- rowSums(counts[,intersect(colnames(counts),long.idx)]);       # 231-320 bp

	ratio.short <- intervals$short/intervals$normal;
	ratio.short[is.nan(ratio.short) | is.infinite(ratio.short)] <- NA;
	intervals$ratio.short <- ratio.short;

	ratio.long <- intervals$long/intervals$normal;
	ratio.long[is.nan(ratio.long) | is.infinite(ratio.long)] <- NA;
	intervals$ratio.long <- ratio.long;

	intervals$nfrags <- rowSums(counts); #apply(intervals[,c('short', 'normal', 'long')],1,sum);
	intervals$coverage <- intervals$nfrags/sum(intervals$nfrags, na.rm = TRUE);

	intervals$Prop.short <- intervals$short / intervals$nfrags;
	intervals$Prop.long <- intervals$long / intervals$nfrags;

	# write output to file
	write.table(
		as.data.frame(intervals),
		file = generate.filename(arguments$sample, 'per_bin_sizes','tsv'),
		row.names = FALSE,
		col.names = TRUE,
		sep = '\t'
		);
	}

### SAVE SESSION INFO ##############################################################################
save.session.profile(generate.filename('ExtractFragmentSizes','SessionProfile','txt'));
