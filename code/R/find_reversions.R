### find_reversions.R ##############################################################################
# Compare SNV/INDEL/SV data to find potential reversions.

### PREPARE SESSION ################################################################################
# import libraries
library(findReversions);
library(argparse);

# source in custom functions
source('~/git/analysis/helper_functions/session.functions.R');

# get command line args
parser <- ArgumentParser();

parser$add_argument('-d', '--directory', type = 'character', help = 'path to working directory',
	default = '/cluster/projects/ovgroup/projects/OV_Superset/EVOLVE/ctDNA/pipeline_suite/Reversions');
parser$add_argument('-p', '--project', type = 'character', help = 'project name');
parser$add_argument('-q', '--query', type = 'character', help = 'query maf (or tsv) file');
parser$add_argument('-m', '--maf', type = 'character', help = 'input maf file (somatic and/or germline)');
parser$add_argument('-v', '--mavis', type = 'character', help = 'input SV file (mavis output)');
parser$add_argument('-t', '--targets', type = 'character', help = 'input bed of target regions',
	default = NULL);
parser$add_argument('-x', '--expected', type = 'character', help = 'input maf of known reversions',
	default = NULL);
parser$add_argument('-s', '--splice_ai', type = 'character', help = 'path to vep+spliceAI output',
	default = NULL);
parser$add_argument('-a', '--triallelic', type = 'character', help = 'path to triallelic snps (maf-like)',
	default = NULL);

arguments <- parser$parse_args();

# move to working directory
setwd(arguments$directory);

### TESTING ######################################
if (is.null(arguments$project)) {
	arguments$project <- 'EVOLVE_ctDNA';
	arguments$expected <- 'known_exome_reversions.maf';
	arguments$targets <- '/cluster/projects/ovgroup/projects/OV_Superset/EVOLVE/ctDNA/pipeline_suite/CHARM-MMR_plus_EVOLVE_hg38.bed';

# ctDNA input
# somatic <- '../Ensemble_calls/20220713_EVOLVE_ctdna_plus_exome__mutation_data_cbioportal.txt';
# germline <- '../HaplotypeCaller/CPSR/2022-06-15_EVOLVE_ctDNA_mutations_for_cbioportal.tsv';

	arguments$maf <- '../HaplotypeCaller/cohort/VCF2MAF/2022-07-15_EVOLVE_ctDNA_mutations_for_cbioportal.tsv';
	arguments$mavis <- 'filtered_mavis_output.txt';
	arguments$splice_ai <- '2022-07-14_EVOLVE_ctDNA_ensemble_mutations_SPLICEAI.txt';
	arguments$triallelic <- '../bam_readcount/2022-07-14_EVOLVE_ctDNA_triallelic_SNPs.tsv';

	# exome input
	arguments$query <- 'EVOLVE_exome_mutations.txt';
	}
##################################################

# what's the date?
date <- Sys.Date();

# read in list of known reversions (ie, true positives)
key.maf.fields <- unlist(strsplit('Tumor_Sample_Barcode,Hugo_Symbol,Chromosome,Start_Position,End_Position,Variant_Classification,Variant_Type,Reference_Allele,Tumor_Seq_Allele1,Tumor_Seq_Allele2',','));

if (!is.null(arguments$expected)) {
	known.reversions.full <- read.delim(arguments$expected, comment.char = '#');
	known.reversions <- known.reversions.full[,key.maf.fields];
	}

# indicate key gene lists
hr.geneset <- c('ATM', 'ATR', 'EMSY', 'BRIP1', 'CHEK1', 'CHEK2', 'RAD50', 'RAD51C', 'RAD51D', 'BRIP1', 'BARD1', 'RAD51B', 'FAM175A', 'NBN', 'PALB2', 'PTEN', 'MRE11', 'BRCA1', 'BRCA2');
target.geneset <- c('BRCA1','BRCA2','PALB2','ABCB1','TP53','CCNE1');

## Format target regions
# read data
target_bed <- read.delim(arguments$targets, header = FALSE);
colnames(target_bed)[1:4] <- c('Chromosome','Start','End','ID');

# for now, do not include mutations in MSI and Agena (genotyping) sites
target_bed <- subset(target_bed,!grepl("MSI|Agena",ID) | grepl("HBOC",ID))

# and add some padding
target_bed$Start <- target_bed$Start - 100;
target_bed$End <- target_bed$End + 100;

### MAIN ###########################################################################################
# format the query dataset (exome SNV/INDELs)
query.data <- format_queryset(
	data = arguments$query,
	genes = unique(c(hr.geneset, target.geneset)),
	is.maf = TRUE
	);

# format the test dataset (ctDNA SNV/INDEL/SVs)
test.data <- format_testset(
	maf = arguments$maf,
	mavis = arguments$mavis,
	genes = unique(c(hr.geneset, target.geneset))
	);

### this section will be project specific!
# correct sample IDs (mavis changes all underscores (_) to dashes (-))
query.data$Tumor_Sample_Barcode <- gsub('_','-', query.data$Tumor_Sample_Barcode);
test.data$Tumor_Sample_Barcode <- gsub('_','-', test.data$Tumor_Sample_Barcode);

# fill in patient IDs (since tumour ID's don't match between exome and ctDNA)
query.data$Patient <- substr(query.data$Tumor_Sample_Barcode, 0, 11);
test.data$Patient <- substr(test.data$Tumor_Sample_Barcode, 0, 11);

# add in mutation status (ie, germline or somatic)
test.data$Mutation_Status <- 'unknown';

#########################################

# find overlaps; this outputs a list of data objects
overlap.data <- find_candidate_reversions(
	query = query.data[,colnames(test.data)],
	data = test.data,
	targets = target_bed,
	query.sample.field = 'Patient',
	test.sample.field = 'Patient'
	);

reversion.candidates <- overlap.data$reversion.data;

#### SPLICEAI #####
if (!is.null(arguments$splice_ai)) {
	splice.ai <- read.delim(arguments$splice_ai, skip = 42);
	splice.ai <- splice.ai[grepl('SpliceAI_cutoff=PASS', splice.ai$Extra),];

	splice.ai$Chromosome <- sapply(splice.ai$Location, function(i) { unlist(strsplit(as.character(i),':'))[1] } );
	splice.ai$Start_Position <- sapply(splice.ai$Location, function(i) { unlist(strsplit(as.character(i),':'))[2] } );

	splice.ai$SPLICEAI <- sapply(
		splice.ai$Extra,
		function(i) {
			ai.stats <- unlist(strsplit(as.character(i),';'));
			ai.stats <- ai.stats[grepl('SpliceAI_pred', ai.stats)];
			ai.stats <- unlist(strsplit(ai.stats,'\\|'))[2:5];
			idx <- which(as.numeric(ai.stats) > 0.2);
			if (length(idx) == 0) { return(NA); } else {
				labels <- c('AG','AL','DG','DL')[idx];
				ai.stats <- ai.stats[idx];
				to.return <- paste(paste(labels,ai.stats, sep = '='), collapse = ';');
				return(paste0('SpliceAI:',to.return));
				}
			}
		);

	alt.test.data <- merge(
		test.data, 
		unique(splice.ai[!is.na(splice.ai$SPLICEAI),c('Chromosome','Start_Position','SPLICEAI')])
		);

	alt.test.data$Variant_Classification <- apply(
		alt.test.data[,c('Variant_Classification','SPLICEAI')],
		1,
		paste,
		collapse = ','
		);

	# find overlaps; this outputs a list of data objects
	alt.overlap.data <- find_candidate_reversions(
		query = query.data[which(query.data$Mutation_Status == 'germline'),colnames(test.data)],
		data = alt.test.data[,colnames(test.data)],
		targets = target_bed,
		query.sample.field = 'Patient',
		test.sample.field = 'Patient',
		window.size = 10000
		)$reversion.data;

	if (exists('alt.overlap.data')) {
		reversion.candidates <- rbind(
			overlap.data$reversion.data,
			alt.overlap.data
			);
		}
	}

#### TRIALLELIC VARIANTS #####
if (!is.null(arguments$triallelic)) {

	tri.snps <- read.delim(arguments$triallelic);
	tri.snps$Tumor_Sample_Barcode <- gsub('_','-', tri.snps$Tumor_Sample_Barcode);
	if (is.logical(tri.snps$Tumor_Seq_Allele2)) {
		tri.snps$Tumor_Seq_Allele2 <- 'T';
		}

	alt.test.data <- merge(
		test.data,
		tri.snps
		);

	alt.test.data$Event_ID <- paste0(
		alt.test.data$Event_ID,
		'>',
		alt.test.data$Tumor_Seq_Allele3
		);

	alt.test.data$Variant_Classification <- 'triallelic_pos';
	alt.test.data$Mutation_Status <- 'somatic';

	# find overlaps; this outputs a list of data objects
	alt.overlap.data <- find_candidate_reversions(
		query = query.data[which(query.data$Mutation_Status == 'germline'),colnames(test.data)],
		data = alt.test.data[,colnames(test.data)],
		targets = target_bed,
		query.sample.field = 'Patient',
		test.sample.field = 'Patient',
		window.size = 10000
		)$reversion.data;

	if (exists('alt.overlap.data')) {
		reversion.candidates <- rbind(
			overlap.data$reversion.data,
			alt.overlap.data
			);
		}
	}

# check for known status
if (!is.null(arguments$expected)) {
	reversion.candidates$Status <- 'novel';
	for (i in 1:nrow(known.reversions)) {
		patient <- substr(known.reversions[i,]$Tumor_Sample_Barcode, 0, 11);
		chr <- as.character(known.reversions[i,]$Chromosome);
		pos <- known.reversions[i,]$Start_Position;

		idx.test <- which(reversion.candidates$Patient == patient &
			reversion.candidates$Chromosome == chr &
			reversion.candidates$Start_Position.test <= pos &
			reversion.candidates$End_Position.test >= pos
			);

		if (length(idx.test) > 0) {
			reversion.candidates[idx.test,]$Status <- 'previous evidence';
			}
		}
	}

# annotate candidates with potential functions
annotated.data <- annotate_candidates(reversion.candidates);
if (any(grepl('SpliceAI',annotated.data$Variant_Classification.test))) {
	annotated.data[grepl('SpliceAI',annotated.data$Variant_Classification.test),]$Comment <- 'SpliceAI';
	}
if (any(grepl('triallelic',annotated.data$Variant_Classification.test))) {
	annotated.data[grepl('triallelic',annotated.data$Variant_Classification.test),]$Comment <- 'triallelic';
	}

# add to final output list
this.order <- c(1,7,2:6,8:10,14,11:13,15:ncol(annotated.data));
overlap.data$candidates.annotated <- annotated.data[,this.order];

# save reversions list to file
write.table(
	overlap.data$candidates.annotated,
	file = generate.filename(arguments$project,'possible_reversions','tsv'),
	row.names = FALSE,
	col.names = TRUE,
	sep = '\t'
	);

save(
	overlap.data,
	file = generate.filename(arguments$project,'findReversions','RData')
	);

### SAVE SESSION INFO ##############################################################################
save.session.profile(generate.filename(arguments$project, 'FindReversions_SessionInfo', 'txt'));
