### find_triallelic_snps.R #########################################################################
# Use germline SNP positions as input for bam-readcount to identify tri-allelic positions

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

### PREPARE SESSION ################################################################################
# import libraries
library(argparse);

# import command line arguments
parser <- ArgumentParser();

parser$add_argument('-d', '--directory', type = 'character', help = 'path to data directory');
parser$add_argument('-p', '--project', type = 'character', help = 'project name');
parser$add_argument('-t', '--targets', type = 'character', help = 'path to target regions', default = NULL);

arguments <- parser$parse_args();

setwd(arguments$directory);

### READ DATA ######################################################################################
# read in target positions (+ any annotations; ie, Sample or Gene names)
if (!is.null(arguments$targets)) {
	anno.data <- read.delim(arguments$targets, header = FALSE);
	} else {
	anno.data <- NULL;
	}

# find files and read in
readcount.files <- list.files(pattern = 'readcounts.tsv', recursive = TRUE);
my.data <- list();
for (i in readcount.files) {
	smp <- sub('_readcounts.tsv','',basename(i));
	tmp.data <- read.delim(i, header = FALSE)[,c(1:4,6:9)];
	my.data[[smp]] <- tmp.data[!grepl(':',tmp.data[,1]),];
	}

# create a single data object
combined.data <- do.call(rbind, my.data);
combined.data$Sample <- sapply(rownames(combined.data),function(i) { unlist(strsplit(i,'\\.'))[1] } );
rownames(combined.data) <- 1:nrow(combined.data);
colnames(combined.data) <- c('Chromosome','Position','Ref','Depth','Base.A','Base.C','Base.G','Base.T','Sample');

combined.data$Depth <- as.numeric(combined.data$Depth);

# extract basecounts
combined.data$Count.A <- sapply(
	combined.data$Base.A, function(i) { unlist(strsplit(as.character(i),':'))[2] } );
combined.data$Count.C <- sapply(
	combined.data$Base.C, function(i) { unlist(strsplit(as.character(i),':'))[2] } );
combined.data$Count.G <- sapply(
	combined.data$Base.G, function(i) { unlist(strsplit(as.character(i),':'))[2] } );
combined.data$Count.T <- sapply(
	combined.data$Base.T, function(i) { unlist(strsplit(as.character(i),':'))[2] } );

combined.data$Count.A <- as.numeric(combined.data$Count.A);
combined.data$Count.T <- as.numeric(combined.data$Count.T);
combined.data$Count.C <- as.numeric(combined.data$Count.C);
combined.data$Count.G <- as.numeric(combined.data$Count.G);

# calculate variant allele frequencies
combined.data$vaf.A <- (combined.data$Count.A / combined.data$Depth);
combined.data$vaf.T <- (combined.data$Count.T / combined.data$Depth);
combined.data$vaf.C <- (combined.data$Count.C / combined.data$Depth);
combined.data$vaf.G <- (combined.data$Count.G / combined.data$Depth);

# try to determine if this is likely a tri-allelic variant (or random sequence error)
combined.data$Allele_Count <- apply(
	combined.data[,grepl('vaf.',colnames(combined.data))],
	1,
	function(i) {
		length(i[which(i > 0.01)]);
		}
	);

# format for output
key.fields <- c('Sample','Chromosome','Position','Allele_Count',colnames(combined.data)[grepl('Count.|vaf.', colnames(combined.data))]);

to.write <- merge(
	anno.data[,c(1,5:7,9:13,16,35:45)],
	droplevels(combined.data[,key.fields]),
	by.x = c('Tumor_Sample_Barcode','Chromosome','Start_Position'),
	by.y = c('Sample','Chromosome','Position'),
	all = FALSE
	);

# save to file
write.table(
	to.write[order(to.write$Chromosome, to.write$Start_Position, to.write$Tumor_Sample_Barcode),],
	file = generate.filename(arguments$project, 'germline_SNP_AlleleCounts','tsv'),
	row.names = FALSE,
	col.names = TRUE,
	sep = '\t'
	);

# extract and format triallelic sites
triallelic <- to.write[which(to.write$Allele_Count == 3),];

triallelic$Tumor_Seq_Allele3 <- NA;
for (i in 1:nrow(triallelic)) {

	vafs <- triallelic[i,grepl('vaf\\.', colnames(triallelic))];
	colnames(vafs) <- gsub('vaf\\.','',colnames(vafs));
	known.alleles <- unique(as.character(
		unlist(triallelic[i,c('Reference_Allele','Tumor_Seq_Allele1','Tumor_Seq_Allele2')])
		));

	vafs <- sort(vafs[,setdiff(colnames(vafs), known.alleles)], decreasing = TRUE);
	triallelic[i,]$Tumor_Seq_Allele3 <- colnames(vafs)[1];
	}

# save to file
write.table(
	triallelic[,c(1:10,31,22:30)],
	file = generate.filename(arguments$project, 'triallelic_SNPs','tsv'),
	row.names = FALSE,
	col.names = TRUE,
	sep = '\t'
	);

