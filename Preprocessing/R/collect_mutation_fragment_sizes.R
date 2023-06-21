### collect_mutation_fragment_sizes.R ##############################################################
# Collect fragment sizes from cfDNA at known mutation sites

### PREPARE SESSION ################################################################################
# import libraries
library(BoutrosLab.plotting.general);

source('/cluster/home/sprokope/git/analysis/helper_functions/session.functions.R');

working.dir <- '/cluster/projects/ovgroup/projects/OV_Superset/EVOLVE/ctDNA/pipeline_suite/fragment_size';

setwd(working.dir);

### READ DATA ######################################################################################
# look at fragment size at known mutation sites
maf <- read.delim('../cBioportal/EVOLVE_ctDNA_20220908/mutation_data_extended.txt');

# also get purity estimates
purity <- read.delim('../Ensemble_calls/tumour_content/2022-09-06_EVOLVE_ctDNA__estimated_tumour_content.tsv');

# and results from above analysis (global fragment sizes)
global.fs <- read.delim('2022-09-07_EVOLVE_ctDNA_fragment_size_summary.tsv');
global.fs$Patient <- substr(global.fs$Sample, 0, 11);

# finally, get fragment sizes for each variant
input.files <- list.files(pattern = 'alt.reads', recursive = TRUE);
input.files <- input.files[!grepl('Screening2', input.files)]; # junk sample (replicate bam for no reason)
input.files <- input.files[!grepl('reversion', input.files)]; # ignore these ones
input.files <- input.files[!grepl('germline', input.files)]; # ignore these ones

### FORMAT DATA ####################################################################################
# format mutation data
maf <- maf[which(maf$Mutation_Status == 'somatic'),];
maf <- maf[which(maf$Variant_Type == 'SNP'),];

maf$Patient <- substr(maf$Tumor_Sample_Barcode,0,11);

# first, let's only look at KNOWN somatic mutation sites
snp.data <- unique(maf[!grepl('ctDNA',maf$Tumor_Sample_Barcode),c('Patient','Hugo_Symbol','Chromosome','Start_Position','End_Position','Reference_Allele','Tumor_Seq_Allele2')]);
snp.data <- snp.data[which(snp.data$Hugo_Symbol %in% c('TP53','BRCA1','BRCA2')),];
snp.data$MUTATION <- paste(
	snp.data$Chromosome,
	snp.data$Start_Position,
	snp.data$Reference_Allele,
	snp.data$Tumor_Seq_Allele2,
	sep = '_'
	);

# compare fragment sizes (alt vs ref)
stats.data <- merge(
	purity[,c('Patient.ID','Sample','Group.ctdna')],
	snp.data[,c('Patient','Hugo_Symbol','MUTATION')],
	by.x = 'Patient.ID',
	by.y = 'Patient',
	all.x = TRUE
	);

# indicate fields to fill in
stats.data[,c('Total.Fragments','Median.Total')] <- NA;
stats.data[,c('Total.REF','Total.ALT','Median.REF','Median.ALT')] <- NA;
stats.data[,c('Median.FC','p.value','short.ratio','Mean.REF','Mean.ALT')] <- NA;

# read in and process data
fragment.data <- list();

for (file in input.files) {

	smp <- unlist(strsplit(basename(file), '_allUNIQUE'))[1];
	variant <- sub('_alt.reads','',unlist(strsplit(basename(file), '_allUNIQUE__'))[2]);

	idx <- which(stats.data$Sample == smp & stats.data$MUTATION == variant);
	if (length(idx) == 0) { next; }

	# get fragment sizes for each mutation
	ref <- as.numeric(read.delim(sub('alt','ref',file), header = FALSE)$V1);
	alt <- tryCatch(
		expr = as.numeric(read.delim(file, header = FALSE)$V1),
		error = function(e) { NA }
		);

	# compare ref/alt at each site
	stats.data[idx,c('Total.Fragments')] <- length(ref) + length(alt);
	stats.data[idx,c('Total.REF','Total.ALT')] <- c(length(ref), length(alt));
	stats.data[idx,c('Median.Total')] <- median(c(ref,alt));
	stats.data[idx,c('Median.REF','Median.ALT')] <- c(median(ref), median(alt));
	stats.data[idx,c('Mean.REF','Mean.ALT')] <- c(mean(ref), mean(alt));
	stats.data[idx,]$Median.FC <- median(alt) / median(ref);
	stats.data[idx,]$p.value <- tryCatch(
		expr = wilcox.test(ref, alt)$p.value,
		error = function(e) { NA }
		);

	# find ratio of short:normal reads (overall at site)
	all.reads <- c(ref,alt);
	stats.data[idx,]$short.ratio <- length(all.reads[which(all.reads >= 90 & all.reads < 151)]) /
		length(all.reads[which(all.reads >= 151 & all.reads < 231)]);

	# save data to per-sample ref/alt list
	if (smp %in% names(fragment.data)) {
		fragment.data[[smp]]$alt <- c(fragment.data[[smp]]$alt, alt);
		fragment.data[[smp]]$ref <- c(fragment.data[[smp]]$ref, ref);
		} else {
		fragment.data[[smp]]$alt <- alt;
		fragment.data[[smp]]$ref <- ref;
		}

	gc();
	}

na.counts <- apply(stats.data[,c('Median.REF','Median.ALT')], 1, function(i) { all(is.na(i)) } );
stats.data <- stats.data[!na.counts,];

write.table(
	stats.data,
	file = generate.filename('EVOLVE_ctDNA','_known_mutation_fragment_size_summary','tsv'),
	row.names = FALSE,
	col.names = TRUE,
	sep = '\t'
	);

# do some quick tests
ref.data <- list();
alt.data <- list();

for (smp in names(fragment.data)) {
	ref.data[[smp]] <- fragment.data[[smp]]$ref;
	alt.data[[smp]] <- fragment.data[[smp]]$alt;
	}

# format data to save
mutation.data <- stats.data[,-which(colnames(stats.data) == 'MUTATION')];
colnames(mutation.data)[which(colnames(mutation.data) == 'Group.ctdna')] <- 'Group';
mutation.data$Sample <- gsub('_ctDNA', '', mutation.data$Sample);

names(ref.data) <- gsub('_ctDNA','', names(ref.data));
names(alt.data) <- gsub('_ctDNA','', names(alt.data));

save(
	ref.data,
	alt.data,
	mutation.data,
	file = 'EVOLVE_ctDNA__mutation_fragment_size.RData'
	);

### SAVE SESSION INFO ##############################################################################
save.session.profile(generate.filename('FSbyMutationSummary','SessionProfile','txt'));
