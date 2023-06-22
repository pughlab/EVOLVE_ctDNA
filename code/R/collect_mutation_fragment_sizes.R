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
maf$Exome <- 0;
maf[!grepl('ctDNA', maf$Tumor_Sample_Barcode),]$Exome <- 1;

# make sure to note which are KNOWN somatic mutation sites
snp.data <- aggregate(
	Exome ~ Patient + Hugo_Symbol + Chromosome + Start_Position + End_Position + Reference_Allele + Tumor_Seq_Allele2,
	maf[,c('Patient','Hugo_Symbol','Chromosome','Start_Position','End_Position','Reference_Allele','Tumor_Seq_Allele2','Exome')],
	max
	);

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
	snp.data[,c('Patient','Hugo_Symbol','MUTATION','Exome')],
	by.x = 'Patient.ID',
	by.y = 'Patient',
	all.x = TRUE
	);

# indicate fields to fill in
stats.data[,c('Total.Fragments','Median.Total')] <- NA;
stats.data[,c('Total.REF','Total.ALT','Median.REF','Median.ALT')] <- NA;
stats.data[,c('Median.FC','p.value','short.ratio','Mean.REF','Mean.ALT')] <- NA;

# read in and process data
all.fragment.data <- list();
known.fragment.data <- list();

for (file in input.files) {

	smp <- unlist(strsplit(basename(file), '_allUNIQUE'))[1];
	variant <- sub('_alt.reads','',unlist(strsplit(basename(file), '_allUNIQUE__'))[2]);

	idx <- which(stats.data$Sample == smp & stats.data$MUTATION == variant);
	if (length(idx) == 0) { next; }

	is.exome <- (stats.data[idx,]$Exome == 1);

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
	if (smp %in% names(all.fragment.data)) {
		all.fragment.data[[smp]]$alt <- c(all.fragment.data[[smp]]$alt, alt);
		all.fragment.data[[smp]]$ref <- c(all.fragment.data[[smp]]$ref, ref);
		} else {
		all.fragment.data[[smp]]$alt <- alt;
		all.fragment.data[[smp]]$ref <- ref;
		}

	# save exome-specific mutation data to per-sample ref/alt list
	if (is.exome & (stats.data[idx,]$Hugo_Symbol %in% c('TP53','BRCA1','BRCA2'))) {
		if (smp %in% names(known.fragment.data)) {
			known.fragment.data[[smp]]$alt <- c(known.fragment.data[[smp]]$alt, alt);
			known.fragment.data[[smp]]$ref <- c(known.fragment.data[[smp]]$ref, ref);
			} else {
			known.fragment.data[[smp]]$alt <- alt;
			known.fragment.data[[smp]]$ref <- ref;
			}
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

# do some quick reformatting for easier plotting later
ref.data <- list();
alt.data <- list();

for (smp in names(known.fragment.data)) {
	ref.data[[smp]] <- known.fragment.data[[smp]]$ref;
	alt.data[[smp]] <- known.fragment.data[[smp]]$alt;
	}

# summarize mutation data (is ALT different from REF?)
ref.v.alt.results <- data.frame(
	Sample = names(all.fragment.data),
	Median.FC.all = NA,
	short.ratio.all = NA,
	p.value.all = NA,
	Median.FC.known = NA,
	short.ratio.known = NA,
	p.value.known = NA
	);

for (i in 1:nrow(ref.v.alt.results)) {

	smp <- ref.v.alt.results[i,]$Sample;

	# all mutations detected in ctDNA
	median.REF <- median(all.fragment.data[[smp]]$ref);
	median.ALT <- median(all.fragment.data[[smp]]$alt);
	ref.v.alt.results[i,]$Median.FC.all <- median.ALT / median.REF;
	all.reads <- c(all.fragment.data[[smp]]$ref,all.fragment.data[[smp]]$alt);
	ref.v.alt.results[i,]$short.ratio.all <- length(all.reads[which(all.reads >= 90 & all.reads < 151)]) /
		length(all.reads[which(all.reads >= 151 & all.reads < 231)]);
	ref.v.alt.results[i,]$p.value.all <- tryCatch(
		expr = wilcox.test(all.fragment.data[[smp]]$ref, all.fragment.data[[smp]]$alt)$p.value,
		error = function(e) { NA }
		);

	# only known mutations (from WES of tissue)
	if (smp %in% names(known.fragment.data)) {
		median.REF <- median(known.fragment.data[[smp]]$ref);
		median.ALT <- median(known.fragment.data[[smp]]$alt);
		ref.v.alt.results[i,]$Median.FC.known <- median.ALT / median.REF;
		all.reads <- c(known.fragment.data[[smp]]$ref,known.fragment.data[[smp]]$alt);
		ref.v.alt.results[i,]$short.ratio.known <- length(all.reads[which(all.reads >= 90 & all.reads < 151)]) /
			length(all.reads[which(all.reads >= 151 & all.reads < 231)]);
		ref.v.alt.results[i,]$p.value.known <- tryCatch(
			expr = wilcox.test(known.fragment.data[[smp]]$ref, known.fragment.data[[smp]]$alt)$p.value,
			error = function(e) { NA }
			);
		}
	}

# format data to save
mutation.data <- stats.data[,-which(colnames(stats.data) == 'MUTATION')];
colnames(mutation.data)[which(colnames(mutation.data) == 'Group.ctdna')] <- 'Group';
mutation.data$Sample <- gsub('_ctDNA', '', mutation.data$Sample);

names(ref.data) <- gsub('_ctDNA','', names(ref.data));
names(alt.data) <- gsub('_ctDNA','', names(alt.data));

ref.v.alt.results$Sample <- gsub('_ctDNA','', ref.v.alt.results$Sample);

save(
	ref.data,
	alt.data,
	mutation.data,
	ref.v.alt.results,
	file = 'EVOLVE_ctDNA__mutation_fragment_size.RData'
	);

### SAVE SESSION INFO ##############################################################################
save.session.profile(generate.filename('FSbyMutationSummary','SessionProfile','txt'));
