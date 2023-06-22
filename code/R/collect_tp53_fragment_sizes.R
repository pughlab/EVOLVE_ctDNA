### collect_tp53_fragment_sizes.R ##################################################################
# Examine differences in fragment sizes in cfDNA at somatic TP53 mutation sites.

### PREPARE SESSION ################################################################################
# import libraries
library(BoutrosLab.plotting.general);

source('/cluster/home/sprokope/git/analysis/helper_functions/session.functions.R');

working.dir <- '/cluster/projects/ovgroup/projects/OV_Superset/EVOLVE/ctDNA/pipeline_suite/fragment_size';

setwd(working.dir);

### READ DATA ######################################################################################
# get fragment size summary (per mutation)
fragment.summary <- read.delim('2023-02-24_EVOLVE_ctDNA__tp53_mutation_fragment_size_summary__annotated.tsv');

# finally, get fragment sizes for each variant
input.files <- list.files(pattern = 'alt.reads', recursive = TRUE);
input.files <- input.files[!grepl('Screening2', input.files)]; # junk sample (replicate bam for no reason)
input.files <- input.files[!grepl('reversion', input.files)]; # ignore these ones
input.files <- input.files[!grepl('germline', input.files)]; # ignore these ones

### MAIN ###########################################################################################
# store each type of data for a combined/cumulative plot
known.variants <- list();
chip.variants <- list();
other.variants <- list();

# check/plot each TP53 mutation
for (i in 1:nrow(fragment.summary)) {

	smp <- fragment.summary[i,]$Sample;
	variant <- fragment.summary[i,]$MUTATION;

	alt.file <- input.files[grepl(smp, input.files) & grepl(variant, input.files)];
	ref.file <- sub('alt', 'ref', alt.file);

	# get fragment sizes for each mutation
	ref <- read.delim(ref.file, header = FALSE);
	alt <- tryCatch(
		expr = read.delim(alt.file, header = FALSE),
		error = function(e) { NA }
		);

	if (fragment.summary[i,]$STATUS == 'known') {
		known.variants$ref <- c(known.variants$ref, as.numeric(ref$V1));
		known.variants$alt <- c(known.variants$alt, as.numeric(alt$V1));
		} else if (fragment.summary[i,]$STATUS == 'chip') {
		chip.variants$ref <- c(chip.variants$ref, as.numeric(ref$V1));
		chip.variants$alt <- c(chip.variants$alt, as.numeric(alt$V1));
		} else {
		other.variants$ref <- c(other.variants$ref, as.numeric(ref$V1));
		other.variants$alt <- c(other.variants$alt, as.numeric(alt$V1));
		}

	gc();
	}
	
save(
	known.variants, chip.variants, other.variants,
	file = 'EVOLVE_ctDNA__TP53_fragment_sizes.RData'
	);

### SAVE SESSION INFO ##############################################################################
save.session.profile(generate.filename('FSbyMutationSummary','SessionProfile','txt'));
