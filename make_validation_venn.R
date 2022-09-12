library(VennDiagram);

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


setwd('/cluster/projects/ovgroup/projects/OV_Superset/EVOLVE/ctDNA/pipeline_suite/Ensemble_calls/wxs_comparison');

overlap.data <- read.delim('2022-09-08_EVOLVE_ctDNA__validation_data__overlaps.tsv');

venn.diagram(
	x = list(
		exome = which(overlap.data$exome == 1),
		eot = which(overlap.data$eot == 1),
		baseline = which(overlap.data$baseline == 1),
		"on-trial" = which(overlap.data$cycle2 == 1)
		),
	cex = 1.8,
	cat.cex = 1.8,
	filename = generate.filename('EVOLVE_ctDNA','_validation_venn','png'),
	resolution = 200,
	units = 'in',
	height = 6,
	width = 6
	);



