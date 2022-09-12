### assess_contamination.R #########################################################################
# Compare somatic SNV/INDELs between possible contamination pairs

### PREAMBLE #######################################################################################
library(VennDiagram);
library(BoutrosLab.plotting.general);

setwd('/cluster/projects/ovgroup/projects/OV_Superset/EVOLVE/ctDNA/pipeline_suite/Ensemble_calls');

### READ DATA ######################################################################################
# load in combined mutation calls
load('2022-07-13_EVOLVE_ctDNA_CombinedMutationData.RData');

### FORMAT DATA ####################################################################################
# extract sample IDs
combined.data$Sample <- combined.data$Tumor_Sample_Barcode;
combined.data$Sample <- gsub('_allUNIQUE','',gsub('_DCS','', gsub('_SSCS','', gsub('_DCS_SSCS','',combined.data$Sample))));

combined.data <- combined.data[!grepl('Screening2',combined.data$Tumor_Sample_Barcode),];

# indicate ensemble method
combined.data$ENSEMBLE <- 0;
combined.data$ENSEMBLE.mod <- 0;

# current method takes ALL MuTect2 calls
combined.data[which(combined.data$FILTER == 'PASS'),]$ENSEMBLE.mod <- 1;

# original method takes MuTect2 + low VAF
combined.data[which(combined.data$Count >= 3),]$ENSEMBLE <- 1;
combined.data[which(combined.data$MuTect2 == 1 & combined.data$VAF < 0.1),]$ENSEMBLE <- 1;

passed.variants <- unique(combined.data[which(combined.data$ENSEMBLE == 1),c('Chromosome','Start_Position','End_Position','Reference_Allele','Tumor_Seq_Allele2','Variant_Type','Patient')]);
passed.variants$C4 <- 1;

combined.data <- merge(
	combined.data[,!grepl('C4',colnames(combined.data))],
	passed.variants,
	all.x = TRUE
	);

combined.data[which(combined.data$C4 == 1),]$ENSEMBLE <- 1;

### EXAMINE OVERLAP ################################################################################
# collapse to 1 call per sample per position
tmp <- aggregate(
	ENSEMBLE.mod ~ Sample + Chromosome + Start_Position + End_Position + Variant_Type,
	combined.data,
	FUN = max
	);

tmp$VarID <- paste0(tmp$Chromosome, '_', tmp$Start_Position, '_', tmp$End_Position,'_', tmp$Variant_Type);
tmp <- tmp[which(tmp$ENSEMBLE.mod == 1),];

# compare calls between EVO-009-009_C13D1 and EVO-009-020
venn.diagram(
	x = list(
		'009_bl' = tmp[which(tmp$Sample == 'EVO-009-009_ctDNA_Screening'),]$VarID,
		'020_C7' = tmp[which(tmp$Sample == 'EVO-009-020_ctDNA_C7D1'),]$VarID,
		'009_C13' = tmp[which(tmp$Sample == 'EVO-009-009_ctDNA_C13D1'),]$VarID,
		'020_bl' = tmp[which(tmp$Sample == 'EVO-009-020_ctDNA_Screening'),]$VarID
		),
	cex = 1.8,
	cat.cex = 1.8,
	filename = '2022-07-25_minor_contamination__EVO-009-009-C13D1__EVO-009-020.png',
	fill = c('darkgreen','green','darkgreen','green'),  #'darkorchid4','darkorchid1')
	col = c('darkgreen','green','darkgreen','green'),  #'darkorchid4','darkorchid1')
	alpha = 0.5,
	margin = 0.05,
	resolution = 200,
	units = 'in',
	height = 6,
	width = 6
	);

# compare calls between EVO-009-010_C10D1 and EVO-009-016
venn.diagram(
	x = list(
		'010_bl' = tmp[which(tmp$Sample == 'EVO-009-010_ctDNA_Screening'),]$VarID,
		'016_C2' = tmp[which(tmp$Sample == 'EVO-009-016_ctDNA_C2D1'),]$VarID,
		'010_C10' = tmp[which(tmp$Sample == 'EVO-009-010_ctDNA_C10D1'),]$VarID,
		'016_bl' = tmp[which(tmp$Sample == 'EVO-009-016_ctDNA_Screening'),]$VarID
		),
	cex = 1.8,
	cat.cex = 1.8,
	fill = c('darkorchid4','darkorchid1','darkorchid4','darkorchid1'),
	col = c('darkorchid4','darkorchid1','darkorchid4','darkorchid1'),
	alpha = 0.5,
	filename = '2022-07-25_minor_contamination__EVO-009-010-C10D1__EVO-009-016.png',
	margin = 0.05,
	resolution = 200,
	units = 'in',
	height = 6,
	width = 6
	);

### ALTERNATE
# get germline variants for affected samples
germline <- rbind(
	read.delim('/cluster/projects/ovgroup/projects/OV_Superset/Combined/EVOLVE/Exome/HaplotypeCaller/cohort/VCF2MAF/EVO-009-020-BC_HaplotypeCaller_annotated.maf', skip = 1, stringsAsFactors = FALSE),
	read.delim('/cluster/projects/ovgroup/projects/OV_Superset/Combined/EVOLVE/Exome/HaplotypeCaller/cohort/VCF2MAF/EVO-009-009-BC_HaplotypeCaller_annotated.maf', skip = 1, stringsAsFactors = FALSE),
	read.delim('/cluster/projects/ovgroup/projects/OV_Superset/Combined/EVOLVE/Exome/HaplotypeCaller/cohort/VCF2MAF/EVO-009-010-BC_HaplotypeCaller_annotated.maf', skip = 1, stringsAsFactors = FALSE),
	read.delim('/cluster/projects/ovgroup/projects/OV_Superset/Combined/EVOLVE/Exome/HaplotypeCaller/cohort/VCF2MAF/EVO-009-016-BC_HaplotypeCaller_annotated.maf', skip = 1, stringsAsFactors = FALSE)
	);

germline$Patient <- substr(germline$Tumor_Sample_Barcode, 0, 11);
germline$n_vaf <- germline$t_alt_count / germline$t_depth;
germline <- unique(germline[,c('Patient','Chromosome','Start_Position','End_Position','Hugo_Symbol','Variant_Type','Allele','n_vaf')]);

# compare germline first (are they related perhaps?)
tmp <- merge(
	merge(
		germline[which(germline$Patient == 'EVO-009-009'),-1],
		germline[which(germline$Patient == 'EVO-009-020'),-1],
		by = c('Chromosome','Start_Position','End_Position','Hugo_Symbol','Variant_Type','Allele'),
		suffixes = c('.009','.020'),
		all = TRUE
		),
	merge(
		germline[which(germline$Patient == 'EVO-009-010'),-1],
		germline[which(germline$Patient == 'EVO-009-016'),-1],
		by = c('Chromosome','Start_Position','End_Position','Hugo_Symbol','Variant_Type','Allele'),
		suffixes = c('.010','.016'),
		all = TRUE
		),
	all = TRUE
	);

cor(tmp[,grepl('n_vaf',colnames(tmp))], use = 'pairwise');


# get ctDNA mutations
ctdna <- read.delim('ensemble_mutation_data.tsv');

ctdna <- ctdna[!grepl('Screening2', ctdna$Tumor_Sample_Barcode),]
ctdna$Sample <- gsub('_allUNIQUE','',gsub('_DCS','', gsub('_SSCS','', gsub('_DCS_SSCS','',ctdna$Tumor_Sample_Barcode))));
ctdna$t_vaf <- ctdna$t_alt_count / ctdna$t_depth;

ctdna <- aggregate(
	t_vaf ~ Sample + Chromosome + Start_Position + End_Position + Hugo_Symbol + Variant_Type + Allele,
	ctdna,
	max
	);

combined <- merge(ctdna, germline);




