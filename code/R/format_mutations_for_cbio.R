### format_mutations_for_cbio.R ####################################################################

### PREAMBLE #######################################################################################
setwd('/cluster/projects/ovgroup/projects/OV_Superset/EVOLVE/ctDNA/pipeline_suite/Ensemble_calls');

source('/cluster/home/sprokope/git/analysis/helper_functions/session.functions.R');

# get and format clinvar data
clinvar <- read.delim('/cluster/projects/ovgroup/projects/VENUS/Analysis//ClinVar_20220416__pathogenic_candidates.tsv');
clinvar$Chromosome <- factor(
	clinvar$Chromosome,
	levels = c(1:22,'X','Y'),
	labels = paste0('chr',c(1:22,'X','Y'))
	);
clinvar <- clinvar[!is.na(clinvar$Chromosome),];

clinvar2 <- aggregate(clinvar$CLINSIG, by = list(Chromosome = clinvar$Chromosome, Position = clinvar$Position, dbSNP_RS = clinvar$dbSNP_RS), FUN = function(i) { paste(unique(tolower(i)), collapse = ',') } );
colnames(clinvar2)[4] <- 'CLINVAR';

### READ DATA ######################################################################################
# read in mutation data
ctdna <- read.delim('2022-07-13_EVOLVE_ctDNA_ensemble_mutations__oncokb_annotated.tsv')
exome <- read.delim('/cluster/projects/pughlab/projects/OV_Superset/cBioportal/EVOLVE/mutation_data_extended.txt');

# ensure same columns
maf.fields <- colnames(exome);
ctdna[,setdiff(maf.fields, colnames(ctdna))] <- NA;

# fix this one duplicate sample ID that keeps messing up cbioportal!
ctdna <- ctdna[!grepl('Screening2', ctdna$Tumor_Sample_Barcode),];

# clarify patient/sample IDs
ctdna$Patient <- substr(ctdna$Tumor_Sample_Barcode,0,11);
exome$Patient <- substr(exome$Tumor_Sample_Barcode,0,11);

ctdna$Sample <- ctdna$Tumor_Sample_Barcode;
ctdna$Tumor_Sample_Barcode <- gsub('_allUNIQUE','', gsub('_SSCS','', gsub('_DCS','', gsub('_DCS_SSCS','', ctdna$Sample))));

ctdna$Group <- sapply(ctdna$Sample, function(i) {
        if (grepl('Screening',i)) { unlist(strsplit(as.character(i),'Screening_'))[2] 
		} else if (grepl('EOT', i)) {
		unlist(strsplit(as.character(i),'EOT_'))[2]
		} else {
                unlist(strsplit(as.character(i),'D1_'))[2]
		}
        });

ctdna$Group <- factor(ctdna$Group, levels = c('allUNIQUE','DCS_SSCS','DCS','SSCS'));
ctdna <- ctdna[order(ctdna$Tumor_Sample_Barcode, ctdna$Group, -ctdna$t_depth),];

# quick ctdna filter (remove known germline)
germline <- do.call(rbind, lapply(list.files(path = '/cluster/projects/ovgroup/projects/OV_Superset/Combined/EVOLVE/Exome/HaplotypeCaller/cohort/VCF2MAF/', pattern = 'BC_HaplotypeCaller_annotated.maf$', full.names = T), function(i) { tmp <- read.delim(i, skip = 1); tmp$n_vaf <- tmp$t_alt_count / tmp$t_depth; tmp$Patient <- substr(tmp$Tumor_Sample_Barcode, 0, 11); return(tmp[,c('Patient','Hugo_Symbol','Chromosome','Start_Position','End_Position','n_vaf')]); } ) );

germline <- germline[which(germline$n_vaf > 0.1),];

ctdna <- merge(ctdna, germline, all.x = TRUE);
ctdna$Mutation_Status <- 'somatic';
ctdna[!is.na(ctdna$n_vaf),]$Mutation_Status <- 'germline';

somatic <- ctdna[which(ctdna$Mutation_Status == 'somatic'),c(colnames(exome),'Sample','Group')];

# split out known germline from exome data
germline.keep <- exome[which(exome$Mutation_Status == 'germline'),];
exome <- exome[which(exome$Mutation_Status == 'somatic'),];

### check/apply clinvar updates to germline
x <- merge(germline.keep, clinvar2[,c('dbSNP_RS','CLINVAR')]);
if (nrow(x) > 0) {
	# check each manually and add if match
	x[,c('Hugo_Symbol','Start_Position','Tumor_Seq_Allele1','Tumor_Seq_Allele2','dbSNP_RS','CLIN_SIG','CLINVAR')]
	}
rm(x);
###

# filter exome mutations to minimum VAF
t_vaf <- exome$t_alt_count / exome$t_depth;
exome <- exome[which(t_vaf >= 0.05),];

to.remove <- setdiff(
        which(exome$FLAG.high_pop == TRUE),
        which(exome$VARIANT_IN_ONCOKB == 'True')
        );

exome <- exome[-to.remove,];

# get recurrence data
exome.counts <- aggregate(exome$Tumor_Sample_Barcode, by = list(Patient = exome$Patient, Chromosome = exome$Chromosome, Start_Position = exome$Start_Position, End_Position = exome$End_Position, Variant_Type = exome$Variant_Type), FUN = length);
colnames(exome.counts)[6] <- 'Exome.Count';

ctdna.counts <- aggregate(somatic$Tumor_Sample_Barcode, by = list(Patient = somatic$Patient, Chromosome = somatic$Chromosome, Start_Position = somatic$Start_Position, End_Position = somatic$End_Position, Variant_Type = somatic$Variant_Type), FUN = function(i) { length(unique(i)) });
colnames(ctdna.counts)[6] <- 'ctDNA.Count';

var.counts <- merge(exome.counts, ctdna.counts, all = TRUE);

validated <- merge(somatic, var.counts[!is.na(var.counts$Exome.Count),1:4]);
discovery <- merge(somatic, var.counts[is.na(var.counts$Exome.Count),1:4]);

validated <- validated[order(validated$Chromosome, validated$Start_Position, validated$Sample, validated$Group),];
discovery <- discovery[order(discovery$Chromosome, discovery$Start_Position, discovery$Sample, discovery$Group),];

# add variants detected in both TUMOUR and ctDNA
filtered <- validated[!duplicated(validated[,c('Chromosome','Start_Position','End_Position','Tumor_Sample_Barcode','Variant_Type')]),];

to.write <- rbind(
        exome[,maf.fields],
        filtered[,maf.fields]
        );

# now add all new pathogenic variants detected in ctDNA
filtered <- discovery[!duplicated(discovery[,c('Chromosome','Start_Position','End_Position','Tumor_Sample_Barcode','Variant_Type')]),];

###
x <- merge(filtered, clinvar2[,c('dbSNP_RS','CLINVAR')], all.x = TRUE);
if (nrow(x) > 0) {
	# check each manually and add if match
	#x[which(x$CLIN_SIG == '' & !is.na(x$CLINVAR)),c('Hugo_Symbol','Start_Position','Tumor_Seq_Allele1','Tumor_Seq_Allele2','dbSNP_RS','CLIN_SIG','CLINVAR')];
	x[which(x$CLIN_SIG == '' & x$CLINVAR == 'likely_pathogenic'),]$CLIN_SIG <- 'likely_pathogenic';
	x[which(x$CLIN_SIG == '' & x$CLINVAR == 'pathogenic/likely_pathogenic'),]$CLIN_SIG <- 'pathogenic,likely_pathogenic';
	x[which(x$CLIN_SIG == 'uncertain_significance' & x$CLINVAR == 'pathogenic/likely_pathogenic'),]$CLIN_SIG <- 'pathogenic,likely_pathogenic';
	}
filtered <- x[,colnames(filtered)];
rm(x);
###

pathogenic <- which(sapply(filtered$CLIN_SIG, function(i) {
        j <- unlist(strsplit(as.character(i),','));
        if (any(j %in% c('pathogenic','likely_pathogenic'))) { TRUE } else { FALSE }
        }));

oncokb <- which(filtered$VARIANT_IN_ONCOKB == 'True');

filtered <- filtered[unique(c(pathogenic,oncokb)),];

# keep any position with a VAF > LOD (within each patient)
filtered$t_vaf <- filtered$t_alt_count / filtered$t_depth;
keep.pos <- aggregate(
	t_vaf ~ Patient + Hugo_Symbol + Start_Position + Variant_Classification + Variant_Type + HGVSp_Short,
	filtered,
	max
	);
keep.pos <- keep.pos[which(
	(keep.pos$t_vaf >= 0.005 & keep.pos$Variant_Type == 'SNP') |
	(keep.pos$t_vaf >= 0.05 & keep.pos$Variant_Type != 'SNP')
	),];

filtered <- merge(filtered, keep.pos[,1:6]);

## remove some false positives
remove.fp <- c(
	# MSH3:80675096 large mix of INS and DEL (region is repeating A's)
	which(filtered$Hugo_Symbol == 'MSH3' & filtered$Start_Position == 80675096),
	# MSH3:80675095 large mix of INS and DEL (region is repeating A's)
	which(filtered$Hugo_Symbol == 'MSH3' & filtered$Start_Position == 80675095),
	# TP53:7673796 in EVO-009-011 (very low VAF, called in single BAM per sample, called by M2 only)
	which(filtered$Hugo_Symbol == 'TP53' & filtered$Start_Position == 7673796),
	# MSH2:47414411 in EVO-009-011 (complex region; high VAF for all bases)
	which(filtered$Hugo_Symbol == 'MSH2' & filtered$Start_Position == 47414411),
	# MLH1:37025636 in EVO-009-016 (complex region; high VAF for all bases)
	which(filtered$Hugo_Symbol == 'MLH1' & filtered$Start_Position == 37025636),
	# MSH2:47414282 in EVO-009-020 (very low VAF, called in single BAM per sample, called by M2 only)
	which(filtered$Hugo_Symbol == 'MSH2' & filtered$Start_Position == 47414282),
	# BRCA2:32379886 large mix of INS and DEL (region is repeating A's)
	which(filtered$Hugo_Symbol == 'BRCA2' & filtered$Start_Position == 32379886),
	# BRCA2:32337452 in EVO-400-003 (very low VAF, called by M2 only)
	which(filtered$Hugo_Symbol == 'BRCA2' & filtered$Start_Position == 32337452),
	# BRCA2:32363511 in EVO-400-007 (very low VAF, called by M2 only)
	which(filtered$Hugo_Symbol == 'BRCA2' & filtered$Start_Position == 32363511),
	# BARD1:214792327 in EVO-009-016 (very low VAF, called by M2 only)
	which(filtered$Hugo_Symbol == 'BARD1' & filtered$Start_Position == 214792327),
	# PALB2:23635000 in EVO-400-007 (low VAF, called by M2 only; high recurrence with VERY low vaf)
	which(filtered$Hugo_Symbol == 'BRCA2' & filtered$Start_Position == 32363511)
	);

filtered <- filtered[-remove.fp,];

# contamination between EVO-009-016_C2 and EVO-009-010_C10
remove.contam <- c(
	which(filtered$Hugo_Symbol == 'TP53' & filtered$Start_Position == 7670694 & grepl('EVO-009-016', filtered$Tumor_Sample_Barcode)),
	which(filtered$Hugo_Symbol == 'TP53' & filtered$Start_Position == 7673802 & grepl('EVO-009-010', filtered$Tumor_Sample_Barcode)),
	which(filtered$Hugo_Symbol == 'BRCA1' & filtered$Start_Position == 43047666 & grepl('EVO-009-010', filtered$Tumor_Sample_Barcode))
	);

filtered <- filtered[-remove.contam,];

to.write <- rbind(
	to.write,
        filtered[,maf.fields]
        );

# now add in known pathogenic germline variants
pathogenic <- sapply(germline.keep$CLIN_SIG, function(i) {
        j <- unlist(strsplit(as.character(i),','));
        if (any(j %in% c('pathogenic','likely_pathogenic'))) { TRUE } else { FALSE }
        });

germline.keep <- germline.keep[which(pathogenic),];

ctdna.germline <- merge(
	ctdna[which(ctdna$Mutation_Status == 'germline'),],
	germline.keep[,c('Patient','Chromosome','Start_Position','End_Position','Variant_Type')]
	);

ctdna.germline <- ctdna.germline[order(ctdna.germline$Chromosome, ctdna.germline$Start_Position, ctdna.germline$Sample, ctdna.germline$Group),];

filtered <- ctdna.germline[!duplicated(ctdna.germline[,c('Chromosome','Start_Position','End_Position','Tumor_Sample_Barcode','Variant_Type')]),];

to.write <- rbind(
        to.write,
        germline.keep[,maf.fields],
	filtered[,maf.fields]
        );

# fix for cBioportal
for (i in c('n_depth','n_ref_count','n_alt_count','t_depth','t_ref_count','t_alt_count')) {
	to.write[is.na(to.write[,i]),i] <- '';
	}

# save to file
write.table(
        to.write,
        file = generate.filename('EVOLVE_ctdna_plus_exome_','mutation_data_cbioportal','txt'),
        row.names = FALSE,
        col.names = TRUE,
        sep = '\t',
        quote = FALSE
        );

### SAVE SESSION INFO ##############################################################################
save.session.profile(generate.filename('format_cbio_mutations','SessionProfile','txt'));
