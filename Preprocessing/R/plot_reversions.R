### JIGV ###
# set up
output.dir <- '/cluster/projects/ovgroup/projects/OV_Superset/EVOLVE/ctDNA/pipeline_suite/Reversions/JIGV';
bam.list <- '/cluster/projects/ovgroup/projects/OV_Superset/EVOLVE/ctDNA/pipeline_suite/configs/sample_info.txt';
site.list <- '/cluster/projects/ovgroup/projects/OV_Superset/EVOLVE/ctDNA/pipeline_suite/Reversions/2021-09-07_EVOLVE_ctDNA_possible_reversions.tsv';

# write placeholder command
jIGV_command <- "/cluster/projects/pughlab/bin/jigv_v0.1.8/jigv --sites BEDFILE --fasta /cluster/tools/data/genomes/human/hg38/iGenomes/Sequence/WholeGenomeFasta/genome.fa --annotation /cluster/projects/pughlab/references/jigv_refs/hg38.refGene.bed.gz BAMFILE > OUTFILE";

# read in data
sample.info <- read.delim(bam.list);
reversions <- read.delim(site.list);

# move to output directory
setwd(output.dir);

# initiate file to write to
write('### Commands for JIGV ###', file = 'jIGV_possible_reversions.txt');

# loop over each sample
for (smp in unique(reversions$Sample.test)) {

	patient <- substr(smp,0,11);

	# find the bam
	smp.type <- unlist(strsplit(smp, '-'))[5];
	bam.type <- if (grepl('allUNIQUE', smp)) { 'ALL_UNIQUE';
		} else if (grepl('SSCS', smp) & !grepl('DCS', smp)) { 'SSCS';
		} else if (grepl('DCS', smp) & !grepl('SSCS', smp)) { 'DCS';
		} else if (grepl('SSCS', smp) & grepl('DCS', smp)) { 'DCS_SSCS';
		}

	bam.file <- sample.info[which(sample.info$Patient == patient &
		sample.info$Tumor.Type == smp.type),bam.type];

	# extract the query (ie, exome) site(s)
	query.pos.fields <- c('Chromosome','Start_Position.query','End_Position.query');
	query.site <- reversions[which(reversions$Sample.test == smp),query.pos.fields];

	# change from VCF to BED format
	query.site$Start_Position.query <- query.site$Start_Position.query - 1;
	colnames(query.site) <- c('CHROM','START','END');

	# extract the test (ie, ctDNA) site(s)
	test.pos.fields <- c('Chromosome','Start_Position.test','End_Position.test');
	test.site <- reversions[which(reversions$Sample.test == smp),test.pos.fields];

	# change from VCF to BED format
	test.site$Start_Position.test <- test.site$Start_Position.test - 1;
	colnames(test.site) <- c('CHROM','START','END');

	# combine and sort sites
	site.list <- unique(rbind(query.site, test.site));
	site.list$CHROM <- factor(site.list$CHROM, levels = paste0('chr',c(1:22,'X','Y')));
	site.list <- site.list[order(site.list$CHROM, site.list$START, site.list$END),];

	# write to bed file
	write.table(
		site.list,
		file = paste0(smp, '_reversion_sites.bed'),
		row.names = FALSE,
		col.names = FALSE,
		quote = FALSE,
		sep = '\t'
		);

	# update the command
	new.cmd <- sub('BAMFILE', bam.file, jIGV_command);
	new.cmd <- sub('BEDFILE', paste0(smp, '_reversion_sites.bed'), new.cmd);
	new.cmd <- sub('OUTFILE', paste0(smp, '_reversion_sites_jigv.html'), new.cmd);

	# write it to file
	write(new.cmd, file = 'jIGV_possible_reversions.txt', append = TRUE);

	# run it
	system(new.cmd);
	}


### GVIZ ###
# set up
output.dir <- '/cluster/projects/ovgroup/projects/OV_Superset/EVOLVE/ctDNA/pipeline_suite/Reversions/GVIZ';
bam.list <- '/cluster/projects/ovgroup/projects/OV_Superset/EVOLVE/ctDNA/pipeline_suite/configs/sample_info.txt';
site.list <- '/cluster/projects/ovgroup/projects/OV_Superset/EVOLVE/ctDNA/pipeline_suite/Reversions/2021-09-07_EVOLVE_ctDNA_possible_reversions.tsv';

# write placeholder command
bamSnapShot_command <- "Rscript ~/git/WGS_pipeline/utils/bamSnapShot.R --bam BAMFILE --output_dir . --sites SITE";

# read in data
sample.info <- read.delim(bam.list);
reversions <- read.delim(site.list);

# move to output directory
setwd(output.dir);

# initiate file to write to
write('### Commands for bamSnapShot.R ###', file = 'bamSnapShot_possible_reversions.txt');

# loop over each potential reversion
for (i in 1:nrow(reversions)) {

	# extract patient/sample info
	patient <- as.character(reversions[i,]$Patient);
	smp <- as.character(reversions[i,]$Sample.test);

	# extract the query (ie, exome) site
	query.site <- paste0(
		reversions[i,]$Chromosome,':',
		reversions[i,]$Start_Position.query,'-',
		reversions[i,]$End_Position.query
		);

	# extract the test (ie, ctDNA) site
	reversion.site <- paste0(
		reversions[i,]$Chromosome,':',
		reversions[i,]$Start_Position.test,'-',
		reversions[i,]$End_Position.test
		);

	# find the bam
	smp.type <- unlist(strsplit(smp, '-'))[5];
	bam.type <- if (grepl('allUNIQUE', smp)) { 'ALL_UNIQUE';
		} else if (grepl('SSCS', smp) & !grepl('DCS', smp)) { 'SSCS';
		} else if (grepl('DCS', smp) & !grepl('SSCS', smp)) { 'DCS';
		} else if (grepl('SSCS', smp) & grepl('DCS', smp)) { 'DCS_SSCS';
		}

	bam.file <- sample.info[which(sample.info$Patient == patient &
		sample.info$Tumor.Type == smp.type),bam.type];

	# update the command
	new.cmd <- sub('BAMFILE', bam.file, bamSnapShot_command);
	new.cmd <- sub('SITE', paste0(query.site, ',', reversion.site), new.cmd);

	# write it to file
	write(new.cmd, file = 'bamSnapShot_possible_reversions.txt', append = TRUE);

	# run it
	system(new.cmd);

	# write command to rename files (bam names are long and unclear)
	rename.cmd <- 'rename BAMFILE SAMPLE *pdf';
	rename.cmd <- sub('BAMFILE', basename(bam.file), rename.cmd);
	rename.cmd <- sub('SAMPLE', smp, rename.cmd);

	# write it to file
	write(paste0(rename.cmd, "\n"), file = 'bamSnapShot_possible_reversions.txt', append = TRUE);

	# run it
	system(rename.cmd);
	}
