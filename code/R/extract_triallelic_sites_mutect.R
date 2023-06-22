
mutect.stats <- list.files(pattern = 'stats.gz', recursive = TRUE);

key.stats.fields <- c('contig','position','tumor_name');

triallelic.sites <- list();

for (file in mutect.stats) {

	print(paste("Reading file", file, '...'));
	smp <- as.character(unlist(strsplit(basename(file),'_'))[1]);
	data <- read.delim(file, skip = 1);
	data.subset <- data[grepl('triallelic_site',data$failure_reasons),];
	print(paste(">>Found", nrow(data.subset), "potential triallelic sites."));
	triallelic.sites[[smp]] <- data.subset[,key.stats.fields];
	gc();

	}

combined <- do.call(rbind, triallelic.sites)
colnames(combined) <- c('contig','position','tumor_name');

write.table(
	unique(combined),
	file = '2022-07-13_EVOLVE_ctDNA__possible_triallelic_sites.tsv',
	row.names = FALSE,
	col.names = TRUE,
	sep = '\t'
	);
