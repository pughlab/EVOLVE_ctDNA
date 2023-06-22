setwd('/cluster/projects/ovgroup/projects/OV_Superset/EVOLVE/ctDNA/pipeline_suite/');

# get cnv data (for CCNE1 amplifications)
cnvkit <- read.delim('CNVKit/cnv_calls_with_purity__ci/2022-09-12_EVOLVE_ctDNA__cnvkit_calls.tsv');
cnvkit <- cnvkit[which(cnvkit$gene == 'CCNE1'),];
cnvkit$Sample <- gsub('_ctDNA','', cnvkit$Sample);

mops <- read.delim('panelCNmops/sp_v1/2022-10-20_EVOLVE_ctDNA__summarized_CCNE1_calls.tsv');
mops$Sample <- gsub('_ctDNA','',mops$Sample);

# get fragment info
setwd('fragment_size/per_interval');

short <- read.delim('2022-09-29_EVOLVE_ctDNA__ratio_short_per_interval.tsv')
long <- read.delim('2022-09-29_EVOLVE_ctDNA__ratio_long_per_interval.tsv')

wxs.with.ccne1.amps <- c('EVO-009-004','EVO-009-006','EVO-009-007','EVO-009-009','EVO-009-011','EVO-400-007','EVO-400-008');

# overall
tmp <- apply(short[grepl('CCNE1',short$Name), grep('EVO',colnames(short))],2,median)
tmp2 <- apply(long[grepl('CCNE1',long$Name), grep('EVO',colnames(long))],2,median)
my.data <- merge(tmp,tmp2,by = 'row.names')

colnames(my.data) <- c('Sample','short','long')
my.data$Sample <- gsub('_ctDNA','',my.data$Sample)
my.data$Sample <- gsub('\\.','-',my.data$Sample)

my.data$Expected <- NA;
my.data[grepl('Screening|C1D1',my.data$Sample),]$Expected <- 0;
for (smp in wxs.with.ccne1.amps) {
	idx <- which(grepl(smp, my.data$Sample) & grepl('Screening|C1D1',my.data$Sample))
	if (length(idx) > 0) { my.data[idx,]$Expected <- 1; }
	}

my.data[,c('CNVKit','MOPS')] <- 0;
my.data[which(my.data$Sample %in% mops[which(mops$Call == 1),]$Sample),]$MOPS <- 1;
my.data[which(my.data$Sample %in% cnvkit[which(cnvkit$log2 > 2),]$Sample),]$CNVKit <- 1;


## do some tests
# compare known AMP vs other (WES)
wilcox.test(
	my.data[which(my.data$Expected == 1),]$short,
	my.data[which(my.data$Expected == 0),]$short
	);	# no difference

wilcox.test(
	my.data[which(my.data$Expected == 1),]$long,
	my.data[which(my.data$Expected == 0),]$long
	);	# no difference

# compare cnvkit results
wilcox.test(
	my.data[which(my.data$CNVKit == 1),]$short,
	my.data[which(my.data$CNVKit == 0),]$short
	);	# p = 0.003406

wilcox.test(
	my.data[which(my.data$CNVKit == 1),]$long,
	my.data[which(my.data$CNVKit == 0),]$long
	);	# no difference

# compare mops results
wilcox.test(
	my.data[which(my.data$MOPS == 1),]$short,
	my.data[which(my.data$MOPS == 0),]$short
	);	# 0.048

wilcox.test(
	my.data[which(my.data$MOPS == 1),]$long,
	my.data[which(my.data$MOPS == 0),]$long
	);	# no difference

# per interval
colnames(short) <- gsub('_ctDNA','',colnames(short))
colnames(short) <- gsub('\\.','-',colnames(short))

tmp <- short[grepl('CCNE1',short$Name),]
rownames(tmp) <- tmp$Name
tmp <- tmp[,grepl('EVO',colnames(tmp))]

ratio.data <- data.frame(t(tmp))
test.data <- merge(my.data[,c('Sample','Expected','CNVKit','MOPS')], ratio.data, by.x = 'Sample', by.y = 'row.names');

df <- data.frame(Interval = colnames(test.data)[5:ncol(test.data)], p = NA, Median.A = NA, Median.B = NA);
for (i in 1:nrow(df)) {
	interval <- df[i,]$Interval;
	df[i,]$p <- wilcox.test(test.data[which(test.data$CNVKit == 1),interval], test.data[which(test.data$CNVKit == 0),interval])$p.value;
	df[i,]$Median.A <- median(test.data[which(test.data$CNVKit == 1),interval]);
	df[i,]$Median.B <- median(test.data[which(test.data$CNVKit == 0),interval]);
	}

df$Diff <- df$Median.B - df$Median.A




