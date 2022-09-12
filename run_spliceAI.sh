# load required modules
module load perl
module load samtools
module load tabix
module load vep

# move to directory
cd /cluster/projects/ovgroup/projects/OV_Superset/EVOLVE/ctDNA/pipeline_suite/spliceAI

# first, convert MAF to VCF
perl /cluster/projects/pughlab/bin/vcf2maf-1.6.17/maf2vcf.pl --input-maf 2022-07-13_EVOLVE_ctDNA_ensemble_mutation_data.tsv --output-dir . --output-vcf 2022-07-13_EVOLVE_ctDNA_ensemble_mutation_data.vcf --per-tn-vcfs --ref-fasta /cluster/tools/data/genomes/human/hg38/iGenomes/Sequence/WholeGenomeFasta/genome.fa

# run VEP with SpliceAI plugin
vep -i 2022-07-13_EVOLVE_ctDNA_ensemble_mutation_data.vcf --plugin SpliceAI,snv=/cluster/projects/pughlab/references/SpliceAI/spliceai_scores.raw.snv.hg38.vcf.gz,indel=/cluster/projects/pughlab/references/SpliceAI/spliceai_scores.raw.indel.hg38.vcf.gz,cutoff=0.2 --fasta /cluster/tools/data/genomes/human/hg38/iGenomes/Sequence/WholeGenomeFasta/genome.fa --offline --assembly GRCh38 --dir /cluster/projects/pughlab/references/VEP_cache/98 --force_overwrite -o 2022-07-14_EVOLVE_ctDNA_ensemble_mutations_SPLICEAI.txt

