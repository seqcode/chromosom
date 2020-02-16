# Downloads and processes Hi-C datasets

## Get 250kb resolution Hi-C data
sh get_hic_intra.sh GM12878_combined 250000
sh get_hic_inter.sh GM12878_combined 250000

## Process the data to make observed/expected matrices
python ../scripts/inter_oe.py GM12878_combined 250kb hic_GM12878_combined_250kb

perl ../scripts/processMatrix.pl ~/group/genomes/hg19/hg19.info 250000 hic_GM12878_combined_250kb/GM12878_combined_250kb_interchromosomal.oe.tsv hic_GM12878_combined_250kb/GM12878_combined_250kb.coords hic_GM12878_combined_250kb/GM12878_combined_250kb.oe.data 0.1 10
perl ../scripts/processMatrix.pl ~/group/genomes/hg19/hg19.info 250000 hic_GM12878_combined_250kb/GM12878_combined_250kb_interchromosomal.contact.tsv hic_GM12878_combined_250kb/GM12878_combined_250kb.coords hic_GM12878_combined_250kb/GM12878_combined_250kb.contact.data 1 10000

