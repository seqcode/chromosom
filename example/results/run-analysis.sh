
#Chromosomes
java -Xmx38G -jar chromosom.jar use --som ../som/GM12878_250kb_50x50_pearson.som --searchdir ../projection-data/GM12878/chromosomes --out chromosomes

#A vs B compartments
java -Xmx38G -jar chromosom.jar use --pairwise --som ../som/GM12878_250kb_50x50_pearson.som --searchdir ../projection-data/GM12878/compartments/superonly --out compartments_AB

#A vs B compartments
java -Xmx38G -jar chromosom.jar use --som ../som/GM12878_250kb_50x50_pearson.som --searchdir ../projection-data/GM12878/compartments --out subcompartments_red
java -Xmx38G -jar chromosom.jar use --som ../som/GM12878_250kb_50x50_pearson.som --searchdir ../projection-data/GM12878/compartments --out subcompartments_blue --flipcolors

#Chromatin-marks
java -Xmx38G -jar chromosom.jar use --som ../som/GM12878_250kb_50x50_pearson.som --searchdir ../projection-data/GM12878/chromatin-marks --out chromatin-marks

#TFs
java -Xmx38G -jar chromosom.jar use --som ../som/GM12878_250kb_50x50_pearson.som --searchdir ../projection-data/GM12878/tfs --out tfs
perl ../scripts/histomaker.pl tfs/lorenz_analysis.txt 10 0 1 2 1 1 >tfs/histogram.txt

#Chromatin-states
java -Xmx38G -jar chromosom.jar use --som ../som/GM12878_250kb_50x50_pearson.som --searchdir ../projection-data/GM12878/chromatin-state --out chromatin-state

#Chromatin-marks vs compartments
mkdir marks-vs-compartments
java -Xmx38G -jar chromosom.jar use --correlation --som ../som/GM12878_250kb_50x50_pearson.som --searchdir ../projection-data/GM12878/marks-vs-compartments --out marks-vs-compartments

#TFs vs TFs
mkdir tfs-vs-tfs
java -Xmx38G -jar chromosom.jar use --correlation --som ../som/GM12878_250kb_50x50_pearson.som --searchdir ../projection-data/GM12878/tfs --out tfs-vs-tfs

#all vs all
mkdir all-vs-all
java -Xmx38G -jar chromosom.jar use --correlation --som ../som/GM12878_250kb_50x50_pearson.som --searchdir ../projection-data/GM12878/all --out all-vs-all
