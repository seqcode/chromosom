wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525_GM12878_subcompartments.bed.gz
gunzip GSE63525_GM12878_subcompartments.bed.gz

grep A1 GSE63525_GM12878_subcompartments.bed >A1.bed
grep A2 GSE63525_GM12878_subcompartments.bed >A2.bed
grep B1 GSE63525_GM12878_subcompartments.bed >B1.bed
grep B2 GSE63525_GM12878_subcompartments.bed >B2.bed
grep B3 GSE63525_GM12878_subcompartments.bed >B3.bed
grep B4 GSE63525_GM12878_subcompartments.bed >B4.bed
grep NA GSE63525_GM12878_subcompartments.bed >NA.bed
rm GSE63525_GM12878_subcompartments.bed

cat A1.bed A2.bed >A.bed
cat B1.bed B2.bed B3.bed B4.bed >B.bed

mkdir superonly
cp A.bed superonly/
cp B.bed superonly/
