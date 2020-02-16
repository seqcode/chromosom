set -e

CELL_TYPE="K562"
RES=$1

mkdir -p hic_${CELL_TYPE}_${RES}

cd hic_${CELL_TYPE}_${RES}


if [ ! -e GSE63525_K562_interchromosomal_contact_matrices.tar.gz ]
then

    curl ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525_K562_interchromosomal_contact_matrices.tar.gz -o GSE63525_K562_interchromosomal_contact_matrices.tar.gz

fi

RES_KB=$(($RES/1000))

if [ $RES_KB -lt 1000 ]
then

    RES_STRING=$RES_KB"kb"
else

    RES_STRING=$(($RES_KB/1000))"mb"
fi

DIR=$RES_STRING"_resolution_interchromosomal"

if [ ! -e K562_interchromosomal/$DIR ]
then

    tar xzf GSE63525_K562_interchromosomal_contact_matrices.tar.gz K562_interchromosomal/$DIR

fi

cd ..

CHROMS=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X)

for i in `seq 0 $((${#CHROMS[@]}-1))`
do

    CHROM1=${CHROMS[$i]}

    j=0
    while (( $j < $i ))
    do
	CHROM2=${CHROMS[$j]}
	
	python normalize.py hic_${CELL_TYPE}_${RES}/K562_interchromosomal $RES $CHROM1 --chrom2 $CHROM2
	
	mv hic_${CELL_TYPE}_${RES}/K562_interchromosomal_${CHROM2}_${CHROM1}_${RES_STRING}.bed hic_${CELL_TYPE}_${RES}/K562_${CHROM2}_${CHROM1}_${RES_STRING}.bed
	j=$((j+1))

    done
done
