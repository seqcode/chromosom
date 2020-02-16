set -e

CELL_TYPE=$1
RES=$2
RES_KB=$(($RES/1000))
if [ $RES_KB -lt 1000 ]
then

    RES_STRING=$RES_KB"kb"
else

    RES_STRING=$(($RES_KB/1000))"mb"
fi
OUTDIR=${RES_STRING}"_resolution_interchromosomal"

mkdir -p hic_${CELL_TYPE}_${RES_STRING}

cd hic_${CELL_TYPE}_${RES_STRING}

if [ ! -d ${CELL_TYPE}/${OUTDIR} ]
then
	if [ ! -e ../GSE63525_$CELL_TYPE"_interchromosomal_contact_matrices".tar.gz ]
        then
               curl ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525_$CELL_TYPE"_interchromosomal_contact_matrices".tar.gz -o ../GSE63525_$CELL_TYPE"_interchromosomal_contact_matrices".tar.gz
        fi
	tar xzf ../GSE63525_$CELL_TYPE"_interchromosomal_contact_matrices".tar.gz ${CELL_TYPE}"_interchromosomal"/${OUTDIR}
	mkdir -p ${CELL_TYPE}
	mv ${CELL_TYPE}"_interchromosomal"/${OUTDIR} ${CELL_TYPE}/
	rm -R ${CELL_TYPE}"_interchromosomal"
fi


#if [ $CELL_TYPE == "K562" ] || [ $CELL_TYPE == "KBM7" ]
#        then
#                CHROMS=(1 2 3 4 5 6 7 8 10 11 12 13 14 15 16 17 18 19 20 21)    #skip translocated                                               
#        else
#                CHROMS=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X)
#fi
CHROMS=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X)

for i in `seq 0 $((${#CHROMS[@]}-1))`
do

    CHROM1=${CHROMS[$i]}

    j=0
    while (( $j < $i ))
    do
	CHROM2=${CHROMS[$j]}

	echo "$CHROM1 - $CHROM2"	
	python ../normalize.py ${CELL_TYPE} $RES $CHROM1 --chrom2 $CHROM2
	
	j=$((j+1))

    done
done
cd ..
