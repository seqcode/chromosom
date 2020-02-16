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

mkdir -p hic_${CELL_TYPE}_${RES_STRING}
cd hic_${CELL_TYPE}_${RES_STRING}

if [ ! -d ${CELL_TYPE}/${RES_STRING}_resolution_intrachromosomal ]
	then
		if [ ! -e ../GSE63525_$CELL_TYPE"_intrachromosomal_contact_matrices".tar.gz ]
			then
				curl ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525_$CELL_TYPE"_intrachromosomal_contact_matrices".tar.gz -o ../GSE63525_$CELL_TYPE"_intrachromosomal_contact_matrices".tar.gz
		fi
	tar xzf ../GSE63525_$CELL_TYPE"_intrachromosomal_contact_matrices".tar.gz ${CELL_TYPE}/${RES_STRING}_resolution_intrachromosomal
fi

#if [ $CELL_TYPE == "K562" ] || [ $CELL_TYPE == "KBM7" ]
#	then
#		CHROMS=(1 2 3 4 5 6 7 8 10 11 12 13 14 15 16 17 18 19 20 21)	#skip translocated
#	else
#		CHROMS=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X)
#fi
CHROMS=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X)

for CHROM in ${CHROMS[*]}
do
	if [ -d ${CELL_TYPE}/${RES_STRING}_resolution_intrachromosomal/chr$CHROM ] && [ ! -e ${CELL_TYPE}_${CHROM}_${RES_STRING}.bed ]
		then
			echo $CHROM
			python ../normalize.py $CELL_TYPE $RES $CHROM
	fi
done
	
cd ..
