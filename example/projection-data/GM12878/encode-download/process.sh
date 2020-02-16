perl process.pl dnase.data ../chromatin-marks
perl process.pl histone-marks.data ../chromatin-marks
perl process.pl tf.data ../tfs

#Remove some files due to extreme counts

#This DNase-seq set has too many peaks (>400k)
rm ../chromatin-marks/DNase-seq__1_ENCFF804BNU.bed
#These TF sets have <1000 peaks
rm ../tfs/ChIP-seq_CBX3-human_1-2_ENCFF607ZXG.bed
rm ../tfs/ChIP-seq_BRCA1-human_1-2_ENCFF686JJM.bed
rm ../tfs/ChIP-seq_ZZZ3-human_1-2_ENCFF078ONC.bed
rm ../tfs/ChIP-seq_NR2C2-human_1-2_ENCFF208TMB.bed
rm ../tfs/ChIP-seq_CHD4-human_1-2_ENCFF428HSL.bed
rm ../tfs/ChIP-seq_YBX1-human_1-2_ENCFF506ECO.bed
rm ../tfs/ChIP-seq_UBTF-human_1-2_ENCFF499CTR.bed
rm ../tfs/ChIP-seq_UBTF-human_1-2_ENCFF791WCX.bed

mkdir ../marks-vs-compartments
cp ../chromatin-marks/*bed ../marks-vs-compartments/
cp ../compartments/A.bed ../marks-vs-compartments/
cp ../compartments/B.bed ../marks-vs-compartments/


mkdir ../all
cp ../tfs/*bed ../all
cp ../compartments/*bed ../all
cp ../chromatin-marks/*bed ../all
