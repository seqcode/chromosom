bigBedToBed ideas_roadmap127_2019.114.bb ideas_roadmap127_2019.114.bed

for ((i=1; i<=41; i++)); do perl ../../../scripts/grepCol.pl ideas_roadmap127_2019.114.bed 3 $i -w >"state"$i".bed"; done

#The following states have fewer than 1000 instances
rm state23.bed state24.bed state18.bed state10.bed state14.bed state11.bed

rm ideas_roadmap127_2019.114.bed
