# Train a SOM

# resolution and SOM size
res="250kb";
size=50;
# metric: 0=pearson, 1=cosine
metric=1;
metricname="pearson";
if [ $metric = "1" ]
then
	metricname="cosine";
fi

time(java -Xmx78G -jar ../chromosom.jar train --np 20 --xnodes $size --ynodes $size --kvarmax 1.2 --kvarmin 0.2 --iter 1000 --numsom 10 --sim $metric --out "GM12878_"$res"_"$size"x"$size"_"$metricname --in "../hic-data/hic_GM12878_combined_"$res"/GM12878_combined_"$res".oe.data" )
