# ChromoSOM example analysis

## This is the top-level execution script. 
## See the README for prerequisites. All execution scripts assume that chromosom.jar has been added to the PATH. 

## Note that you can skip Hi-C data download and SOM training by unzipping som/GM12878_250kb_50x50_pearson.som.gz and then skipping to step 3) below

## 1) Download the data
cd hic-data
sh run-hic-data.sh
## Make the HiC matrix plots
cd plots
sh make-hic-heatmaps.sh
cd ..

## 2) Train the SOM
cd som
sh trainSOM.sh
cd ..

## 3) Download projection data
cd projection-data/GM12878
cd compartments
sh download-compartments.sh
cd ..
cd encode-download
sh download.sh
sh process.sh
cd ..
cd chromatin-state
sh convert.sh
cd ..
cd ../..

## 4) Project data to the trained SOM and analyze
cd results
sh run-analysis.sh
cd ..
