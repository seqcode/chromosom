# ChromoSOM
__A Self-Organizing Map for analyzing chromatin interactions__

## Installation

_Dependencies:_ To run ChromoSOM, you need Java version 1.8+. To build the code from scratch, you will additionally need ant installed. 

_Building/getting executables:_ 
To download the code and build your own executable jar file:
```{r, engine='sh', count_lines}
git clone https://github.com/seqcode/chromosom.git
cd chromosom
ant build
ant makechromosomjar
```
Alternatively, you can download a pre-built JAR file from the Releases page on the github repo. 


## Running ChromoSOM

To run ChromoSOM, use the following:
```{r, engine='sh', count_lines}
java -Xmx20G -jar chromosom.jar {mode} {options}
```
In the above, the “-Xmx20G” argument tells java to use up to 20GB of memory. The {mode} argument can be _train_, _view_, or _use_, reflecting ChromoSOM's three main functionalities: training new maps, viewing the maps (GUI), and analyzing the distribution of data on a trained map, respectively. 

## Train

Trains a new SOM using a Hi-C interaction matrix.  

*Arguments:*

* --in <filename> : input Hi-C interaction matrix filename
* --out <filename> : output SOM filename
* --xnodes <int> : number of nodes on X-axis of SOM grid (default=50)
* --ynodes <int> : number of nodes on Y-axis of SOM grid (default=50)
* --iter <int> : number of SOM training iterations (default=1000)
* --kvarmax <double> : kernel variance max (default=1.2)
* --kvarmin <double> : kernel variance min (default=0.2)
* --numsom <int> : number of SOMs to train (default=10)
* --np <int> : number of processors to use during training (default=1)
* --sim <pearson/cosine> : similarity metric to use during training (default=pearson)

An example command line prompt for training a new map is as follows:

```{r, engine='sh', count_lines}
java -Xmx78G -jar chromosom.jar train --np 20 --xnodes 50 --ynodes 50 --kvarmax 1.2 --kvarmin 0.2 --iter 1000 --numsom 10 --sim pearson --in "GM12878_combined_250kb.oe.data" --out "GM12878_250kb_50x50_pearson"
```

The above example will train a 50x50 map (2500 total nodes) named "GM12878_250kb_50x50_pearson" based on the data in the "GM12878_combined_250kb.oe.data" file. The input file must take the form of a path from the current directory to the matrix file. Training will proceed for 1000 iterations, with a kernel variance shrinking from 1.2 to 0.2, and using Pearson similarity. Ten maps will be trained for quality control (the best SOM is kept), and training will occupy 20 processors (or the maximum available). 

The input Hi-C interaction matrix is expected to be a tab-separate flat-file containing a genome-wide interaction matrix. The data should be in a format similar to the following, with valid chromosomal coordinates as the column and row labels:

```{r, engine='sh', count_lines}
	chr1:500001-750000	chr1:750001-1000000	chr1:1000001-1250000	...
chr1:500001-750000	0	-0.103371206779113	-0.426276348661793	...
chr1:750001-1000000	-0.103371206779113	0	0.618833415223183	...
chr1:1000001-1250000	-0.426276348661793	0.618833415223183	0	...
...
```

## View

View a trained SOM using an interactive GUI.

*Arguments:*

* --useweights: use weights defined in BED file score column

```{r, engine='sh', count_lines}
java -Xmx5G -jar chromosom.jar view
```

This prompt will launch a GUI interface which will allow the user to find and open SOM file. The map display will then open allowing for interaction with the map.


## Use

Project genomic activities (defined by BED coordinate files) to a trained SOM. Generates modified Gini coefficients and other statistics about the distribution of data on the SOM. 

*Arguments:*

* --som <filename>: trained SOM filename
* --searchdir <directory name>: directory containing search files
* --out <output directory>: output directory
* --correlation: run correlation analysis between all pairs of datasets
* --pairwise: run pairwise projections for all pairs of datasets
* --useweights: use weights defined in BED file score column

```{r, engine='sh', count_lines}
java -Xmx38G -jar chromosom.jar use --som GM12878_250kb_50x50_pearson.som --searchdir data --out results
```

In the above, "GM12878_250kb_50x50_pearson.som" is a path to the file containing a trained SOM. The "data" directory contains one or many BED files for analysis. All BED files in the directory will be analyzed. Analysis will produce several files in the "results" directory, including a text file summarizing the modified Gini coefficients, and an image of the output grid data distribution, for each data file in the indicated directory.


## Example analysis

You can reproduce most of the figures and analyses from the manuscript by navigating to the "example" directory in this repo and executing the run.sh script. Note that you will need the following prerequisties for aspects of the analyses and production of figures:
* Java v1.8+
* ChromoSOM JAR file, added to PATH
* MultiMDS [https://github.com/seqcode/multimds]
* Python v2.7+
* Python libraries:
	* numpy
	* matplotlib
	* seaborn
	
