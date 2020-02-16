# ChromoSOM
## ChromoSOM: a Self-Organizing Map for analyzing chromatin interactions

To run ChromoSOM, use the following:
```{r, engine='sh', count_lines}
java -Xmx20G -jar chromosom.jar {mode} {options}
```
In the above, the “-Xmx20G” argument tells java to use up to 20GB of memory. The {mode} argument can be _train_, _view_, or _use_, reflecting ChromoSOM's three main functionalities: training new maps, viewing the maps (GUI), and analyzing the distribution of data on a trained map, respectively. 


## Train

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

java -Xmx38G MultiMap Train 50 50 1.2 0.2 1000 10 1 "My Map" "My Matrix" 40 0

MultiMap is the package in which the code for training resides. The -Xmx38G argument is required to handle memory load. This call will produce a 50x50 map (2500 total nodes) named “My Map” based on the data in the “My Matrix” file. "My Matrix" must take the form of a path from the current directory to the matrix file. Training will proceed for 1000 iterations, with a kernel variance shrinking from 1.2 to 0.2. Ten maps will be trained for quality control, and training will occupy 40 processors (or the maximum available). The final parameter indicates whether a second file should be saved containing the weight vectors for each, this file can be used for further map analysis. A zero in this parameter will save no such file, a 1 will save the file.

## View
An example command line prompts for viewing a trained map is as follows:

• java -Xmx38G BatchTrainer View 0 1

The prompt has two additional integer parameters. The first int refers to the column the program will look at for each projected data set to find the genomic loci. By default, 0 should be used for the output of Somatic. The second int refers to which column represents the weights in projected data sets; a value of -1 will assign equal weights to each referenced loci. Executing either of these prompts will launch a GUI interface which will allow the user to find and open SOM file. The map display will then open allowing for interaction with the map.

## Using
An example command line prompts for using a trained map without opening the GUI component is as follows:

• java -Xmx38G BatchTrainer Use “My Map” “My Data Directory” 0 1 1

The "My Map" argument represents the path to the file containing a trained SOM. The "My Data Directory" argument is a path to a directory containing one or many BED files for analysis. The first int refers to the column the program will look at for each projected data set to find the genomic loci, and the second int refers to which column represents the weights in projected data sets; a value of -1 will assign equal weights to each referenced loci. The last int takes a value of 1 or 0. A value of 0 will create a text file output names "Lorenz Analysis" in which the Gini coefficients calculated for each data file in the indicated directory will be printed. A value of one will generate a new directory name "Lorenz Analysis, in which the same text file will be generated; additional, an image of the output grid will be saved for each data file.
