# somatic
## Self-Organizing Map for Analysis of Trends in Interacting Chromatin

There are three main functionalities: training new maps, viewing the maps, and using the maps without GUI


## Training:
An example command line prompt for training a new map is as follows:

• java -Xmx38G MultiMap Train 50 50 1.2 0.2 1000 10 1 "My Map" "My Matrix" 40

MultiMap is the package in which the code for training resides. The -Xmx38G argument is
required to handle memory load. This call will produce a 50x50 map (2500 total nodes) named
“My Map” based on the data in the “My Matrix” file. "My Matrix" must take the form of a path from the current directory to the matrix file. Training will proceed for 1000 iterations, with a kernel variance shrinking from 1.2 to 0.2. Ten maps will be trained for quality control, and training will occupy 40 processors (or the maximum available)

## Viewing
An example command line prompts for viewing a trained map is as follows:

• java -Xmx38G BatchMap View

• java -Xmx38G BatchMap View 0 1

The second prompt has two additional integer parameters. The first int refers to the column the program will look at for each projected data set to find the genomic loci, and the second int refers to which column repressents the weights in projected data sets; a value of -1 will assign equal weights to each refernced loci. Executing either of these prompts will launch a GUI interface which will allow the user to find and open SOM file. The map display will then open allowing for interaction with the map.

## Using
An example command line prompts for using a trained map without openning the GUI component is as follows:

• java -Xmx38G BatchMap Use “My Map” “My Data Directory” 0 1 1

The "My Map" argument repressents the path to the file containing a trained SOM. The "My Data Directory" argument is a path to a directory containing one or many data files for analysis. The first int refers to the column the program will look at for each projected data set to find the genomic loci, and the second int refers to which column repressents the weights in projected data sets; a value of -1 will assign equal weights to each refernced loci. The last int takes a value of 1 or 0. A value of 0 will create a text file output names "Lorenze Analysis" in which the Gini coefficients calculated for each data file in the indicated directory will be printed. A value of one will generate a new directory name "Lorenz Analysis, in which the same text file will be generated; additional, an image of the output grid will be saved for each data file.
