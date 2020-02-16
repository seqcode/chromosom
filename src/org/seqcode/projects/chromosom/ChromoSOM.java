package org.seqcode.projects.chromosom;

import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.TimeZone;

public class ChromoSOM {

	public static String version = "0.1";
	
	String mode=null;
	String trainInFilename = null;
	String trainOutFilename = null;
	int trainXnodes = 50;
	int trainYnodes = 50;
	int trainIter = 1000;
	double trainKvarMax = 1.2;
	double trainKvarMin = 0.2;
	int trainNumSOM = 10;
	int numProc = 1;
	boolean useCosineSim = false;
	String trainSimilarity = "pearson";
	boolean includeIntraChromosomal = true;
	
	boolean useWeights = false;
	
	String somFilename = null;
	String searchDirectory = null;
	String outDirectory = null;
	boolean nodes = false;
	boolean flipColors = false;
	boolean runCorrelation = false;
	boolean runPairwise = false;
	
	
	public ChromoSOM(String [] args, boolean noHelp) {
		ArgParser ap = new ArgParser(args);
		System.err.println("ChromoSOM version: "+ChromoSOM.version+"\n\n");
		
		if(args.length==0 || ap.hasKey("h")){
			if(!noHelp)
				System.err.println(getHelpHeader()+"\n"+getHelpTrain()+"\n"+getHelpView()+"\n"+getHelpUse());			
		}else{
			
			if(args[0].equalsIgnoreCase("train")){		//Train the SOM
				if(args.length==1 || ap.hasKey("h")){
					if(!noHelp)
						System.err.println(getHelpHeader()+"\n"+getHelpTrain());
				}else{
					
					setMode("train");
					
					//Input Matrix file
					trainInFilename = Args.parseString(args, "in", null);
					if(trainInFilename==null){
						System.err.println("Input Hi-C interaction matrix required");
						System.exit(1);
					}
					//Output path
					DateFormat df = new SimpleDateFormat("yyyy-MM-dd-hh-mm-ss");  
				    df.setTimeZone(TimeZone.getTimeZone("EST"));
				    trainOutFilename = Args.parseString(args, "out", "SOM_"+df.format(new Date()));
				    //SOM size
				    trainXnodes =  Args.parseInteger(args, "xnodes", trainXnodes);
				    trainYnodes =  Args.parseInteger(args, "ynodes", trainYnodes);
				    //Iterations
				    trainIter = Args.parseInteger(args, "iter", trainIter);
				    //kvar
				    trainKvarMax = Args.parseDouble(args, "kvarmax", trainKvarMax);
				    trainKvarMin = Args.parseDouble(args, "kvarmin", trainKvarMin);
				    //numSOM
				    trainNumSOM =  Args.parseInteger(args, "numsom", trainNumSOM);
				    //numProc
				    numProc =  Args.parseInteger(args, "np", numProc);
				    numProc = Math.min(numProc, java.lang.Runtime.getRuntime().availableProcessors());
				    //similarity
				    trainSimilarity = Args.parseString(args, "sim", trainSimilarity);
				    if(trainSimilarity.equals("pearson"))
						useCosineSim = false;
					if(trainSimilarity.equals("cosine"))
						useCosineSim = true;
					
					train();
					
				}
				
			}else if(args[0].equalsIgnoreCase("view")){		//SOM viewer (GUI)
				
				setMode("view");
				
				useWeights = ap.hasKey("useweights");
				
				view();
			
			}else if(args[0].equalsIgnoreCase("compartments")){		//SOM viewer (GUI) - compartments mode
				if(args.length==1 || ap.hasKey("h")){
					if(!noHelp)
						System.err.println(getHelpHeader()+"\n"+getHelpView());
				}else{
					
					setMode("compartments");
					
					somFilename = Args.parseString(args, "som", null);
					if(somFilename==null){
						System.err.println("Trained SOM required");
						System.exit(1);
					}
					searchDirectory = Args.parseString(args, "searchdir", null);
					nodes = ap.hasKey("nodes");
					
					viewcompartments();
				}
				
			}else if(args[0].equalsIgnoreCase("use")){		//Analyze a trained SOM
				if(args.length==1 || ap.hasKey("h")){
					if(!noHelp)
						System.err.println(getHelpHeader()+"\n"+getHelpView());
				}else{
					
					setMode("use");
					
					somFilename = Args.parseString(args, "som", null);
					if(somFilename==null){
						System.err.println("Trained SOM required");
						System.exit(1);
					}
					searchDirectory = Args.parseString(args, "searchdir", null);
					outDirectory = Args.parseString(args, "out", null);
					useWeights = ap.hasKey("useweights");
					flipColors = ap.hasKey("flipcolors");
					runCorrelation = ap.hasKey("correlation");
					runPairwise = ap.hasKey("pairwise");
					
					if(runCorrelation)
						corSOM();
					else if(runPairwise)
						useSOMpairwise();
					else
						useSOM();
				}
			}
		}
	}
	
	//Settors
	public void setMode(String m){mode = m;}
	public void setTrainInFile(String f){trainInFilename = f;}
	public void setTrainOutFile(String f){trainOutFilename = f;}
	public void setTrainXnodes(int i) {trainXnodes=i;}
	public void setTrainYnodes(int i) {trainYnodes=i;}
	public void setTrainIter(int i) {trainIter=i;}
	public void setTrainKvarMax(double i) {trainKvarMax=i;}
	public void setTrainKvarMin(double i) {trainKvarMin=i;}
	public void setTrainNumSOM(int i) {trainNumSOM=i;}
	public void setNumProc(int i) {numProc=i;}
	public void setIncludeIntraChromosomal(boolean s) {includeIntraChromosomal=s;}
	public void setSimilarity(String s){
		if(s.equals("pearson")){
			trainSimilarity = "pearson";
			useCosineSim = false;
		}
		if(s.equals("cosine")){
			trainSimilarity = "cosine";
			useCosineSim = true;
		}
	}
	
	public static void main(String[] args) {
		ChromoSOM som = new ChromoSOM(args, false);
	}

	/**
	 * Train the SOM
	 */
	public void train(){
		SOMTrainer trainer = new SOMTrainer(trainNumSOM, trainXnodes, trainYnodes, trainKvarMax, trainKvarMin, trainIter, useCosineSim, trainOutFilename, trainInFilename, numProc, includeIntraChromosomal);
	    trainer.train();
	}
	/**
	 * SOM viewer GUI
	 */
	public void view(){
		SOMViewer viewer = new SOMViewer(700,700, "view");
		viewer.initViewer(useWeights);
	}
	/**
	 * SOM viewer GUI compartment mode
	 */
	public void viewcompartments(){
		SOMViewer viewer = new SOMViewer(700,700, "view");
		viewer.initCompartments(somFilename, searchDirectory, nodes);
	}
	/**
	 * Analyze a trained SOM
	 */
	public void useSOM(){
		SOMAnalysis analyzer = new SOMAnalysis();
		analyzer.useSOM(somFilename, searchDirectory, outDirectory, !useWeights, flipColors);
	}
	/**
	 * Project pairs of datasets to a trained SOM
	 */
	public void useSOMpairwise(){
		SOMAnalysis analyzer = new SOMAnalysis();
		analyzer.useSOMpairwise(somFilename, searchDirectory, outDirectory, !useWeights, flipColors);
	}
	/**
	 * Run correlation between datasets on a trained SOM
	 */
	public void corSOM(){
		SOMAnalysis analyzer = new SOMAnalysis();
		analyzer.corSOM(somFilename, searchDirectory, outDirectory);
	}
	
	
	
	protected String getHelpHeader(){ 
			return "Copyright (C) Mahony Lab (2017-2020)\n" +
			"\n" +
			"ChromoSOM comes with ABSOLUTELY NO WARRANTY. This is free software, and you are welcome\n"+
			"to redistribute it under certain conditions. See the MIT license for details.\n"+
			"\nUsage:\tchromosom {mode} {options}\n" +
			"\n\tmodes: train/view/use\n";
	}
	protected String getHelpTrain(){ 
			return "chromosom train:\n"+
			"\t--in <filename>: input Hi-C interaction matrix filename\n" +
			"\t--out <filename>: output SOM filename\n" +
			"\t--xnodes <int>: number of nodes on X-axis of SOM grid (default=" + trainXnodes +")\n"+
			"\t--ynodes <int>: number of nodes on Y-axis of SOM grid (default=" + trainYnodes +")\n"+
			"\t--iter <int>: number of SOM training iterations (default=" + trainIter +")\n"+
			"\t--kvarmax <double>: kernel variance max (default=" + trainKvarMax +")\n"+
			"\t--kvarmin <double>: kernel variance min (default=" + trainKvarMin +")\n"+
			"\t--numsom <int>: number of SOMs to train (default=" + trainNumSOM +")\n"+
			"\t--np <int>: number of processors to use during training (default=" + numProc +")\n"+
			//"\t--nointra: flag to not use intrachromosomal data during similarity calculations (default=save files)\n"+
			"\t--sim <pearson/cosine>: similarity metric to use during training (default=" + trainSimilarity +")\n";
	}
	protected String getHelpView(){ 
			return "chromosom view:\n"+
			"\t--useweights: use weights defined in BED file score column\n";
	}
	protected String getHelpViewCompartments(){ 
		return "chromosom compartments:\n"+
			"\t--som <filename>: trained SOM filename\n"+
			"\t--searchdir <directory name>: directory containing search files\n"+
			"\t--nodes: nodes flag\n";
	}
	protected String getHelpUse(){ 
		return "chromosom use:\n"+
			"\t--som <filename>: trained SOM filename\n"+
			"\t--searchdir <directory name>: directory containing search files\n"+
			"\t--out <output directory>: output directory\n" +
			"\t--correlation: run correlation analysis between all pairs of datasets\n" +
			"\t--pairwise: run pairwise projections for all pairs of datasets\n" +
			"\t--useweights: use weights defined in BED file score column\n";
	}
	
}
