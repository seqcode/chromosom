package org.seqcode.projects.somatic;


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
	boolean trainSaveWeights = true;
	String trainSimilarity = "pearson";
	
	public ChromoSOM() {
		
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
	public void setTrainSaveWeights(boolean s) {trainSaveWeights=s;}
	public void setTrainSimilarity(String s){trainSimilarity = s;}
	
	
	public static void main(String[] args) {
		System.err.println("ChromoSOM version: "+ChromoSOM.version+"\n\n");
		
	}

	public void printHelp() {
		
		System.err.println("" +
				"Copyright (C) Mahony Lab (2017-2020)\n" +
				"\n" +
				"ChromoSOM comes with ABSOLUTELY NO WARRANTY. This is free software, and you are welcome\n"+
				"to redistribute it under certain conditions. See the MIT license for details.\n"+
				"\nUSAGE:\tchromosom {mode} {options}\n" +
				"\n\tmodes: train/view/use\n" +
				"\nchromosom train:\n"+
				"\t--in <filename>: input Hi-C interaction matrix filename\n" +
				"\t--out <filename>: output SOM filename\n" +
				"\t--xnodes <int>: number of nodes on X-axis of SOM grid (default=" + trainXnodes +")\n"+
				"\t--ynodes <int>: number of nodes on Y-axis of SOM grid (default=" + trainYnodes +")\n"+
				"\t--iter <int>: number of SOM training iterations (default=" + trainIter +")\n"+
				"\t--kvarmax <double>: kernel variance max (default=" + trainKvarMax +")\n"+
				"\t--kvarmin <double>: kernel variance min (default=" + trainKvarMin +")\n"+
				"\t--numsom <int>: number of SOMs to train (default=" + trainNumSOM +")\n"+
				"\t--np <int>: number of processors to use during training (default=" + numProc +")\n"+
				"\t--noweightsfile: flag to not save the SOM weights file (default=save files)\n"+
				"\t--sim <pearson/cosine>: similarity metric to use during training (default=" + trainSimilarity +")\n"+
				"");

	}
}
