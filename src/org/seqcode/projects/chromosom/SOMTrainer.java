package org.seqcode.projects.chromosom;

import java.util.ArrayList;


public class SOMTrainer
{
	protected int numRounds, xNode, yNode;
	protected double sigma, sigmaStop;
	protected int iterations;
	protected boolean useCosineSim, countIntraChrom;
	protected String outPrefix, matrixFilePath;
	protected int availibleThreads;
	
	public SOMTrainer(int numRounds, int xNode, int yNode, double sigma, double sigmaStop, int iterations, boolean cosineSim, String outPrefix, String matrixFilePath, int availibleThreads, boolean countIntraChrom)
	{
		this.numRounds = numRounds;
		this.xNode = xNode;
		this.yNode = yNode;
		this.sigma = sigma;
		this.sigmaStop = sigmaStop;
		this.iterations = iterations;
		this.useCosineSim = cosineSim;
		this.outPrefix = outPrefix;
		this.matrixFilePath=matrixFilePath;
		this.availibleThreads = availibleThreads;
		this.countIntraChrom = countIntraChrom;
	}
	
	protected void train(){
		ArrayList<ThreadMap> op = new ArrayList<ThreadMap>();
	    for(int i = 0; i< numRounds; i++)
	    {
    	    System.out.println("training");
	    	ThreadMap o = new ThreadMap(xNode, yNode, sigma, sigmaStop, iterations, useCosineSim, outPrefix, matrixFilePath, availibleThreads, countIntraChrom);
		 		o.go();
		 		op.add(o);
	    }
	    qualityControl(op, useCosineSim).writeFile(numRounds);
	}
	
	protected ThreadMap qualityControl(ArrayList<ThreadMap> op, boolean similarity)
	{

		ThreadMap best = op.get(0);
		if(!similarity)
		{
			for(int i = 1; i < op.size(); i++)
			{
				if(op.get(i).pearsonQuality>best.pearsonQuality)
				{
					best = op.get(i);
				}
			}
			op.remove(best);
		}
		else
		{
			for(int i = 1; i < op.size(); i++)
			{
				if(op.get(i).cosineQuality>best.cosineQuality)
				{
					best = op.get(i);
				}
			}
			op.remove(best);
		}
		double affirmatives = 0;
		double total = 0;
		for(int i =0; i< best.nodeNum; i++)
		{
			for(int opi = 0; opi < best.assignments.length; opi++)
			{
				for(int opj = 0; opj < best.assignments.length; opj++)
				if(best.assignments[opi][i]&&best.assignments[opj][i])
				{
					for(int l = 0; l < op.size(); l++)
					{
						total++;
						if(op.get(l).findStability(opi, opj))
						{
							affirmatives++;
						}
					}
				}
			}
		}
		best.stability = affirmatives / total;
		return best;
	}
	public class MultiMapThread extends Thread
	{
		int threadsFor, numberOfSoms;
		ArrayList<ThreadMap> t;
		MultiMapThread(int tt, int nnnn)
		{
			threadsFor = tt;
			numberOfSoms = nnnn;
			t = new ArrayList<ThreadMap>();
		}

	}
	
	
	public class UniMapThread extends Thread
	{
		int threadsFor, numberOfSoms;
		ArrayList<ThreadMap> t;
		UniMapThread(int tt, int nnnn)
		{
			threadsFor = tt;
			numberOfSoms = nnnn;
			t = new ArrayList<ThreadMap>();
		}
	}
	
	//Replacing this main method functionality with function in ChromoSOM main class
	public static void main(String[] args)
	{
		 //to train: args = "train" 'int x' 'int y' 'int sigma' 'double learn rate' 'int iterations' 'int n' 'cosine sim' 'weight by dps' 'exponential decay' 'name' 'matrix file' 'number of threads'
    	if(args[0].equalsIgnoreCase(("Train")))
    	{
    		int xArg = 0;
	    	int yArg = 0;
	    	double sig = 1.2;
	    	int iter = 5000;
	    	double sigStop = 0.4;
	    	int n = 10;
    	    
	    	try {
    	        xArg = Integer.parseInt(args[1]);
    	    } catch (NumberFormatException e) {
    	        System.err.println("Argument" + args[1]);
    	    }
    	    try {
    	        yArg = Integer.parseInt(args[2]);
    	    } catch (NumberFormatException e) {
    	        System.err.println("Argument" + args[2]);
    	    }
    	    try {
    	        sig = Double.parseDouble(args[3]);
    	    } catch (NumberFormatException e) {
    	        System.err.println("Argument" + args[3] + " is not an double. Sigma set to 1 by default");
    	    }
    	    try {
    	        sigStop = Double.parseDouble(args[4]);
    	    } catch (NumberFormatException e) {
    	        System.err.println("Argument" + args[4] + " is not an integer. sigma stop set to 0.4 by default");
    	    }
    	    try {
    	        iter = Integer.parseInt(args[5]);
    	    } catch (NumberFormatException e) {
    	        System.err.println("Argument" + args[5] + " is not an integer. Iterations set to 5000 by default");
    	    }
    	    
    	    try {
    	        n = Integer.parseInt(args[6]);
    	    } catch (NumberFormatException e) {
    	        System.err.println("Argument" + args[6] + " is not an integer. n set to 10 by default");
    	    }
    	    
    	    boolean cosine = Integer.parseInt(args[7])==1;
    	    String outPrefix = args[8];
    	    String mat = args[9];
    	    int availThreads = Integer.parseInt(args[10]);
    	    boolean countIntra = Integer.parseInt(args[11])==1;
    	    
    	    SOMTrainer trainer = new SOMTrainer(n, xArg,yArg, sig, sigStop, iter, cosine, outPrefix, mat, availThreads, countIntra);
    	    trainer.train();
    	    
    	}
    	
	}
	 
}


