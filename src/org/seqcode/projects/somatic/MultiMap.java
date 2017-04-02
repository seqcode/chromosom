package org.seqcode.projects.somatic;
import java.awt.BorderLayout;
import java.awt.event.*;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;

import javax.swing.*;
import javax.swing.filechooser.FileNameExtensionFilter;



public class MultiMap extends JFrame
{
    static MultiMap window2;
    int winW, winH;
	public DrawHex d;
	public UseMap u;
	private static final long serialVersionUID = 1L;
	public ArrayList<String> s;
	public String find;
	public MultiMap(int x, int y, String p)
	  {
		 
		//title the window
	    super(p);
	    
	    BorderLayout gui = new BorderLayout();
	    setLayout(gui);
	    winW = x; winH = y;
	    
	    s = new ArrayList<String>();
	  }
	
	
	
	/** will this being static cause problems? */
	public static ThreadMap qualityControl(ArrayList<ThreadMap> op, boolean similarity)
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
//		public void run()
//		{
//			for(int i = 0; i < numberOfSoms; i++)
//			{
//    	    	ThreadMap o = new ThreadMap(xArg,yArg,sigmaArg, lr, iter, Integer.parseInt(args[7]), Integer.parseInt(args[8]), Integer.parseInt(args[9]), args[10], args[11], threadsFor));
//   		 		o.go();
//   		 		t.add(o);
//    	    }
//		}
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
//		public void run()
//		{
//			for(int i = 0; i < numberOfSoms; i++)
//			{
//    	    	ThreadMap o = new ThreadMap(xArg,yArg,sigmaArg, lr, iter, Integer.parseInt(args[7]), Integer.parseInt(args[8]), Integer.parseInt(args[9]), args[10], args[11], threadsFor));
//   		 		o.go();
//   		 		t.add(o);
//    	    }
//		}
	}
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
    	    ArrayList<ThreadMap> op = new ArrayList<ThreadMap>();
    	    for(int i = 0; i< n; i++)
    	    {
        	    System.out.println("working");
    	    	ThreadMap o = new ThreadMap(xArg,yArg, sig, sigStop, iter, Integer.parseInt(args[7]) , args[8], args[9], Integer.parseInt(args[10]));
   		 		o.go();
   		 		op.add(o);
   		 		System.out.println(((double)(i+1))/((double)n));
    	    }
    	    boolean boo = Integer.parseInt(args[7]) == 1;
    	    qualityControl(op, boo).writeFile(n);
    	}
    	//train 10 10 15 LR iterations n 1 1 1 "Name" "File Name" 4 2
    	if(args[0].equalsIgnoreCase(("Test")))
    	{
    		ArrayList<ThreadMap> op = new ArrayList<ThreadMap>();
    	    int j = 20;
    	    for(int lr = 1; lr <=3; lr++)
    	    {   
	    	    while (j>4)
	    	    {
		    		for(int i = 1; i< 5; i++)
		    	    {
		    			//"Test" "Name" "FileName" Threads
		    	    	//ThreadMap o = new ThreadMap(50,50, j, lr, 500, 1, 1, 1, args[1], args[2], Integer.parseInt(args[3]), i);
		   		 		//o.go();
		   		 		//op.add(o);
		   		 		qualityControl(op, true).writeFile(1);
		    	    }
		    	    j-=5;
		    	}
    		}
    	}
	}
	 
}


