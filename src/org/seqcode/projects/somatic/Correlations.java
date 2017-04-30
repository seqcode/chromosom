package org.seqcode.projects.somatic;
import java.awt.BorderLayout;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Scanner;

import org.apache.commons.math3.distribution.NormalDistribution;



public class Correlations
{
	public DrawHex d;
	public ArrayList<String> fold;
	public String mapper, der;
	public double gini, pVal;
	public int nodes,xs,ys,xNodes,yNodes, maxBins, minBins, colorNum, winW, winH, pValnVal, locCol, weightCol, binSize;
	public ArrayList<Node> coords;
	public ArrayList<DataPoint> bins;
	public ArrayList<MiniNode> nodeList;
	public ArrayList<String> searchers,names;
	MiniSystem nodeSystem;
	public boolean weighting, showPVal, addToOldFile, pics, equalWeight;
	public double[][] g;
	public double[] sep;
	
	public Correlations(String map, ArrayList<String> fillle)
	{
		equalWeight = true;
		mapper = map;
		fold = fillle;
		der =System.getProperty("user.dir")+"/"+"Lorenz_Area_"+fold;
		locCol = 0;
		weightCol = 1;
		pValnVal= 1000;
		showPVal = true;
		gini=0;
		weighting = true;
		String mapp = System.getProperty("user.dir")+"/"+map;
		yNodes=0;
		xNodes=0;
		searchers = new ArrayList<String>();
		names = new ArrayList<String>();
		mapReader(mapp);
		pVal = 0;
		binSize = bins.get(0).maxLocus - bins.get(0).minLocus;
		for(String fs: fold)
		{
			File folder = new File(System.getProperty("user.dir")+"/"+fs);
			File[] listOfFiles = folder.listFiles();
			for (File file : listOfFiles) 
			{
			    if (file.isFile() && (file.getName().indexOf(".txt") != -1 || file.getName().indexOf(".domains")!=-1 || file.getName().indexOf(".peaks")!=-1 || file.getName().indexOf(".narrowpeaks")!=-1)|| file.getName().indexOf(".bed")!=-1) 
			    {
			    	searchers.add(file.getPath());
			    	names.add(file.getName());
			    }
			}
		}
		nodes = yNodes*xNodes;
		for(int i = 0; i < searchers.size(); i++)
		{
			System.out.println(searchers.get(i));
		}
		heatMapping();
	}
	public void searchSystem()
	{
		double[][] correlations = new double[searchers.size()][searchers.size()];
		for(int i = 0; i < searchers.size(); i++)
		{
			double[] array1 = new double[nodeList.size()];
			search(searchers.get(i));
			for(int k = 0; k<nodeList.size(); k++)
			{
				array1[k] = nodeList.get(k).counting.size();
			}
			for(int j = i; j < searchers.size(); j++)
			{
				if(i == j)
					correlations[i][j] = 1;
				else
				{
					double[] array2 = new double[nodeList.size()];
					search(searchers.get(j));
					for(int k = 0; k<nodeList.size(); k++)
					{
						array2[k] = nodeList.get(k).counting.size();
					}
					double corr = pearson(array1,array2);
					correlations[i][j] = corr;
					correlations[j][i] = corr;
				}
			}
		}
		writeFile(correlations);
		System.out.println("done");
	}
	public double pearson(double[] nodal, double[] datal)
	{
		double[] nnn = nodal;
		double[] ddd = datal;
		double maxdd = 0;
		double maxnn = 0;
		for(int i = 0; i<nodal.length; i++)
		{
			if(nodal[i]>maxnn)
				maxnn = nodal[i];
		}
		for(int i = 0; i<datal.length; i++)
		{
			if(datal[i]>maxdd)
				maxdd=datal[i];
		}
		for(int i = 0; i<nodal.length; i++) {nodal[i]/=maxnn;}
		for(int i = 0; i<datal.length; i++) {datal[1]/=maxdd;}
		
		double sx = 0.0;
	    double sy = 0.0;
	    double sxx = 0.0;
	    double syy = 0.0;
	    double sxy = 0.0;

	    int n = nodal.length;
	    for(int i = 0; i < n; i++) 
	    {
	      double x = nodal[i];
	      double y = datal[i];

	      sx += x;
	      sy += y;
	      sxx += x * x;
	      syy += y * y;
	      sxy += x * y;
	    }

	    // covariation
	    double cov = sxy / n - sx * sy / n / n;
	    // standard error of x
	    double sigmax = Math.sqrt(sxx / n -  sx * sx / n / n);
	    // standard error of y
	    double sigmay = Math.sqrt(syy / n -  sy * sy / n / n);

	    // correlation is just a normalized covariation
	    
	    double rr = cov / sigmax / sigmay;
	    if(rr > 1 && rr< 1.0000001)
	    	rr=1;
	    if(rr>1.0000001)
	    	System.out.println("whoops");
	    nodal = nnn;
	    datal = ddd;
	    return rr;
	}
	public void search(String file)
	{
		for(int i = 0; i < nodeSystem.size(); i++)
		{
			MiniNode mini = nodeSystem.get(i);
			mini.counting.clear();
			mini.weight = 0;
		}
		int col = 0 ;
		int www = 0;
		ArrayList<String> strings = inputRead(file);
		for(String whole: strings)
		{
			int chr = 0;
			int locus1 = 0; int locus2 = 0;
			double weight = 0;
			String[] wholes = whole.split("\t");
			//System.out.println(whole);
			if(whole.contains("\t") && !whole.contains(":") && whole.contains("chr")&& !wholes[col].contains("-") && !whole.contains("_") )
			{
				if(wholes[col].substring(wholes[col].indexOf("chr")+3).equalsIgnoreCase("X")||wholes[col].substring(wholes[col].indexOf("chr")+3).equalsIgnoreCase("Y"))
					chr = 23;
				else 
					chr  = Integer.parseInt(wholes[col].substring(whole.indexOf("chr")+3));
				String loc1 = wholes[1];
				locus1 = Integer.parseInt(loc1);
				locus2 = Integer.parseInt(wholes[2]);
				if(equalWeight)
					weight = 1.0;
				else
					weight = Double.parseDouble(wholes[www]);     											/** This needs to be taken as an argument somehow, not hard coded*/
				int locus = locus1;
				ArrayList<Integer> inds = new ArrayList<Integer>();
				int count = 0;
				while(locus<locus2)
				{
					for(int i = 0; i< bins.size(); i++)
					{
						if(bins.get(i).chrome == chr && bins.get(i).minLocus <= locus && bins.get(i).maxLocus>=locus)
						{
							count++;
							inds.add(i);
						}
					}
					locus+=binSize;
				}
				weight/=count;
				for(Integer i: inds)
				{
					bins.get(i).myMini.counting.add(bins.get(i));
					bins.get(i).myMini.weight += weight;
				}
			}
			else if(whole.contains("\t") && whole.contains(":")&&whole.contains("chr")&& !wholes[col].contains("-") && !whole.contains("_")&& !whole.substring(whole.indexOf("chr")+3,whole.indexOf(":")).equalsIgnoreCase("M"))
			{
				weighting = true;
				if(wholes[col].substring(whole.indexOf("chr")+3,whole.indexOf(":")).equalsIgnoreCase("X")||wholes[col].substring(whole.indexOf("chr")+3,whole.indexOf(":")).equalsIgnoreCase("Y"))
						chr = 23;
				else 
					chr  = Integer.parseInt(wholes[col].substring(whole.indexOf("chr")+3,whole.indexOf(":")));
				locus1 = Integer.parseInt(wholes[col].substring(whole.indexOf(":")+1));
				//locus2 = Integer.parseInt(whole.substring(whole.indexOf("-")+1,whole.indexOf("\t")));
				int locus = locus1;//(locus1 +locus2)/2;
				if(equalWeight)
					weight = 1.0;
				else
					weight = Double.parseDouble(wholes[www]);     											/** This needs to be taken as an argument somehow, not hard coded*/
				for(int i = 0; i< bins.size(); i++)
				{
					if(bins.get(i).chrome == chr && bins.get(i).minLocus <= locus && bins.get(i).maxLocus>=locus)
					{
						bins.get(i).myMini.counting.add(bins.get(i));
						bins.get(i).myMini.weight += weight;
					}
				}
			}
			else if (whole.contains("\t") && whole.contains(":")&&whole.contains("chr")&& wholes[col].contains("-") && !whole.contains("_") && !whole.substring(whole.indexOf("chr")+3,whole.indexOf(":")).equalsIgnoreCase("M"))
			{
				weighting = true;
				if(wholes[col].substring(whole.indexOf("chr")+3,whole.indexOf(":")).equalsIgnoreCase("X")||wholes[col].substring(whole.indexOf("chr")+3,whole.indexOf(":")).equalsIgnoreCase("Y"))
						chr = 23;
				else 
					chr  = Integer.parseInt(wholes[col].substring(wholes[col].indexOf("chr")+3,wholes[col].indexOf(":")));
				
				String loc1 = wholes[col].substring(wholes[col].indexOf(":")+1, wholes[col].indexOf("-"));
				locus1 = Integer.parseInt(loc1);
				locus2 = Integer.parseInt(wholes[col].substring(wholes[col].indexOf("-")+1));
				if(equalWeight)
					weight = 1.0;
				else
					weight = Double.parseDouble(wholes[www]);     											/** This needs to be taken as an argument somehow, not hard coded*/
				int locus = locus1;
				ArrayList<Integer> inds = new ArrayList<Integer>();
				int count = 0;
				while(locus<locus2)
				{
					for(int i = 0; i< bins.size(); i++)
					{
						if(bins.get(i).chrome == chr && bins.get(i).minLocus <= locus && bins.get(i).maxLocus>=locus)
						{
							count++;
							inds.add(i);
						}
					}
					locus+=binSize;
				}
				weight/=count;
				for(Integer i: inds)
				{
					bins.get(i).myMini.counting.add(bins.get(i));
					bins.get(i).myMini.weight += weight;
				}
			}
		}
	}
	public String coClustering(String file1, String file2)
	{
		return "";
	}
	public ArrayList<String> inputRead(String file)
	{
		ArrayList<String> StringMat = new ArrayList<String>();
		try 
		{
			Scanner in = new Scanner(new FileReader(file));
			in.next();
			//System.out.println(xo + "  x  "+ yo);
			in.next();
			while(in.hasNextLine())
			{
				StringMat.add(in.nextLine());
			}
		} catch (FileNotFoundException e) {e.printStackTrace();}
		return StringMat;

	}
	
	public void mapReader(String f)
	{
		ArrayList<String> StringMat = new ArrayList<String>();
		try 
		{
			Scanner in = new Scanner(new FileReader(f));
			String sizer = in.next();
			int xo = Integer.parseInt(sizer.substring(0,sizer.indexOf("x")));
			int yo = Integer.parseInt(sizer.substring(sizer.indexOf("x")+1,sizer.length()));
			//System.out.println(xo + "  x  "+ yo);
			in.next();
			while(in.hasNext())
			{
				StringMat.add(in.next());
			}
			createNoder(StringMat, xo, yo);
		} catch (FileNotFoundException e) {e.printStackTrace();}
	}
	
	public void createNoder(ArrayList<String> strings, int x, int y)
	{
		nodeList = new ArrayList<MiniNode>(x*y);
		bins = new ArrayList<DataPoint>();
		MiniNode current = null;
		xNodes = x;
		yNodes = y;
		for(String whole: strings)
		{
			if(whole.contains("["))
			{
				current = new MiniNode(whole);
				nodeList.add(current);
			}
			else
				current.addDP(whole);
		}
		for(int i = 0; i<nodeList.size(); i++)
		{
			for(int j = 0; j < nodeList.get(i).bins.size(); j++)
			{
				bins.add(nodeList.get(i).bins.get(j));
			}
		}
		//System.out.println(bins.size());
		nodeSystem = new MiniSystem(nodeList,xNodes,yNodes);
	}
	
	public void heatMapping()
	{
		minBins = (int) (nodeList.get(0).bins.size());
		maxBins = (int) (nodeList.get(0).bins.size());
		for(int i = 0; i<nodeList.size(); i++)
		{
			if(nodeList.get(i).bins.size()>maxBins)
				maxBins = (int) (nodeList.get(i).bins.size());
			if(nodeList.get(i).bins.size()<minBins)
				minBins = (int) (nodeList.get(i).bins.size());
		}
		if(maxBins == minBins) maxBins++;
	}
	public void writeFile(double[][] cors)
	{
		BufferedWriter b;
		try 
		{
			FileWriter ff;
			ff = new FileWriter("cors.txt",true);
			b = new BufferedWriter(ff);
			PrintWriter printer = new PrintWriter(b);
			printer.print("\t");
			for(String s: names)
			{
				printer.print(s+"\t");
			}
			printer.println();
			for(int i = 0; i<cors.length; i++)
			{
				printer.print(names.get(i)+"\t");
				for(int j = 0; j<cors[i].length;j++)
				{
					printer.print(cors[i][j] + "\t");
				}
				printer.println();
			}
			printer.close();
		}
		catch (IOException e) 
		{
			e.printStackTrace();
			System.out.print("no way");
		}
	}
}

