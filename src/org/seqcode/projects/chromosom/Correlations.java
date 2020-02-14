package org.seqcode.projects.chromosom;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Scanner;



public class Correlations
{
	public DrawHex d;
	public ArrayList<String> fold;
	public String mapper;
	public double gini, pVal;
	public int nodes,xs,ys,xNodes,yNodes, maxBins, minBins, colorNum, winW, winH, pValnVal, locCol, weightCol, binSize;
	public ArrayList<DataPoint> bins;
	public ArrayList<MiniNode> nodeList;
	public ArrayList<String> names;
	public ArrayList<File> searchers;
	public MiniSystem nodeSystem;
	public boolean weighting, showPVal, addToOldFile, pics, equalWeight;
	public double[][] g;
	public double[] sep;
	
	public Correlations(String map, ArrayList<String> fillle)
	{
		equalWeight = true;
		mapper = map;
		fold = fillle;
		locCol = 0;
		weightCol = 1;
		pValnVal= 1000;
		showPVal = true;
		gini=0;
		weighting = true;
		String mapp = map;
		yNodes=0;
		xNodes=0;
		searchers = new ArrayList<File>();
		names = new ArrayList<String>();
		mapReader(mapp);
		pVal = 0;
		binSize = bins.get(0).maxLocus - bins.get(0).minLocus;
		for(String fs: fold)
		{
			File folder = new File(fs);
			File[] listOfFiles = folder.listFiles();
			for (File file : listOfFiles) 
			{
			    if (file.isFile() && file.getName().endsWith(".bed")) 
			    {
			    	searchers.add(file);
			    	names.add(file.getName());
			    }
			}
		}
		nodes = yNodes*xNodes;
		
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
					System.out.println(searchers.get(i).getName()+" vs. "+searchers.get(j).getName());
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
		double totaldd = 0;
		double totalnn = 0;
		for(int i = 0; i<nodal.length; i++)
		{
			totalnn += nodal[i];
		}
		for(int i = 0; i<datal.length; i++)
		{
			totaldd+=datal[i];
		}
		for(int i = 0; i<nodal.length; i++) {nodal[i]/=totalnn;}
		for(int i = 0; i<datal.length; i++) {datal[1]/=totaldd;}
		
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

	    nodal = nnn;
	    datal = ddd;
	    return rr;
	}
	public void search(File file)
	{
		for(int i = 0; i < nodeSystem.size(); i++)
		{
			MiniNode mini = nodeSystem.get(i);
			mini.counting.clear();
			mini.weight = 0;
		}

		BEDReader reader;
		try {
			reader = new BEDReader(file);
		
			while(reader.hasNext())
			{
				Region reg = reader.next();
				double weight;
				if(equalWeight)
					weight = 1.0;
				else
					weight = reg.score / 1000; //Note this use of the BED score definition in the manual. 
					
				ArrayList<Integer> inds = new ArrayList<Integer>();
				double count = 0;
				
				//Iterate over bins (datapoints in SOM)
				for(int i = 0; i< bins.size(); i++){
					//Overlap assessment 
					if(bins.get(i).chrome.equals(reg.chr) && 
							((bins.get(i).minLocus <= reg.start && bins.get(i).maxLocus >= reg.start) ||
							(reg.start <= bins.get(i).minLocus && reg.end >= bins.get(i).minLocus )) )
							{
								count++;
								inds.add(i);
								
								//Check for contained relationship to shortcut
								if(bins.get(i).chrome.equals(reg.chr) && 
										(bins.get(i).minLocus <= reg.start && bins.get(i).maxLocus >= reg.end) ) {
									break;
								}
							}
				}
				
				if(!equalWeight)
					weight/=count;
				for(Integer i: inds){
					bins.get(i).myMini.counting.add(bins.get(i));
					bins.get(i).myMini.weight += weight;
				}
			}
		} catch (IOException e) {
			e.printStackTrace();
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
			@SuppressWarnings("resource")
			Scanner in = new Scanner(new FileReader(file));
			in.next();
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
			@SuppressWarnings("resource")
			Scanner in = new Scanner(new FileReader(f));
			String sizer = in.next();
			int xo = Integer.parseInt(sizer.substring(0,sizer.indexOf("x")));
			int yo = Integer.parseInt(sizer.substring(sizer.indexOf("x")+1,sizer.length()));
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
		}
	}
}

