package org.seqcode.projects.somatic;
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



public class UseMap
{
	public DrawHex d;
	public String fold, mapper, der;
	public double gini, pVal;
	public int nodes,xs,ys,xNodes,yNodes, maxBins, minBins, colorNum, winW, winH, pValnVal, locCol, weightCol, binSize;
	public ArrayList<DataPoint> bins;
	public ArrayList<MiniNode> nodeList;
	public ArrayList<String> searchers;
	public MiniSystem nodeSystem;
	public boolean weighting, showPVal, addToOldFile, pics, equalWeight;
	public double[][] g;
	public double[] sep;
	
	public UseMap(String map, String fillle, int locCols, int weightCols, boolean pic, boolean weightin)
	{
		pics = pic;
		mapper = map;
		fold = fillle;
		der =System.getProperty("user.dir")+"/"+"Lorenz_Area_"+fold;
		locCol = locCols;
		weightCol = weightCols;
		pValnVal= 1000;
		showPVal = true;
		gini=0;
		weighting = true;
		String mapp = System.getProperty("user.dir")+"/"+map;
		equalWeight = weightin;
		yNodes=0;
		xNodes=0;
		searchers = new ArrayList<String>();
		mapReader(mapp);
		pVal = 0;
		binSize = bins.get(0).maxLocus - bins.get(0).minLocus;
		File folder = new File(System.getProperty("user.dir")+"/"+fillle);
		File[] listOfFiles = folder.listFiles();

		for (File file : listOfFiles) 
		{
		    if (file.isFile() && (file.getName().indexOf(".txt") != -1 || file.getName().indexOf(".domains")!=-1 || file.getName().indexOf(".peaks")!=-1 || file.getName().indexOf(".narrowpeaks")!=-1)|| file.getName().indexOf(".bed")!=-1) 
		    {
		    	searchers.add(file.getName());
		    }
		}
		nodes = yNodes*xNodes;
		for(int i = 0; i < searchers.size(); i++)
		{
			System.out.println(searchers.get(i));
		}
		heatMapping();
	}
	public void searchSystemDouble()
	{
		if(pics)
		{
			new File(der).mkdir();
			d = new DrawHex(mapper);
			d.nodeBuild(700,700);
	    	d.heatMapping();
	    	d.countingDPS(null);
	    	d.saveImg("whole_map", der,  700, 700);
		}
		ArrayList<String> writer = new ArrayList<String>();
		writer.add("File \t gini \t giniW/WeightedBaseline(n="+pValnVal+")\tzScore\n");
		for (int i = 0; i < searchers.size()-1; i++)
		{
			for(int j = i; j < searchers.size(); j++)
			{
				g = new double[nodes+2][pValnVal+1];
				gini = 0;
				sep = new double[pValnVal+1];
				multiSearch(searchers.get(i), searchers.get(j));
				writer.add(searchers.get(i) + "_" + searchers.get(j) + coClustering(searchers.get(i), searchers.get(j)));
			}
		}
		writeFile(writer);
		System.out.println("done");
	}
	public void searchSystem()
	{
		if(pics)
		{
			new File(der).mkdir();
			d = new DrawHex(mapper);
			d.nodeBuild(700,700);
	    	//d.colors();
	    	d.heatMapping();
	    	d.countingDPS(null);
	    	d.saveImg("whole_map", der,  700, 700);
		}
		ArrayList<String> writer = new ArrayList<String>();
		writer.add("File \t gini \t giniW/WeightedBaseline(n="+pValnVal+")\tzScore\n");
		for (String s:searchers)
		{
			g = new double[nodes+2][pValnVal+1];
			gini = 0;
			sep = new double[pValnVal+1];
			search(s);
			writer.add(s + "\t" + g[nodes+1][0] +"\t" + g[nodes][0] + "\t" + g[nodes+1][1]);
		}
		writeFile(writer);
		System.out.println("done");
	}
	public void multiSearch(String file1, String file2)
	{
		if(pics)
		{
			//d.search(fold+"/"+file, locCol, weightCol);
			d.multiSearch(fold+"/"+file1, fold+"/"+file2);
			d.repaint();
			d.saveImg(file1+"_"+file2, der, 700, 700);
		}
	}
	public void search(String file)
	{
		if(pics)
		{
			d.equalWeight = (weightCol == -1);
			d.search(fold+"/"+file, locCol, weightCol);
			d.repaint();
			d.saveImg(file, der, 700, 700);
		}
		for(int i = 0; i < nodeSystem.size(); i++)
		{
			MiniNode mini = nodeSystem.get(i);
			mini.counting.clear();
			mini.weight = 0;
		}
		ArrayList<String> strings = inputRead(file);
		for(String whole: strings)
		{
			String[] wholes = whole.split("\t");
			String chr;
			int locus1 = 0; int locus2 = 0;
			double weight = 0;
			if(whole.contains("\t") && !whole.contains(":") && whole.contains("chr")&& !wholes[locCol].contains("-") && !whole.contains("_") )
			{
				chr=wholes[0];
				String loc1 = wholes[1];
				locus1 = Integer.parseInt(loc1);
				locus2 = Integer.parseInt(wholes[2]);
				if(equalWeight)
					weight = 1.0;
				else
					weight = Double.parseDouble(wholes[weightCol]);     											/** This needs to be taken as an argument somehow, not hard coded*/
				int locus = locus1;
				ArrayList<Integer> inds = new ArrayList<Integer>();
				int count = 0;
				while(locus<locus2)
				{
					for(int i = 0; i< bins.size(); i++)
					{
						if(bins.get(i).chrome.equals(chr) && bins.get(i).minLocus <= locus && bins.get(i).maxLocus>=locus)
						{
							
							count++;
							inds.add(i);
						}
					}
					locus+=binSize;
				}
				if(showPVal)
				{
					for(int p = 1; p<g[0].length; p++)
					{
						int ree = (int) (Math.random()*bins.size());
						ree = nodeList.indexOf(bins.get(ree).myMini);
						g[ree][p] += weight;
					}
				}
				weight/=count;
				for(Integer i: inds)
				{
					bins.get(i).myMini.counting.add(bins.get(i));
					bins.get(i).myMini.weight += weight;
					g[nodeList.indexOf(bins.get(i).myMini)][0] += weight;
					//System.out.println(weight);
				}
			}
			else if(whole.contains("\t") && whole.contains(":")&&whole.contains("chr")&& !wholes[locCol].contains("-") && !whole.contains("_"))
			{
				weighting = true;
				chr = wholes[locCol].split(":")[0];
				locus1 = Integer.parseInt(wholes[locCol].substring(whole.indexOf(":")+1));
				//locus2 = Integer.parseInt(whole.substring(whole.indexOf("-")+1,whole.indexOf("\t")));
				int locus = locus1;//(locus1 +locus2)/2;
				if(equalWeight)
					weight = 1.0;
				else
					weight = Double.parseDouble(wholes[weightCol]);     											/** This needs to be taken as an argument somehow, not hard coded*/
				for(int i = 0; i< bins.size(); i++)
				{
					if(bins.get(i).chrome.equals(chr) && bins.get(i).minLocus <= locus && bins.get(i).maxLocus >=locus)
					{
						//System.out.println(bins.get(i).name +" "+ bins.get(i).minLocus + "  "+ bins.get(i).maxLocus + "     "+ locus + "   " + whole);
						bins.get(i).myMini.counting.add(bins.get(i));
						bins.get(i).myMini.weight += weight;
						g[nodeList.indexOf(bins.get(i).myMini)][0] += weight;
						//System.out.println(weight);
						if(showPVal)
						{
							for(int p = 1; p<g[0].length; p++)
							{
								int ree = (int) (Math.random()*bins.size());
								ree = nodeList.indexOf(bins.get(ree).myMini);
								g[ree][p] += weight;
							}
						}
					}
				}
			}
			else if (whole.contains("\t") && whole.contains(":")&&whole.contains("chr")&& wholes[locCol].contains("-") && !whole.contains("_") && !whole.substring(whole.indexOf("chr")+3,whole.indexOf(":")).equalsIgnoreCase("M"))
			{
				weighting = true;
				chr = wholes[locCol].split(":")[0];
				
				String loc1 = wholes[locCol].substring(wholes[locCol].indexOf(":")+1, wholes[locCol].indexOf("-"));
				locus1 = Integer.parseInt(loc1);
				locus2 = Integer.parseInt(wholes[locCol].substring(wholes[locCol].indexOf("-")+1));
				if(equalWeight)
					weight = 1.0;
				else
					weight = Double.parseDouble(wholes[weightCol]);     											/** This needs to be taken as an argument somehow, not hard coded*/
				int locus = locus1;
				ArrayList<Integer> inds = new ArrayList<Integer>();
				int count = 0;
				while(locus<locus2)
				{
					for(int i = 0; i< bins.size(); i++)
					{
						if(bins.get(i).chrome.equals(chr) && bins.get(i).minLocus <= locus && bins.get(i).maxLocus>=locus)
						{
							
							count++;
							inds.add(i);
						}
					}
					locus+=binSize;
				}
				if(showPVal)
				{
					for(int p = 1; p<g[0].length; p++)
					{
						int ree = (int) (Math.random()*bins.size());
						ree = nodeList.indexOf(bins.get(ree).myMini);
						g[ree][p] += weight;
					}
				}
				weight/=count;
				for(Integer i: inds)
				{
					bins.get(i).myMini.counting.add(bins.get(i));
					bins.get(i).myMini.weight += weight;
					g[nodeList.indexOf(bins.get(i).myMini)][0] += weight;
					//System.out.println(weight);
				}
			}
		}
		giniWithWeightedBaseline();
		gini();
		pVal();
		/**degreesOfSep();*/
	}
	public String coClustering(String file1, String file2)
	{
		return "";
	}
	public void pVal()
	{
		double[] pValFind = new double[g[nodes].length-1];
		double mean = 0;
		for(int i = 0; i< pValFind.length; i++)
		{
			mean+= g[nodes][i+1];
			pValFind[i] = g[nodes][i+1];
		}
		mean/=pValnVal;
		double std = 0;
		for(int i = 0; i< pValFind.length; i++)
		{
			std += (g[nodes][i+1] - mean)*(g[nodes][i+1] - mean);
		}
		if(std!=0)
		{
			std/=pValnVal;
			std = Math.sqrt(std);
			g[nodes+1][1] = (g[nodes][0]-mean)/std;
		}
	}
	public void degreesOfSep()
	{
		double[][] dist = new double[nodes][nodes];
		for (int i = 0; i<nodes; i++)
		{
			int xCoordA = i%xNodes;
			int yCoordA = i/yNodes;
			for (int j = 0; j<nodes; j++)
			{
				int xCoordB = j%xNodes;
				int yCoordB = j/yNodes;
				
				double ddd = hexDist(xCoordA, yCoordA, xCoordB, yCoordB);
				dist[i][j] = ddd;
			}
		}
		
		for(int i = 0; i<g[0].length; i++)
		{
			double count = 0;
			for(int j = 0; j<nodes; j++)
			{
				for(int k = 0; k<nodes; k++)
				{
					sep[i] += dist[j][k]*g[j][i]*g[k][i];
					count+= g[j][i]*g[k][i];
				}
			}
			sep[i] /= count;
		}
	}
	public double degreesP()
	{
		if(showPVal)
		{
			double slit = 0;
			double avg = 0;
			for(int i = 1; i < sep.length; i++)
			{
				avg += sep[i];
			}
			avg /= pValnVal;
			double std = 0;
			for(int i = 1; i < pValnVal+1; i++)
			{
				std += (sep[i]-avg)*(sep[i]-avg);
			}
			if(std==0)
				std = 0.00001;
			std /= pValnVal;
			std = Math.sqrt(std);
			NormalDistribution norm = new NormalDistribution(avg,std);
			slit = norm.cumulativeProbability(sep[0]);
			return slit;
		}
		else return 0;
	}
	public double hexDist(int x1, int y1, int x2, int y2)
	{
		if(x1==x2 && y2 == y1)
			return 0;
		else
		{
			int hm = 0;
			int vm = 0;
			int cm = 0;
			boolean right = x1-x2<0;
			hm = Math.abs(x1-x2);
			vm = Math.abs(y1-y2);
			if(xNodes-x1<xNodes-x2)
				hm = Math.min(xNodes-x1 + x2,hm);
			else
				hm = Math.min(xNodes-x2 +x1, hm);
			
			if(yNodes-y1<yNodes-y2)
				vm = Math.min(yNodes-y1 + y2,vm);
			else
				vm = Math.min(yNodes-y2 +y1, vm);
			
			if(vm %2 !=0)
			{
				if(!right && y1%2 == 0)
				{
					cm = Math.min(hm, (int)(((double)vm)/2 +.5));
					//System.out.println("(" + x1+ ", " + y1+ ")   ->   ("  +x2 + ", " +  y2 +")  = " + (hm-cm+vm-cm+cm));
				}
				else if(right && y1%2 == 1)
				{
					cm = Math.min(hm, (int)(((double)vm)/2 +.5));
					//System.out.println("(" + x1+ ", " + y1+ ")   ->   ("  +x2 + ", " +  y2 +")  = " + (hm-cm+vm-cm+cm));
				}
			}
			else
				cm = Math.min(hm, (int)(((double)vm)/2));
			return hm+vm-cm; 
		}
	}
	public void giniWithWeightedBaseline()
	{
		// ~equality baseline
		double eArea = 0;
		for(int num = 1; num<=pValnVal; num++)
		{
			int maxWeight = 0;
			for(int i = 0; i<nodeList.size(); i++)
			{
				if(g[i][num]>maxWeight)
					maxWeight = (int) g[i][num];
			}
			double[] equal = new double[maxWeight+1];
			for(int i=0; i<nodeList.size(); i++)
			{
				int o = (int) g[i][num];
				equal[o]++;
			}
			
			int eQTotalWealth = 0;
			for(int i=0; i<equal.length; i++)
			{

				if(num==1 && equal[i]>0)
					System.out.println(i+"\t"+equal[i]);
				eQTotalWealth += i*equal[i];
			}
			
			double[] eQCumPop = new double[maxWeight+1];
			double[] eQCumWealth = new double[maxWeight+1];
			
			eQCumPop[0] = equal[0]/(double)nodeList.size();
			eQCumWealth[0] = 0;
			for(int i = 1; i<equal.length; i++)
			{
				eQCumPop[i] = equal[i]/(double)nodeList.size() + eQCumPop[i-1];
				eQCumWealth[i] = (equal[i]*i)/eQTotalWealth + eQCumWealth[i-1];
			}		
			//Equaliy area
			for(int i=0; i<equal.length-1; i++)
			{
				g[nodes][num]+=(((eQCumPop[i+1]-eQCumPop[i])*eQCumWealth[i])+((eQCumPop[i+1]-eQCumPop[i])*eQCumWealth[i+1]))/2;
				eArea += (((eQCumPop[i+1]-eQCumPop[i])*eQCumWealth[i])+((eQCumPop[i+1]-eQCumPop[i])*eQCumWealth[i+1]))/2;
			}
		}
		eArea /= pValnVal;
		//System.out.println(eArea);
		
		
		//inequality
		double maxWeight = 0;
		for(int j =0; j<g.length; j++)
		{
			if(g[j][0] > maxWeight)
				maxWeight=g[j][0];
		}
		
		/**
		//casting weights to ints may put all to zero in some data types?
		**/
		
		int[] unequal = new int[(int)maxWeight+1];
		for(int i1=0; i1<g.length-1; i1++)
		{
			int o = (int)g[i1][0];
			unequal[o] ++;
		}
		
		double uETotalWealth = 0;
		
		for(int i1=0; i1<unequal.length; i1++)
		{
			uETotalWealth += i1*unequal[i1];
		}

		double[] uECumPop = new double[unequal.length];
		double[] uECumWealth = new double[unequal.length];
		
		uECumPop[0] = unequal[0]/(double)nodeList.size();
		uECumWealth[0] = 0;
		
		for(int i1 = 1; i1<unequal.length; i1++)
		{
			uECumPop[i1] = unequal[i1]/(double)nodeList.size() + uECumPop[i1-1];
			
			uECumWealth[i1] = (unequal[i1]*i1)/uETotalWealth + uECumWealth[i1-1];
		}
		
		//Lorenz area
		double lArea = 0;
		for(int i1=0; i1<unequal.length-1; i1++)
		{
			lArea += (((uECumPop[i1+1]-uECumPop[i1])*uECumWealth[i1])+((uECumPop[i1+1]-uECumPop[i1])*uECumWealth[i1+1]))/2;
		}
		//difference in area
		double dArea = eArea - lArea;
		gini = dArea / (lArea + dArea);
		g[nodes][0] = gini;
		gini=0;
	}
	public void gini()
	{
		// ~equality baseline
		int mdp = nodeList.get(0).bins.size();
		for(int i = 0; i<nodeList.size(); i++)
		{
			if(nodeList.get(i).bins.size()>mdp)
				mdp = (int) (nodeList.get(i).bins.size());
		}
		double[] equal = new double[mdp+1];
		for(int i=0; i<nodeList.size(); i++)
		{
			int o = nodeList.get(i).bins.size();
			equal[o]++;
		}
		
		int eQTotalWealth = 0;
		for(int i=0; i<equal.length; i++)
		{
			eQTotalWealth += i*equal[i];
		}
		//System.out.println("use equal "+ eQTotalWealth);
		
		
		double[] eQCumPop = new double[mdp+1];
		double[] eQCumWealth = new double[mdp+1];
		eQCumPop[0] = equal[0]/nodeList.size();
		eQCumWealth[0] = 0;
		for(int i = 1; i<equal.length; i++)
		{
			eQCumPop[i] = equal[i]/nodeList.size() + eQCumPop[i-1];
			eQCumWealth[i] = (equal[i]*i)/eQTotalWealth + eQCumWealth[i-1];
		}		
		//Equaliy area
		double eArea = 0;
		for(int i=0; i<equal.length-1; i++)
		{
			eArea += (((eQCumPop[i+1]-eQCumPop[i])*eQCumWealth[i])+((eQCumPop[i+1]-eQCumPop[i])*eQCumWealth[i+1]))/2 ;
		}
		
		//inequality
		double maxWeight = 0;
		for(int j =0; j<g.length; j++)
		{
			if(g[j][0] > maxWeight)
				maxWeight=g[j][0];
		}
		
		/**
		//casting weights to ints may put all to zero in some data types?
		**/
		
		int[] unequal = new int[(int)maxWeight+1];
		for(int i1=0; i1<g.length-1; i1++)
		{
			int o = (int)g[i1][0];
			unequal[o] ++;
		}
		
		double uETotalWealth = 0;
		
		for(int i1=0; i1<unequal.length; i1++)
		{
			//if(unequal[i1]>0)
				//System.out.println("unequal:" + unequal[i1]);
			uETotalWealth += i1*unequal[i1];
		}
		//System.out.println("use unequal "+ uETotalWealth);
		
		//System.out.println(uETotalWealth);

		double[] uECumPop = new double[unequal.length];
		double[] uECumWealth = new double[unequal.length];
		
		uECumPop[0] = unequal[0]/(double)nodeList.size();
		uECumWealth[0] = 0;
		
		for(int i1 = 1; i1<unequal.length; i1++)
		{
			uECumPop[i1] = unequal[i1]/(double)nodeList.size() + uECumPop[i1-1];
			
			uECumWealth[i1] = (unequal[i1]*i1)/uETotalWealth + uECumWealth[i1-1];
		}
		
		//Lorenz area
		double lArea = 0;
		for(int i1=0; i1<unequal.length-1; i1++)
		{
			lArea += (((uECumPop[i1+1]-uECumPop[i1])*uECumWealth[i1])+((uECumPop[i1+1]-uECumPop[i1])*uECumWealth[i1+1]))/2;
		}	
		//difference in area
		double dArea = eArea - lArea;
		gini = dArea / (lArea + dArea);
		g[nodes+1][0] = gini;
		gini=0;
	}
	public ArrayList<String> inputRead(String file)
	{
		ArrayList<String> StringMat = new ArrayList<String>();
		try 
		{
			@SuppressWarnings("resource")
			Scanner in = new Scanner(new FileReader(System.getProperty("user.dir")+"/"+fold+"/"+file));
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
			@SuppressWarnings("resource")
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
		//System.out.println(dataPoints.size());
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
	public void writeFile(ArrayList<String> writer)
	{
		BufferedWriter b;
		try 
		{
			FileWriter ff;
			if(pics)
				ff = new FileWriter(der+ "/"+"Lorenz_analysis.txt",true);       //File Naming System needed
			else
				ff = new FileWriter("Lorenz_Analysis.txt",true);
			b = new BufferedWriter(ff);
			PrintWriter printer = new PrintWriter(b);
			for(int i = 0; i< writer.size(); i++)
			{
				printer.println(writer.get(i));
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

