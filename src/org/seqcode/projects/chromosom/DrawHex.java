package org.seqcode.projects.chromosom;
import java.awt.AWTException;
import java.awt.Color;
import java.awt.Font;
import java.awt.GradientPaint;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Rectangle;
import java.awt.Robot;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Scanner;
import javax.imageio.ImageIO;
import javax.swing.JPanel;

public class DrawHex extends JPanel
{
	private static final long serialVersionUID = 1L;
	public double gini;
	public int nodes,xs,ys,xNodes,yNodes, maxDataPoints, minDataPoints, winW, winH, binSize;
	public ArrayList<DataPoint> dataPoints;
	public ArrayList<MiniNode> nodeList;
	public MiniSystem nodeSystem;
	public boolean multi, weighting, norm, equalWeight;
	int swap, col, www, xMin, yMin, xMax, yMax, fontSize=10;
	public double sep, cutoff;
	public Color colorA=Color.red, colorB=Color.blue; //Colors for dataset A and dataset B (if plotting two colors) 
	
	public DrawHex(String s)
	{
		cutoff = 0.01;
		multi = false;
		norm = false;//true;
		swap=0;
		gini=0;
		sep = 0;
		weighting = false;
		String ffs = s;
		yNodes=0;
		xNodes=0;
		reader(ffs);
		nodes = yNodes*xNodes;
    	binSize = dataPoints.get(0).maxLocus - dataPoints.get(0).minLocus;
    	col = 0;
    	www = 1;
    	equalWeight = false;
	}
	public DrawHex(String s, int colLoc, int weightLoc, boolean equalW)
	{
		cutoff = 0.01;
		multi = false;
		norm = false;//true;
		swap=0;
		gini=0;
		sep = 0;
		weighting = false;
		String ffs = s;
		yNodes=0;
		xNodes=0;
		reader(ffs);
		nodes = yNodes*xNodes;
    	binSize = dataPoints.get(0).maxLocus - dataPoints.get(0).minLocus;
    	col = colLoc;
    	www = weightLoc;
    	equalWeight = equalW;
	}
	public void setEqualWeight(boolean ew){equalWeight=ew;}
	public void setColorA(Color c) {colorA=c;}
	public void setColorB(Color c) {colorB=c;}
	
	public void saveImg(String name, File outDir)
	{
		BufferedImage image = null;
		try {
			image = new Robot().createScreenCapture(new Rectangle(this.getLocationOnScreen().x+62, this.getLocationOnScreen().y+30, 511, 457));
		} catch (AWTException e1) {
			e1.printStackTrace();
		}
	    try 
	    {
	        File outputfile = new File(outDir.getAbsolutePath()+"/"+name+".png");
	        ImageIO.write(image, "png", outputfile);
	    } catch (IOException e) {}
	}
	public void saveImg(String name, File outDir, int w, int h)
	{
		this.setVisible(true);
		this.setSize(w, h);
	    BufferedImage bi = new BufferedImage(w, h, BufferedImage.TYPE_INT_RGB);
	    Graphics2D g = bi.createGraphics();
	    print(g);
	    try 
	    {
	        File outputfile = new File(outDir.getAbsolutePath()+"/"+name+".png");
	        ImageIO.write(bi, "png", outputfile);
	    } catch (IOException e) {}
	}
	public void swap()
	{
		if(swap >= 3)
		{
			swap = 0;
			if(multi)
				multiHeatMapping();
			else
				heatMapping();
		}
		else
		{
			swap ++;
			if(multi)
				multiHeatMapping();
			else
				heatMapping();
		}
		repaint();
	}
	public void multiSearch(File file1, File file2)
	{
		multi = true;
		if(file1.equals(file2))
		{
			search(file1);
			return;
		}
		for(int i = 0; i < nodeList.size(); i++)
		{
			MiniNode mini = nodeList.get(i);
			mini.counting.clear();
			mini.blueCounting.clear();
			mini.blueWeight=0;
			mini.weight = 0;
		}
		BEDReader reader;
		try {
			
			//File 1
			reader = new BEDReader(file1);
		
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
				for(int i = 0; i< dataPoints.size(); i++){
					//Overlap assessment 
					if(dataPoints.get(i).chrome.equals(reg.chr) && 
							((dataPoints.get(i).minLocus <= reg.start && dataPoints.get(i).maxLocus >= reg.start) ||
							(reg.start <= dataPoints.get(i).minLocus && reg.end >= dataPoints.get(i).minLocus )) ){
								count++;
								inds.add(i);
								
								//Check for contained relationship to shortcut
								if(dataPoints.get(i).chrome.equals(reg.chr) && 
										(dataPoints.get(i).minLocus <= reg.start && dataPoints.get(i).maxLocus >= reg.end) ) {
									break;
								}
							}
				}

				if(!equalWeight)
					weight/=count;
				for(Integer i: inds)
				{
					dataPoints.get(i).myMini.counting.add(dataPoints.get(i));
					dataPoints.get(i).myMini.weight += weight;
				}
			}
			
			//File 2
			reader = new BEDReader(file2);
			
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
				for(int i = 0; i< dataPoints.size(); i++){
					//Overlap assessment 
					if(dataPoints.get(i).chrome.equals(reg.chr) && 
							((dataPoints.get(i).minLocus <= reg.start && dataPoints.get(i).maxLocus >= reg.start) ||
							(reg.start <= dataPoints.get(i).minLocus && reg.end >= dataPoints.get(i).minLocus )) ){
								count++;
								inds.add(i);
								
								//Check for contained relationship to shortcut
								if(dataPoints.get(i).chrome.equals(reg.chr) && 
										(dataPoints.get(i).minLocus <= reg.start && dataPoints.get(i).maxLocus >= reg.end) ) {
									break;
								}
							}
				}

				if(!equalWeight)
					weight/=count;
				for(Integer i: inds)
				{
					dataPoints.get(i).myMini.blueCounting.add(dataPoints.get(i));
					dataPoints.get(i).myMini.blueWeight += weight;
				}
			}
			
			multiHeatMapping();
			gini();
			//degreesOfSep();
			repaint();

		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	public void search(File file)
	{
		multi = false;
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
				for(int i = 0; i< dataPoints.size(); i++){
					//Overlap assessment 
					if(dataPoints.get(i).chrome.equals(reg.chr) && 
							((dataPoints.get(i).minLocus <= reg.start && dataPoints.get(i).maxLocus >= reg.start) ||
							(reg.start <= dataPoints.get(i).minLocus && reg.end >= dataPoints.get(i).minLocus )) ){
								count++;
								inds.add(i);
								
								//Check for contained relationship to shortcut
								if(dataPoints.get(i).chrome.equals(reg.chr) && 
										(dataPoints.get(i).minLocus <= reg.start && dataPoints.get(i).maxLocus >= reg.end) ) {
									break;
								}
							}
				}

				if(!equalWeight)
					weight/=count;
				for(Integer i: inds)
				{
					dataPoints.get(i).myMini.counting.add(dataPoints.get(i));
					dataPoints.get(i).myMini.weight += weight;
				}
			}
			
			heatMapping();
			gini();
			degreesOfSep();
			repaint();

		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	public void search(File file, int locCol, int weightCol)
	{
		multi = false;
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
				for(int i = 0; i< dataPoints.size(); i++){
					//Overlap assessment 
					if(dataPoints.get(i).chrome.equals(reg.chr) && 
							((dataPoints.get(i).minLocus <= reg.start && dataPoints.get(i).maxLocus >= reg.start) ||
							(reg.start <= dataPoints.get(i).minLocus && reg.end >= dataPoints.get(i).minLocus )) ){
								count++;
								inds.add(i);
								
								//Check for contained relationship to shortcut
								if(dataPoints.get(i).chrome.equals(reg.chr) && 
										(dataPoints.get(i).minLocus <= reg.start && dataPoints.get(i).maxLocus >= reg.end) ) {
									break;
								}
							}
				}

				if(!equalWeight)
					weight/=count;
				
				for(Integer i: inds){
					dataPoints.get(i).myMini.counting.add(dataPoints.get(i));
					dataPoints.get(i).myMini.weight += weight;
				}
			}
		
			heatMapping();
			gini();
			degreesOfSep();
			repaint();

		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	public void gini()
	{

		double[] unequal = new double[maxDataPoints+1];;
		if(weighting)
		{
			for(int i=0; i<nodeList.size(); i++)
			{
				int o = (int) nodeList.get(i).weight;
				unequal[o]++;
			}
		}
		else
		{
			for(int i=0; i<nodeList.size(); i++)
			{
				int o = (int)(nodeList.get(i).counting.size());
				unequal[o]++;
			}
		}
		int uETotalWealth = 0;
		for(int i=0; i<unequal.length; i++)
		{
			uETotalWealth += i*unequal[i];
		}
		double[] uECumPop = new double[maxDataPoints+1];
		double[] uECumWealth = new double[maxDataPoints+1];
		uECumPop[0] = unequal[0]/nodeList.size();
		uECumWealth[0] = 0;
		for(int i = 1; i<unequal.length; i++)
		{
			uECumPop[i] = unequal[i]/nodeList.size() + uECumPop[i-1];
			uECumWealth[i] = (unequal[i]*i)/uETotalWealth + uECumWealth[i-1];
		}
		
		//equality baseline
		int mdp = (int) (nodeList.get(0).bins.size());
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
		double[] eQCumPop = new double[mdp+1];
		double[] eQCumWealth = new double[mdp+1];
		eQCumPop[0] = equal[0]/nodeList.size();
		eQCumWealth[0] = 0;
		for(int i = 1; i<equal.length; i++)
		{
			eQCumPop[i] = equal[i]/nodeList.size() + eQCumPop[i-1];
			eQCumWealth[i] = (equal[i]*i)/eQTotalWealth + eQCumWealth[i-1];
		}
		
		//Lorenz area
		double lArea = 0;
		for(int i=0; i<unequal.length-1; i++)
		{
			lArea += (((uECumPop[i+1]-uECumPop[i])*uECumWealth[i])+((uECumPop[i+1]-uECumPop[i])*uECumWealth[i+1]))/2 ;
		}
		
		//Equaliy area
		double eArea = 0;
		for(int i=0; i<equal.length-1; i++)
		{
			eArea += (((eQCumPop[i+1]-eQCumPop[i])*eQCumWealth[i])+((eQCumPop[i+1]-eQCumPop[i])*eQCumWealth[i+1]))/2 ;
		}
				
		//difference in area
		double dArea = eArea - lArea;
		gini = dArea / (lArea + dArea);
	}
	
	public ArrayList<String> inputRead(File file)
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
	public void countingDPS(String chr)
	{
		multi = false;
		weighting = false;
		for(int i =0; i<nodeSystem.size();i++)
		{
			MiniNode mini = nodeSystem.get(i);
			mini.counting.clear();
			for(int j =0; j<mini.bins.size(); j++)
			{
				//if(chr <= -1)
				//{
				//	if(Math.random()>.9)
				//		mini.counting.add(mini.bins.get(j));
				//}
				//else 
				if (chr == null)
					mini.counting.add(mini.bins.get(j));
				else if (mini.bins.get(j).chrome.contentEquals(chr))
					mini.counting.add(mini.bins.get(j));
			}
		}
	    maxDataPoints = 0;
	    minDataPoints = 0;
		heatMapping();
		gini();
		degreesOfSep();
		repaint();
	}
	public void reader(String f)
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
		dataPoints = new ArrayList<DataPoint>();
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
				dataPoints.add(nodeList.get(i).bins.get(j));
			}
		}
		nodeSystem = new MiniSystem(nodeList,xNodes,yNodes);
	}
	public void paintComponent(Graphics g)
	  {
	    super.paintComponent(g);
	    colorBar(g);
	    //setBackground(Color.BLUE);
	    
//	    int centerNode = 2450;
	    for(int i = 0; i<nodeList.size(); i++)
	    {	
//		 	//For determining if neighbor Distance is working	    	
//	    	if(i == centerNode)
//	    	{
//		    	int xCoordA = i%xNodes;
//				int yCoordA = i/yNodes;
//				
//				for (int j = 0; j<nodeList.size(); j++)
//				{
//					int xCoordB = j%xNodes;
//					int yCoordB = j/yNodes;
//					
//					double ddd = hexDist(xCoordA, yCoordA, xCoordB, yCoordB);
//
//					if(ddd == 0)
//						nodeList.get(j).color = Color.BLACK;
//					if(ddd == 1)
//						nodeList.get(j).color = Color.BLUE;
//					if(ddd == 2)
//						nodeList.get(j).color = Color.GREEN;
//					if(ddd == 3)
//						nodeList.get(j).color = Color.RED;
//					if(ddd == 4)
//						nodeList.get(j).color = Color.BLUE;
//					if(ddd == 5)
//						nodeList.get(j).color = Color.GREEN;
//					if(ddd == 6)
//						nodeList.get(j).color = Color.RED;
//					if(ddd == 7)
//						nodeList.get(j).color = Color.BLUE;
//					if(ddd == 8)
//						nodeList.get(j).color = Color.GREEN;
//					if(ddd == 9)
//						nodeList.get(j).color = Color.RED;
//					if(ddd == 10)
//						nodeList.get(j).color = Color.BLUE;
//					if(ddd == 11)
//						nodeList.`get(j).color = Color.GREEN;
//					if(ddd == 12)
//						nodeList.get(j).color = Color.RED;
//				}
//	    	}

	    	g.setColor(Color.WHITE);
	    	g.fillPolygon(nodeList.get(i).p);
	    	if(multi)
	    	{
	    		g.setColor(nodeList.get(i).color2);
	    		g.fillPolygon(nodeList.get(i).p);	
	    		//nodeList.get(i).color2 = null;
	    	}
    		g.setColor(nodeList.get(i).color);
	    	g.fillPolygon(nodeList.get(i).p);
	    	
	    	if(swap==2||swap ==3)
	    	{
	    		g.setColor(Color.BLACK);
		    	g.setFont(new Font("Arial", Font.PLAIN, fontSize));
		    	g.drawString(""+nodeList.get(i).counting.size(),nodeList.get(i).xLoc-4, nodeList.get(i).yLoc+5);
		    	g.drawPolygon(nodeList.get(i).p);
	    	}
	    }
	    g.setColor(Color.BLACK);
    	g.setFont(new Font("Arial", Font.PLAIN, fontSize));
    	String giniStr = String.format("Gini = %.4f", gini);
    	g.drawString(giniStr, (getWidth()/2)-(g.getFontMetrics().stringWidth(giniStr)/2), (yMin/2)-(g.getFontMetrics().getHeight()/2));
    	//g.drawString("Sep = "+ sep, getWidth()-400, 15);
	}
	public void nodeBuild(int winW, int winH)
	{  
		xMin = (int)(.1 * winW);
	    yMin = (int)(.1 * winH);
	    xMax = (int)(.9 * winW);
	    yMax = (int)(.9 * winH);
	    fontSize = winH/85;
	    
		
	    double winWidth = xMax-xMin;
	    double winHeight = yMax-yMin;
	    
	    double height = winHeight/yNodes;
	    double width = winWidth/xNodes;
	    height = height/3;
	    width = width/2;
	    int x = (int) (width);
	    int y = (int) (height);
	    int w = (int) (winWidth/width);
	    int h = (int)(winHeight/height);
	    for(int i = 0; i<h; i++)
		{
	    	for(int m = 0; m<w; m++)
			{
				if((i%6==1 && m%2==0)||(i%6==4 && m%2 == 1))
				{
					MiniNode p = nodeSystem.get(m/2,i/3);
					p.xLoc = m*x+xMin;
					p.yLoc = i*y+yMin;
				}
			}
		}
	   for(int c = 0; c<nodeList.size(); c++)
	   {
		  MiniNode m = nodeList.get(c);
		  
		  int[] xPoints = new int[6];
		  int[] yPoints = new int[6];
		  
		  xPoints[0] = m.xLoc;
		  yPoints[0] = m.yLoc-(2*y);
		  
		  xPoints[1] = m.xLoc+x;
		  yPoints[1] = m.yLoc-y;
		  
		  xPoints[2] = m.xLoc+x;
		  yPoints[2] = m.yLoc+y;
			
		  xPoints[3] = m.xLoc;
		  yPoints[3] = m.yLoc+(2*y);

		  xPoints[4] = m.xLoc-x;
		  yPoints[4] = m.yLoc+y;
			
		  xPoints[5] = m.xLoc-x;
		  yPoints[5] = m.yLoc-y;
		  
		  m.polygonMaker(xPoints, yPoints, 6);
	   }
	    
	    
//	    for(int c = 0; c<nodeList.size(); c++)
//		   {
//			  MiniNode m = nodeList.get(c);
//			  
//			  int[] xPoints = new int[4];
//			  int[] yPoints = new int[4];
//			  
//			  xPoints[0] = m.xLoc+x;
//			  yPoints[0] = m.yLoc-y;
//			  
//			  xPoints[1] = m.xLoc+x;
//			  yPoints[1] = m.yLoc+y;
//
//			  xPoints[2] = m.xLoc-x;
//			  yPoints[2] = m.yLoc+y;
//				
//			  xPoints[3] = m.xLoc-x;
//			  yPoints[3] = m.yLoc-y;
//			  
//			  m.polygonMaker(xPoints, yPoints, 4);
//		   }
	}
	public void degreesOfSep()
	{
		double sepp = 0;
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
			double count = 0;
			for(int j = 0; j<nodes; j++)
			{
				for(int k = 0; k<nodes; k++)
				{
					sepp += dist[j][k]*(nodeList.get(k).counting.size())*(nodeList.get(j).counting.size());
					count+= nodeList.get(k).counting.size()*(nodeList.get(j).counting.size());
				}
			}
			sepp /= count;
			sep = sepp;
			
	}
	public int hexDist(int x1, int y1, int x2, int y2)
	{
		if(x1 == x2 && y1 == y2)
			return 0;
		else
		{
			int hm = 0;
			int vm = 0;
			int cm = 0;
			boolean right = x1-x2<0;
			
			hm = Math.abs(x1-x2);
			vm = Math.abs(y1-y2);
			if(xNodes-x1<xNodes-x2 && xNodes-x1 + x2 < hm)
			{
				right = true;
				hm =  xNodes-x1 + x2;	
			}
			else if (xNodes-x1>xNodes-x2 && xNodes-x2 +x1 < hm)
			{
				right = false;
				hm = xNodes-x2 +x1;
			}
			
			if(yNodes-y1<yNodes-y2 && yNodes-y1 + y2 < vm)
			{
				vm = yNodes-y1 + y2;
			}
			else if(yNodes-y1>yNodes-y2 && yNodes-y2 + y1 < vm)
			{
				vm = yNodes-y2 + y1;
			}
			
			if(vm %2 !=0)
			{
				if(!right && y1%2 == 0)
				{
					cm = (int) Math.min(hm, (((double)vm)/2 +.5));
				}
				else if(!right && y1%2 == 1)
				{
					cm = (int) Math.min(hm, (((double)vm)/2));
				}
				else if(right && y1%2 == 0)
				{
					cm = (int) Math.min(hm, (((double)vm)/2));
				}
				else //if(right && y1%2 == 1)
				{
					cm = (int) Math.min(hm, (((double)vm)/2 +.5));
				}
			}
			else
				cm = (int) Math.min(hm, ((((double)vm)/2)));
			return hm+vm-cm; 
		}
	}
//	public void colors()
//	{
//		for(int k = 0; k<colorNum;k++)
//		{
//			int red = 255;//(int)(k*255/colorNum);
//			int blue = (int)(255-(k*255/colorNum));
//			int green = (int)(255-(k*255/colorNum));
//			Color c = new Color(red, green, blue);
//			colors.add(c);
//		}
//	}
	public void colorBar(Graphics g)
	{
		Graphics2D g2d = (Graphics2D)g;
		int x = (int)(getWidth()*.2);
		int y = (int)(getHeight()*.96);
		int width = (int)(getWidth()*.55);
		int height = (int)(getHeight()*.02);
		
		GradientPaint colorbar = new GradientPaint(x, y, Color.white, x+width, y, colorA, false);
		g2d.setPaint(colorbar);
		g2d.fillRect(x, y, width, height);
		
		g2d.setPaint(Color.black);
		g2d.setColor(Color.BLACK);
		g2d.drawRect(x, y, width, height);
		g2d.drawString(""+minDataPoints, (int)(getWidth()*.17),(int)(getHeight()*.95));
		g2d.drawString(""+maxDataPoints, (int)(getWidth()*.2+(int)(getWidth()*.55)),(int)(getHeight()*.95));
		
	}
	public void heatMapping()
	{
		minDataPoints = 0;
		maxDataPoints = 0;
//		if(weighting == true && norm)
//		{
//			minDataPoints = 0;//(int) (nodeList.get(0).counting.size() * nodeList.get(0).weight);
//			maxDataPoints = 0;//(int) (nodeList.get(0).counting.size() * nodeList.get(0).weight);
//			
//			for(int i = 0; i<nodeList.size(); i++)
//			{
//				double b = 1;
//				if(nodeList.get(i).bins.size()>0)
//					b = nodeList.get(i).bins.size();
//				if(nodeList.get(i).counting.size()* nodeList.get(i).weight/b > maxDataPoints)
//					maxDataPoints = (int) (nodeList.get(i).counting.size() * nodeList.get(i).weight/b);
//				if(nodeList.get(i).counting.size()* nodeList.get(i).weight/b <minDataPoints)
//					minDataPoints = (int) (nodeList.get(i).counting.size()* nodeList.get(i).weight/b);
//			}
//			if(maxDataPoints == minDataPoints) maxDataPoints++;
//			
//			for(int i = 0; i<nodeList.size(); i++)
//			{	   
//				double b = 1;
//				if(nodeList.get(i).bins.size()>0)
//					b = nodeList.get(i).bins.size();
//				MiniNode p = nodeList.get(i);
//				double d = p.counting.size()*p.weight/b;
//				double doop = maxDataPoints;
//				p.color = new Color(255, (int)(255-(d*255/doop)),(int)(255-d*255/doop));
//			}
//		}
//		else 
		if(weighting == true && equalWeight==false)
		{
			minDataPoints = (int) (nodeList.get(0).counting.size() * nodeList.get(0).weight);
			maxDataPoints = (int) (nodeList.get(0).counting.size() * nodeList.get(0).weight);
			for(int i = 0; i<nodeList.size(); i++)
			{
				if(nodeList.get(i).counting.size()* nodeList.get(i).weight>maxDataPoints)
					maxDataPoints = (int) (nodeList.get(i).counting.size() * nodeList.get(i).weight);
				if(nodeList.get(i).counting.size()* nodeList.get(i).weight<minDataPoints)
					minDataPoints = (int) (nodeList.get(i).counting.size()* nodeList.get(i).weight);
			}
			if(maxDataPoints == minDataPoints) maxDataPoints++;
			for(int i = 0; i<nodeSystem.size(); i++)
			{	   
				MiniNode p = nodeSystem.get(i);
				double d = p.counting.size()*p.weight;
				double doop = maxDataPoints;
				p.color = new Color(colorA.getRed(), colorA.getGreen(),colorA.getBlue(),(int)(255*d/doop));
				if((swap == 1||swap == 3)&& p.counting.size()*p.weight<=cutoff*maxDataPoints)
					p.color = Color.GRAY;
			}
		}
		else
		{
			minDataPoints = (int) (nodeList.get(0).counting.size());
			maxDataPoints = (int) (nodeList.get(0).counting.size());
			for(int i = 0; i<nodeList.size(); i++)
			{
				if(nodeList.get(i).counting.size()>maxDataPoints)
					maxDataPoints = (int) (nodeList.get(i).counting.size());
				if(nodeList.get(i).counting.size()<minDataPoints)
					minDataPoints = (int) (nodeList.get(i).counting.size());
			}
			if(maxDataPoints == minDataPoints) maxDataPoints++;
			for(int i = 0; i<nodeSystem.size(); i++)
			{	   
				MiniNode p = nodeSystem.get(i);
				double d = p.counting.size();
				double doop = maxDataPoints;
				p.color = new Color(colorA.getRed(), colorA.getGreen(),colorA.getBlue(),(int)(255*d/doop));
				if((swap == 1||swap == 3) && p.counting.size()==0)
					p.color = Color.GRAY;
			}
		}
	}
	//both files must be weighted
	public void multiHeatMapping()
	{
		if(weighting == true && equalWeight==false)
		{

			maxDataPoints = (int) (nodeList.get(0).counting.size() * nodeList.get(0).weight);
			int maxBlueDataPoints = (int) (nodeList.get(0).blueCounting.size() * nodeList.get(0).blueWeight);

			for(int i = 0; i<nodeList.size(); i++)
			{
				if(nodeList.get(i).counting.size()* nodeList.get(i).weight>maxDataPoints)
					maxDataPoints = (int) (nodeList.get(i).counting.size() * nodeList.get(i).weight);
				if(nodeList.get(i).blueCounting.size()* nodeList.get(i).blueWeight>maxBlueDataPoints)
					maxBlueDataPoints = (int) (nodeList.get(i).blueCounting.size() * nodeList.get(i).blueWeight);
			}
			for(int i = 0; i<nodeList.size(); i++)
			{	   
				MiniNode p = nodeList.get(i);
				double d = p.counting.size()*p.weight;
				double bl = p.blueCounting.size()*p.blueWeight;
				double doop = maxDataPoints;
				double bloop = maxBlueDataPoints;
				double ratTot = d + bl;
				double ratRat = d/ratTot;
				d/=doop;
				bl/=bloop;
				double redRat = d*ratRat;//(d+(.7*bl))/1.7;
				double blueRat = bl*(1-ratRat);//(bl+(.7*d))/1.7;
				p.color = new Color(colorA.getRed(), colorA.getGreen(),colorA.getBlue(), (int)(255*redRat));
				p.color2 = new Color(colorB.getRed(), colorB.getGreen(),colorB.getBlue(),(int)(255*blueRat));
			
				if((swap == 1||swap == 3) && p.counting.size()*p.weight<=cutoff*maxDataPoints && p.blueCounting.size()*p.blueWeight<=cutoff*maxBlueDataPoints)
					p.color = Color.GRAY;
			}
		}
		else
		{
			maxDataPoints = (int) (nodeList.get(0).counting.size());
			int maxBlueDataPoints = (int) (nodeList.get(0).blueCounting.size());
			for(int i = 0; i<nodeList.size(); i++)
			{
				if(nodeList.get(i).counting.size()>maxDataPoints)
					maxDataPoints = (int) (nodeList.get(i).counting.size());
				if(nodeList.get(i).blueCounting.size()>maxBlueDataPoints)
					maxBlueDataPoints = (int) (nodeList.get(i).blueCounting.size());
			}
			for(int i = 0; i<nodeList.size(); i++)
			{	   
				MiniNode p = nodeList.get(i);
				double d = p.counting.size();
				double bl = p.blueCounting.size();
				double doop = maxDataPoints;
				double bloop = maxBlueDataPoints;
				double ratTot = d + bl;
				double ratRat = d/ratTot;
				d/=doop;
				bl/=bloop;
				double redRat = d*ratRat;//(d+(.7*bl))/1.7;
				double blueRat = bl*(1-ratRat);//(bl+(.7*d))/1.7;
				p.color = new Color(colorA.getRed(), colorA.getGreen(),colorA.getBlue(), (int)(255*redRat));
				p.color2 = new Color(colorB.getRed(), colorB.getGreen(),colorB.getBlue(),(int)(255*blueRat));
			
				if((swap == 1||swap == 3) && p.counting.size()*p.weight<=cutoff*maxDataPoints && p.blueCounting.size()*p.blueWeight<=cutoff*maxBlueDataPoints)
					p.color = Color.GRAY;
			}
		}
	}
}
