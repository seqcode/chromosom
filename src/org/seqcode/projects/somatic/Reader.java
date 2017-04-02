package org.seqcode.projects.somatic;

import java.io.*;
import java.util.ArrayList;
import java.util.Scanner;

public class Reader 
{
	public int count;
	public ArrayList<DataPoint> points;
	public String file;
	public String lookup;
	
	public Reader(String fileName)//, String look)
	{
		file = fileName;
		//lookup = look;
		points = new ArrayList<DataPoint>();
		count = 0;
	}
	public void matrixReader()
	{
		BufferedReader br = null;
		
		String line = "";
		try 
		{
			br = new BufferedReader(new FileReader(file));
			
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
	    try {
           // line = br.readLine();
	        while (line != null) 
	        {
	        	line = br.readLine();
	            makePoint(line);
	        }
	    } catch (IOException e)
	    {
			e.printStackTrace();
		} finally 
		{if(br!=null)
			{
	        	try {
	        		br.close();
	        	} catch (IOException e) 
				{
				e.printStackTrace();
				}
			}
		}
	}
	public void makePoint(String s)
	{
		if(s !=null &&s.indexOf("\t")!=-1 &&!s.substring(0, s.indexOf("\t")).isEmpty())
		{
			count++;
			//System.out.println(s);
			DataPoint p = new DataPoint(); p.setName(s.substring(0,s.indexOf("\t")));
			s = s.substring(s.indexOf("\t")+1,s.length());
			//if(s.contains("1")||s.contains("2")||s.contains("3")||s.contains("4")||s.contains("5")||s.contains("6")||s.contains("7")||s.contains("8")||s.contains("9"))
			{
			//System.out.println(s.length());
			String[] temp = s.split("\\t");
			p.gInit(temp.length);
			try {
				for(int t = 0; t< temp.length; t++)
				{
					p.g[t] = ((double)Double.parseDouble(temp[t]));
				}
				points.add(p);
				
			} catch (NumberFormatException e) {
			}
		}
		}
	}
}
