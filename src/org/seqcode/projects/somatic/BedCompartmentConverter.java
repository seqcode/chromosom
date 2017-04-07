package org.seqcode.projects.somatic;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;

public class BedCompartmentConverter 
{
	public static void main(String[] args)
	{
		BufferedReader br = null;
		String file = "GSE63525_GM12878_subcompartments.bed";
		String line = "";
		String der = file+"_compartments";
		try 
		{
			br = new BufferedReader(new FileReader(System.getProperty("user.dir")+"/"+file));
			
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		ArrayList<String[]> s = new ArrayList<String[]>();
	    try {
            line = br.readLine();
	        while (line != null) 
	        {
	        	if(line.contains("\t"))
	        	{
		            s.add(line.split("\t"));
	        	}
	        	line = br.readLine();
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
	    ArrayList<String> compartments = new ArrayList<String>();
	    for(int i = 0; i< s.size(); i++)
	    {
	    	if(!compartments.contains(s.get(i)[3]))
	    	{
	    		compartments.add(s.get(i)[3]);
	    		System.out.println(s.get(i)[3]);
	    	}
	    }
	    new File(der).mkdir();
	    for(int i = 0; i < compartments.size(); i++)
	    {
		    BufferedWriter b;
			try 
			{
				String g = der+"/"+file + "_" + compartments.get(i)+".txt";
				FileWriter ff = new FileWriter(g,true);
				b = new BufferedWriter(ff);
				PrintWriter printer = new PrintWriter(b);
				for(int j = 0; j<s.size(); j++)
				{
					printer.println(s.get(j)[0]+":"+s.get(j)[1]+"-"+s.get(j)[2]);
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
}
