package org.seqcode.projects.chromosom;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Iterator;


public class BEDReader implements Iterator<Region>{
	
	
	private BEDLine nextLine;
	private BufferedReader br=null;
	
	public BEDReader(File f) throws IOException { 
		br = new BufferedReader(new FileReader(f));
		findNextLine();
	}
	
	private void findNextLine() throws IOException { 
		String line = br.readLine();
		if(line != null) { 
			nextLine = new BEDLine(line);
		} else { 
			nextLine = null;
			br.close();
		}
	}
	
	public boolean hasNext() { 
		return nextLine != null;
	}
	
	public Region next() { 
		Region ret = new Region(nextLine);
		try {
			findNextLine();
		} catch (IOException e) {
			e.printStackTrace();
			nextLine = null;
			try {
				br.close();
			} catch (IOException e1) {
				e1.printStackTrace();
			}
		}
		return ret;
	}
	public void close(){
		if(br !=null)
			try {
				br.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
	}
	public void remove() { 
		throw new UnsupportedOperationException();
	}
		
}
