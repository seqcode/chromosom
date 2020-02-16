package org.seqcode.projects.chromosom;
import java.awt.Color;
import java.util.ArrayList;


public class SOMAnalysis
{
	public UseMap u;
	
	public SOMAnalysis()
	{	 

	}
	
	public void useSOM(String mapFileName, String searchDir, String outDir, boolean equalWeights, boolean flipColors){
		System.setProperty("java.awt.headless", "true");
		UseMap u = new UseMap(mapFileName, searchDir, outDir, true, equalWeights);
		if(flipColors) {
			u.setColorA(Color.blue);
			u.setColorB(Color.red);
		}
		u.searchSystem();
	}
	
	public void useSOMpairwise(String mapFileName, String searchDir, String outDir, boolean equalWeights, boolean flipColors){
		System.setProperty("java.awt.headless", "true");
		UseMap u = new UseMap(mapFileName, searchDir, outDir, true, equalWeights);
		if(flipColors) {
			u.setColorA(Color.blue);
			u.setColorB(Color.red);
		}
		u.searchSystemDouble();
	}
	
	public void corSOM(String mapFileName, String searchDir, String outDirectory){
		System.setProperty("java.awt.headless", "true");
		ArrayList<String> folders = new ArrayList<String>();
		folders.add(searchDir);
		Correlations c = new Correlations(mapFileName, folders, outDirectory);
		c.searchSystem();
	}
	
	
}

