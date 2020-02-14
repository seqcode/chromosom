package org.seqcode.projects.chromosom;

public class Region {

	public String chr;
	public int start;
	public int end;
	public char strand='+';
	public double score = 1;
	
	public Region (String c, int s, int e, char str, double sc) {
		chr = c;
		start = s;
		end = e;
		strand = str;
		score = sc;
	}
	public Region (String c, int s, int e, char str) { this(c, s, e, str, 1.0);}
	public Region (String c, int s, int e) { this(c, s, e, '+', 1.0);}
	public Region (BEDLine bed) {
		chr = bed.getChrom();
		start = bed.getChromStart();
		end = bed.getChromEnd();
		strand = bed.getStrand() == null ? '+' : bed.getStrand();
		score = bed.getScore() == null ? 1.0 : bed.getScore();
	}
	
	public int getWidth() {return (end-start+1);}
	
}
