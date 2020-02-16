package org.seqcode.projects.chromosom;

import java.awt.BorderLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.util.ArrayList;

import javax.swing.Box;
import javax.swing.JButton;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JMenuBar;
import javax.swing.JTextField;
import javax.swing.filechooser.FileNameExtensionFilter;

public class SOMViewer extends JFrame
{
	int winW, winH;
	public DrawHex d;
	public DrawCompartment dd;
	public UseMap u;
	private static final long serialVersionUID = 1L;
	public ArrayList<String> s;
	public String find;
	
	public SOMViewer(int x, int y, String p)
	{	 
		//title the window
	    super(p);
	    
	    BorderLayout gui = new BorderLayout();
	    setLayout(gui);
	    winW = x; winH = y;
	    
	    s = new ArrayList<String>();
	}
	
	/**
	 * Set up the SOM GUI
	 * @param useWeights
	 */
	public void initViewer(boolean useWeights){
		String trained = "";
		JFileChooser chooser = new JFileChooser();
	    FileNameExtensionFilter filter = new FileNameExtensionFilter(
	        "SOM files", "som");
	    chooser.setFileFilter(filter);
	    chooser.setCurrentDirectory(new File(System.getProperty("user.dir")));
	    chooser.setVisible(true);
	    int returnVal = chooser.showOpenDialog(this);
	    if(returnVal == JFileChooser.APPROVE_OPTION) 
	    	trained = chooser.getSelectedFile().getPath();
	    else if(returnVal != JFileChooser.CANCEL_OPTION){
	    	return;
	    }else{
	    	System.exit(1);
	    }
	    
	    JMenuBar menubar = new JMenuBar();
	    this.setJMenuBar(menubar);
	    int a = 0;
	    this.drawGrid(trained, a, a, useWeights);
	         
	    this.setUpMenu(menubar, this);
	    
	    this.setBounds(0,0,700,700);
	    this.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
	    this.setVisible(true);
	    this.setResizable(true);
	}
	
	/**
	 * Set up the SOM compartment viewer
	 * @param useWeights
	 */
	public void initCompartments(String somFilename, String searchFiles, boolean nodes){
		this.drawComps(somFilename, searchFiles, nodes);
		this.setBounds(0,0,700,700);
		this.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		this.setVisible(true);
		this.setResizable(true);	
	}
	
	//For viewing
	protected void drawGrid(String s)
	{
		d = new DrawHex(s);
		add(d, BorderLayout.CENTER);
		d.nodeBuild(700,700);
    	//d.colors();
    	d.heatMapping();
	}
	protected void drawComps(String s, String ss, boolean nodes)
	{
		 dd = new DrawCompartment(s, ss, nodes);
		 add(dd, BorderLayout.CENTER);
		 dd.nodeBuild(750,750);
		 dd.heatMapping();
		 dd.swap();
		 dd.repaint();
	}
	protected void drawGrid(String s, int i, int j, boolean boo)
	{
		 d = new DrawHex(s,i,j,boo);
		add(d, BorderLayout.CENTER);
		d.nodeBuild(700,700);
    	//d.colors();
    	d.heatMapping();
    	d.countingDPS(null);
	}
	//Finds and shows chromes
	protected void setUpMenu(JMenuBar menubar, SOMViewer self)
	{
		 JButton save = new JButton("Save");
	     menubar.add(save);
		 JButton search = new JButton("Search");
	     menubar.add(search);
	     JButton search2 = new JButton("Search 2");
	     menubar.add(search2);
	     JButton view = new JButton("View Swap");
	     menubar.add(view);
	     JButton button = new JButton("chr");
	     menubar.add(button);
	     
	     
        final JTextField texter = new JTextField("chr...", 20);
	    menubar.add(texter);
	    menubar.add(Box.createHorizontalGlue());
	    
	    button.addActionListener(new ActionListener() {
	    	 
            public void actionPerformed(ActionEvent e)
            {
                //Execute when button is pressed
                find = texter.getText();

                d.countingDPS(find);
            }
        });
	    save.addActionListener(new ActionListener() {
	    	 
            public void actionPerformed(ActionEvent e)
            {
                d.saveImg(texter.getText(),new File(System.getProperty("user.dir")));
            }
        });
	    search.addActionListener(new ActionListener() {
	    	 
            public void actionPerformed(ActionEvent e)
            {
                //Execute when search is pressed
            	File trained=null;
 	    		JFileChooser chooser = new JFileChooser();
 	    	    FileNameExtensionFilter filter = new FileNameExtensionFilter(
 	    	        "TXT files", "txt");
 	    	    chooser.setFileFilter(filter);
 	    	    chooser.setCurrentDirectory(new File(System.getProperty("user.dir")));
 	    	    int returnVal = chooser.showOpenDialog(self);
 	    	    if(returnVal == JFileChooser.APPROVE_OPTION) {
 	    	    	trained = chooser.getSelectedFile();
 	    	    }
 	    	if(d != null && trained !=null)
 	    		d.search(trained);
            }
        });
	    
	    search2.addActionListener(new ActionListener() {
	    	 
            public void actionPerformed(ActionEvent e)
            {
                //Execute when search is pressed
            	File trained=null, trained2=null;
 	    		JFileChooser chooser = new JFileChooser();
 	    	    FileNameExtensionFilter filter = new FileNameExtensionFilter(
 	    	        "TXT files", "txt");
 	    	    chooser.setFileFilter(filter);
 	    	    chooser.setCurrentDirectory(new File(System.getProperty("user.dir")));
 	    	    int returnVal = chooser.showOpenDialog(self);
 	    	    if(returnVal == JFileChooser.APPROVE_OPTION) {
 	    	    	trained = chooser.getSelectedFile();
 	    	    }
 	    	    //Execute when search is pressed
 	    		chooser = new JFileChooser();
 	    	    filter = new FileNameExtensionFilter(
 	    	        "TXT files", "txt");
 	    	    chooser.setFileFilter(filter);
 	    	    chooser.setCurrentDirectory(new File(System.getProperty("user.dir")));
 	    	    returnVal = chooser.showOpenDialog(self);
 	    	    if(returnVal == JFileChooser.APPROVE_OPTION) {
 	    	    	trained2 = chooser.getSelectedFile();
 	    	    }
 	    	if(d != null && trained!=null && trained2!=null)
 	    		d.multiSearch(trained, trained2);
            }
        });
	    
	    view.addActionListener(new ActionListener() {
	    	 
            public void actionPerformed(ActionEvent e)
            {
                //Execute when button is pressed
                d.swap();
            }
        });  
	}

	public static void main(String[] args) {
		ChromoSOM som = new ChromoSOM(args, true);
		som.view();
	}
	
}

