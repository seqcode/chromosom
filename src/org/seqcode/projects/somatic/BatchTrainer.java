package org.seqcode.projects.somatic;
import java.awt.BorderLayout;
import java.awt.event.*;
import java.io.File;
import java.util.ArrayList;

import javax.swing.*;
import javax.swing.filechooser.FileNameExtensionFilter;


public class BatchTrainer extends JFrame
{
    static BatchTrainer window2;
    int winW, winH;
	public DrawHex d;
	public DrawCompartment dd;
	public UseMap u;
	private static final long serialVersionUID = 1L;
	public ArrayList<String> s;
	public String find;
	public BatchTrainer(int x, int y, String p)
	  {
		 
		//title the window
	    super(p);
	    
	    BorderLayout gui = new BorderLayout();
	    setLayout(gui);
	    winW = x; winH = y;
	    
	    s = new ArrayList<String>();
	  }
	
	
	
	 //For viewing
	 public void drawGrid(String s)
	 {
		d = new DrawHex(s);
		add(d, BorderLayout.CENTER);
		d.nodeBuild(700,700);
    	//d.colors();
    	d.heatMapping();
	 }
	 public void drawComps(String s, String ss, boolean nodes)
	 {
		 dd = new DrawCompartment(s, ss, nodes);
		 add(dd, BorderLayout.CENTER);
		 dd.nodeBuild(750,750);
		 dd.heatMapping();
		 dd.swap();
		 dd.repaint();
	 }
	 public void drawGrid(String s, int i, int j, boolean boo)
	 {
		 d = new DrawHex(s,i,j,boo);
		add(d, BorderLayout.CENTER);
		d.nodeBuild(700,700);
    	//d.colors();
    	d.heatMapping();
    	d.countingDPS(null);
	 }
	 //Finds and shows chromes
	 public void setUpMenu(JMenuBar menubar)
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
 	    	    int returnVal = chooser.showOpenDialog(window2);
 	    	    if(returnVal == JFileChooser.APPROVE_OPTION) {
 	    	    	trained = chooser.getSelectedFile();
 	    	    	//System.out.println("You chose to open this file: " +chooser.getSelectedFile().getName());
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
 	    	    int returnVal = chooser.showOpenDialog(window2);
 	    	    if(returnVal == JFileChooser.APPROVE_OPTION) {
 	    	    	trained = chooser.getSelectedFile();
 	    	    	//System.out.println("You chose to open this file: " +chooser.getSelectedFile().getName());
 	    	    }
 	    	    //Execute when search is pressed
 	    		chooser = new JFileChooser();
 	    	    filter = new FileNameExtensionFilter(
 	    	        "TXT files", "txt");
 	    	    chooser.setFileFilter(filter);
 	    	    chooser.setCurrentDirectory(new File(System.getProperty("user.dir")));
 	    	    returnVal = chooser.showOpenDialog(window2);
 	    	    if(returnVal == JFileChooser.APPROVE_OPTION) {
 	    	    	trained2 = chooser.getSelectedFile();
 	    	    	//System.out.println("You chose to open this file: " +chooser.getSelectedFile().getName());
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
	 public static void main(String[] args)
	  {
		 //to view: args = "view"
	    	
	    	
	    	if(args[0].equalsIgnoreCase("view"))
	    	{
	    		window2 = new BatchTrainer(700,700, "view");

	    	    String trained = "";
	    		JFileChooser chooser = new JFileChooser();
	    		System.out.println("1");
	    	    FileNameExtensionFilter filter = new FileNameExtensionFilter(
	    	        "TXT files", "txt");
	    	    System.out.println("2");
	    	    chooser.setFileFilter(filter);
	    	    System.out.println("3");
	    	    chooser.setCurrentDirectory(new File(System.getProperty("user.dir")));
	    	    System.out.println("4");
	    	    chooser.setVisible(true);
	    	    int returnVal = chooser.showOpenDialog(window2);
	    	    if(returnVal == JFileChooser.APPROVE_OPTION) 
	    	    {
	    	    	trained = chooser.getSelectedFile().getPath();
	    	    	//System.out.println("You chose to open this file: " +chooser.getSelectedFile().getName());
	    	    }
	    	    else if(returnVal != JFileChooser.CANCEL_OPTION)
	    	    {
	    	    	System.out.println("okay, nevermind");
	    	    	return;
	    	    }
	    	    System.out.println("6");
	    	    JMenuBar menubar = new JMenuBar();
	    	    window2.setJMenuBar(menubar);
	    	    try
	    	    {
	    	        int a = Integer.parseInt(args[1]);
	    	        int b = Integer.parseInt(args[2]);
	    	        if(b==-1)
	    	        	window2.drawGrid(trained, a, a, b == -1);
	    	        else
	    	        	window2.drawGrid(trained, a, b, b == -1);
	    	    }catch(NumberFormatException e)
	    	    {
	    	    	window2.drawGrid(trained);   
	    	    }
	    	         
	    	    window2.setUpMenu(menubar);
	    	    
	    	    window2.setBounds(0,0,700,700);
	    	    window2.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
	    	    window2.setVisible(true);
	    	    window2.setResizable(true);	
	    	}
	    	if(args[0].equalsIgnoreCase("use"))
	    	{
	    		String mapFileName = args[1];
	    		String searchFiles = args[2];
	    		UseMap u = new UseMap(mapFileName, searchFiles, Integer.parseInt(args[3]), Integer.parseInt(args[4]),Integer.parseInt(args[5])==1,Integer.parseInt(args[4])==-1 );
	    		u.searchSystem();
	    	}
	    	if(args[0].equalsIgnoreCase("compartments"))
	    	{
	    		window2 = new BatchTrainer(750,750, "view");
	    		String mapFileName = args[1];
	    		String searchFiles = args[2];
	    		window2.drawComps(mapFileName, searchFiles, args[3]!=null &&args[3].contains("nodes"));
	    		window2.setBounds(0,0,700,700);
	    	    window2.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
	    	    window2.setVisible(true);
	    	    window2.setResizable(true);	
	    	}
	    	if(args[0].equalsIgnoreCase("cor"))
	    	{
	    		String map = args[1];
	    		ArrayList<String> folders = new ArrayList<String>();
	    		for(int i = 2; i<args.length; i++)
	    		{
	    			folders.add(args[i]);
	    		}
	    		Correlations c = new Correlations(map, folders);
	    		c.searchSystem();
	    	}
	    	if(args[0].equalsIgnoreCase("pairwise"))
	    	{
	    		String mapFileName = args[1];
	    		String searchFiles = args[2];
	    		UseMap u = new UseMap(mapFileName, searchFiles, Integer.parseInt(args[3]), Integer.parseInt(args[4]),Integer.parseInt(args[5])==1, Integer.parseInt(args[4])==-1);
	    		u.searchSystemDouble();
	    	}
	}
}

