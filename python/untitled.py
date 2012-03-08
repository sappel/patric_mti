/*
 * $Id: CustomContourDemo.java,v 1.1 2008-10-07 19:38:38 rstein Exp $
 * 
 * $Date: 2008-10-07 19:38:38 $ $Revision: 1.1 $ $Author: rstein $
 * 
 * Copyright CERN, All Rights Reserved.
 */
package cern.jdve.demo;

import java.awt.*;
import java.io.BufferedReader;
import java.io.InputStreamReader;
import javax.swing.*;

import cern.jdve.Chart;
import cern.jdve.ChartInteractor;
import cern.jdve.Style;
import cern.jdve.data.DataSet;
import cern.jdve.data.DefaultDataSet;
import cern.jdve.data.DefaultDataSet3D;
import cern.jdve.data.DefaultDataSource;
import cern.jdve.interactor.DataPickerInteractor;
import cern.jdve.renderer.ContourChartRenderer;
import cern.jdve.renderer.PolylineChartRenderer;
import cern.jdve.scale.TimeStepsDefinition;

/**
 * The class is an extension of the simple ContourDemo 
 * 
 * @version $Id: CustomContourDemo.java,v 1.1 2008-10-07 19:38:38 rstein Exp $
 */

public class Resonance extends JPanel {
	private static final long serialVersionUID = 7234889971218316473L;
    private Chart chart;
    ContourChartRenderer contourChartRenderer;
    
    // resonances
    private PolylineChartRenderer   resonanceLineRenderer;   
    private DefaultDataSource       resonances = new DefaultDataSource();
    private Style[]                 linestyles = new Style[1000];
    private static int              resonanceOrder = 4; 
    private int                     lastResonanceOrder = 0;
    
    // resonance line style definitions  
    private final Style rlineStyle1 = new Style(new BasicStroke(1.5f), new Color(255, 0, 0), new Color(255, 0, 0)); // first and second order resonances
    private final Style rlineStyle2 = new Style(new BasicStroke(1.5f, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND, 10, new float[] {5, 5}, 0), new Color(100, 100, 100), new Color(100, 100, 100)); // first and second order resonances
    private final Style rlineStyle3  = new Style(new BasicStroke(1.5f, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND, 10f, new float[] {5f,4f,1f,4f}, 10), new Color(100, 100, 100), new Color(100, 100, 100)); // third order resonances
    private final Style rlineStyle4  = new Style(new BasicStroke(1.5f, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND, 10f, new float[] {1f,8f,1f,8f}, 10), new Color(100, 100, 100), new Color(100, 100, 100)); // fourth order resonances
 
    // scan start and end parameter
    private double 			Qyst=2.95, Qyend=3.65, Qxst=3.95, Qxend=4.35;
    private double			dQy=-0.03, dQx=0.02;
    
    public Resonance() {
        initComponents();
    }

    private void initComponents() {
        this.setLayout(new GridBagLayout());
        this.setPreferredSize(new Dimension(500, 500));
        GridBagConstraints grid_c = new GridBagConstraints();
        grid_c.fill = GridBagConstraints.BOTH;
        grid_c.weightx = 1.0;
        grid_c.gridwidth = 2;
        grid_c.gridx = 0;

        grid_c.weighty = 100.0;
        grid_c.gridy = 0;
        add(createChart(), grid_c);
        
        //grid_c.gridwidth = 1;
        //grid_c.weighty = 1.0;
        //grid_c.gridy = 1;
        //add(createChart(), grid_c);
        //grid_c.gridx = 1;
        //add(new JPanel(), grid_c);
    }


    private Chart createChart() {
        this.chart = new Chart();

        // Set rendering type and data set 
        this.chart.getYAxis().setRange(Qyst, Qyend);	 
        this.chart.getXAxis().setRange(Qxst, Qxend);		           
        
        this.contourChartRenderer = new ContourChartRenderer();
        this.contourChartRenderer.setDataSet(createSpectrumDataSet(Qyst+dQy,Qyend+dQy,Qxst+dQx,Qxend+dQx,"HScan"));
        this.chart.addRenderer(this.contourChartRenderer); 
          
         
        drawResonanceLines(true, true, resonanceOrder,4,3);   // x: 0-0.5; y: 0-0.5
        drawResonanceLines(true, false, resonanceOrder,4,3);  // x: 0-0.5; y: 0.5-1
        drawResonanceLines(false, false, resonanceOrder,4,3); // x: 0.5-1; y: 0.5-1
        drawResonanceLines(false, true, resonanceOrder,4,3);  // x: 0.5-1; y: 0-0.5
        drawResonanceLines(false, false, resonanceOrder,3,3); // x: 0.5-1; y: 0.5-1
        drawResonanceLines(false, true, resonanceOrder,3,3);  // x: 0.5-1; y: 0-0.5
        drawResonanceLines(true, false, resonanceOrder,4,2);  // x: 0-0.5; y: 0.5-1
        drawResonanceLines(false, false, resonanceOrder,4,2); // x: 0.5-1; y: 0.5-1
        drawResonanceLines(false, false, resonanceOrder,3,2); // x: 0.5-1; y: 0.5-1
        
        //chart.setLegendVisible(true);
        //chart.setLegendTitle("Legend");
    
        // Set substeps for both scales
        this.chart.getXScale().setStepUnit(null, null);
        this.chart.getYScale().setStepUnit(null, null);
              
        // Set titles for both scales
        this.chart.setXScaleTitle("Q_x");	
        this.chart.setYScaleTitle("Q_y");
        
        //this.chart.addInteractor(new DataPickerInteractor()); 
              
        return this.chart;
    }
    
     /**
     * Read default spectra sample data o.k.
     * 
     * @return a data set
     */
    public DataSet createSpectrumDataSet(double Qyst, double Qyend, double Qxst, double Qxend, String filename) {
        try {
        	String file=filename+"_1.dat";
            BufferedReader in;
            in = new BufferedReader(
                    new InputStreamReader(ContourDemo.class.getResourceAsStream(file)));
            String buffer = in.readLine();

            String[] parse = buffer.split(" ");
            int ny = Integer.parseInt(parse[0]);
            int nx = Integer.parseInt(parse[1]);
            nx=nx*2;
            double[][] data = new double[nx][ny];
                  
            while ((buffer = in.readLine()) != null) {
                parse = buffer.split(" ");
                final int fft_bin = Integer.parseInt(parse[0]);
                final int row = Integer.parseInt(parse[1]);
                final Double val = Double.parseDouble(parse[2]);
                data[row][fft_bin] = val;//(float) indB(val);            
            }
            in.close();
            
            file=filename+"_2.dat";
            in = new BufferedReader(
                    new InputStreamReader(ContourDemo.class.getResourceAsStream(file)));
        
            buffer = in.readLine();
            parse = buffer.split(" ");
            
              
           
            while ((buffer = in.readLine()) != null) {
                parse = buffer.split(" ");
                final int fft_bin = Integer.parseInt(parse[0]);
                final int row = Integer.parseInt(parse[1]);
                final Double val = Double.parseDouble(parse[2]);
                data[row+nx/2][fft_bin] = val;//(float) indB(val);            
            }
            in.close();

            double[] y = new double[ny];
            for (int i = 0; i < ny ; i++) {
                y[i] = i*(Qyend-Qyst)/ny+Qyst;
            }

            double[] x = new double[nx];
            for (int i = 0; i < nx; i++) {
                x[i] = i* (Qxend-Qxst)/nx+Qxst;
            }

            // create and return data set.
            DefaultDataSet3D ds = new DefaultDataSet3D("Spectrum");
            ds.set(x,y,data, true, true); 
            return ds;
            

        } catch (final Exception e) {
            e.printStackTrace();
        }
        return null;
    }
      
    /**
     * draw tune resonance lines
     * ... is skipped for performance reasons if the resonance order hasn't change
     * 
     * @param Qxbelow
     * @param Qybelow
     * @param resonanceOrder
     */
    private void drawResonanceLines(boolean Qxbelow, boolean Qybelow, int resonanceOrder, int Qx, int Qy) {
            
        // 
        if (lastResonanceOrder == resonanceOrder) {
            return;
        }
        this.resonanceLineRenderer = new PolylineChartRenderer();
       
        
        // overwrites eventual older resonance lines
        this.resonances = new DefaultDataSource();

        // iterate through the various required orders
        // start with the largest to force the lower (more important) orders to draw over
        // the higher order ones.
        int rlineCount=0;
        for (int order = resonanceOrder; order>0; order--) {        
            
            final double resonances[][] = TuneResonanceDefinitionSmall.resonance;
            int count = 0;
            for (int i=0; i < resonances.length; i++) {
            	final double tuple[] = resonances[i];
                final int tupleOrder = (int)(tuple[0]); // first index
                if (tupleOrder==order) {
                    
                    // add line to list
                    DefaultDataSet resonanceLine = new DefaultDataSet("order="+order);
                    final double x0 = Qtrafo(tuple[1], Qxbelow);
                    final double y0 = Qtrafo(tuple[2], Qybelow);
                    final double x1 = Qtrafo(tuple[3], Qxbelow);
                    final double y1 = Qtrafo(tuple[4], Qybelow);
                    resonanceLine.add(x0+Qx, y0+Qy); // x0,y0 mirrored if necessary
                    resonanceLine.add(x1+Qx, y1+Qy); // x1,y1 mirrored if necessary                    
                    this.resonances.addDataSet(resonanceLine);
                    
                    switch (order) {
                    case 1: {
                        this.linestyles[rlineCount] = rlineStyle1.copy();  
                    } break;
                    case 2: {
                        this.linestyles[rlineCount] = rlineStyle2.copy();  
                    } break;
                    case 3: {
                        this.linestyles[rlineCount] = rlineStyle3.copy();   
                    } break;
                    case 4: {
                        this.linestyles[rlineCount] = rlineStyle4.copy();   
                    } 
                    }
                    count++;
                    rlineCount++;
                }
            }    
        }        
        this.resonanceLineRenderer.setStyles(this.linestyles);
        this.resonanceLineRenderer.setDataSource(this.resonances);
        this.chart.addRenderer(this.resonanceLineRenderer);
    }
    
    /**
     * Transforms value to above/below 0.5 
     * @param val     - the tune value
     * @param below05 - whether value is supposed to be 0.5 
     * @return
     */
     private double Qtrafo(double val, boolean below){
         if (below) return val;
         else return (1.0-val);
     }

 

    public static void main(String[] args) {
        JFrame f = new JFrame("Resonance");
        Resonance demo = new Resonance();
        f.getContentPane().add(demo);
        f.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        f.pack();
        f.setVisible(true);
    }
}