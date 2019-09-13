package com.mycompany.imagej;
import ij.*;
import ij.plugin.filter.BackgroundSubtracter;
import ij.plugin.Duplicator;
import ij.measure.Calibration;
import ij.process.*;
import ij.gui.*;
import ij.gui.Overlay;
import java.awt.*;
import ij.plugin.*;
import ij.plugin.filter.PlugInFilter;
import ij.plugin.frame.*;
import ij.measure.ResultsTable;
import java.util.*;
import java.awt.Point;
import java.awt.event.*;
import org.apache.commons.math3.stat.descriptive.moment.*;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import java.awt.image.BufferStrategy;
import java.awt.Rectangle;

public class Haralick_Features extends MouseAdapter implements PlugInFilter, KeyListener{
	ImagePlus img;
	ImagePlus hddisplay;
	ImageCanvas clickcanvas;
	Rectangle roi;
	ImageProcessor ip, roimask;
	boolean donesellectingpeaks = false, flag = true;
	int distanceno = 11, angleno = 45, binno = 8, pixelno = 0, graphsize = 400, pixeldistance;
	int[] bins, distancearr, peakvalleyarr;
	int maskheight, maskwidth;
	byte[] roixels;
	int[] pixarr;
	short[] imarr;
	int[][] validarr, xoffs, yoffs;
	double [][][][] glcmarr;
	double [][] tdarr;
	double [] radanglearr, degranglearr, hdarr, tthetaarr, realdistancearr;
	double score, realdistance = 2;
	int corrno = 10;
	PlotWindow shownplot;
	Plot shownhdplot;
	ArrayList<Double> clickxarray = new ArrayList<Double>(), clickyarray = new ArrayList<Double>();
	Roi roit;
	byte[] hdpixels;
	String[] options = {"Microtubules",  "Sarcomeres", "Microtubules-old"};
	String chosenoption = "Sarcomeres", units;
	Calibration icalibration;
	ResultsTable iresultstable = ResultsTable.getResultsTable();


	public int setup(String arg, ImagePlus img){
		this.img = img;
		IJ.register(Haralick_Features.class);
		return DOES_16+NO_CHANGES;
	}
	public int[][] glcmsub(int[] image, byte[] mask, int xoff, int yoff, int imwidth,
			int maskwidth, int maskheight, int xstart, int ystart) {
		int[][] retval = new int[binno][binno];
		for (int x = 0; x < maskwidth; x++) {
			int xoff1 = x + xoff;
			// int xoff2 = x - xoff;
			for (int y = 0; y < maskheight; y++) {

				int maski = y * maskwidth + x;
				if (mask[maski] != -1) {continue;}
				int imi = (y + ystart) * imwidth + xstart + x;
				int val = image[imi];


				int yoff1 = y + yoff;
				if (yoff1 < maskheight && yoff1 >= 0 && xoff1 < maskwidth && xoff1 >= 0) {
				    int maski1 = yoff1 * maskwidth + xoff1;
				    if (mask[maski1] == -1) {
						int imi1 = (yoff1 + ystart) * imwidth + xstart + xoff1;
						int val1 = image[imi1];
						retval[val][val1] += 1;
				    }
				}

				/* int yoff2 = y - yoff;
				if (yoff2 < maskheight && yoff2 >= 0 && xoff2 < maskwidth && xoff2 >= 0) {
					int maski2 = yoff2 * maskwidth + xoff2;
					if (mask[maski2] == -1) {
						int imi2 = (yoff2 + ystart) * imwidth + xstart + xoff2;
						int val2 = image[imi2];
						retval[val][val2] += 1;
					}
				} */
			}
		}
		return retval;
	}

	public void run(ImageProcessor ip) {
		this.img = new Duplicator().crop(img);
		this.ip = this.img.getProcessor();
		this.img.show();
		this.img.updateAndDraw();
		this.imarr = (short[]) this.ip.getPixels();
		this.roi = ip.getRoi();
		this.roit = new Roi(roi);
		this.roimask = ip.getMask();
		this.roixels = (byte[]) roimask.getPixels();
		this.maskheight = roi.height;
		this.maskwidth = roi.width;
		this.icalibration = img.getCalibration();
		this.units = icalibration.getUnits();
		if (units == "microns") {
			units = IJ.micronSymbol + "m";
		}
		GenericDialog gd1 = new GenericDialog("Select Algorithm Version and parameters");
		gd1.addChoice("Algorithm version", options, chosenoption);
	    gd1.addNumericField("Angles: ", angleno, 0);
	    gd1.addNumericField("Feature size in " + units + ": ", realdistance, 0);
	    gd1.addNumericField("numberofdistances", distanceno, 0);
	    gd1.addNumericField("Bin number", binno, 0);
		gd1.showDialog();
	    if (gd1.wasCanceled()) return;
	    this.angleno = (int) gd1.getNextNumber();
	    this.realdistance = gd1.getNextNumber();
	    this.pixeldistance = (int) icalibration.getRawX(realdistance);
	    this.distanceno = (int) gd1.getNextNumber();
	    this.binno = (int) gd1.getNextNumber();
		this.chosenoption =  ((Choice) gd1.getChoices().get(0)).getSelectedItem();
	    if (chosenoption == "Sarcomeres"|| chosenoption == "Microtubules") {
	    	this.distanceno = pixeldistance * 3;
	    }
		if(chosenoption == "Sarcomeres") {
			int backgroundradius = pixeldistance / 2;
			System.out.println("radius: " + Integer.toString(backgroundradius));
			new BackgroundSubtracter().rollingBallBackground(this.ip, backgroundradius, false, false, false, false, false);
		}
	    this.distancearr = new int[distanceno];
		this.radanglearr = new double[angleno];
		this.degranglearr = new double[angleno];
		this.validarr = new int[angleno][distanceno];
		int di = pixeldistance;
		int db = 1;
		if (chosenoption == "Sarcomeres" || chosenoption == "Microtubules" ) {
			for (int i = 0; i < pixeldistance * 3; i++){
				distancearr[i] = i+1;
			}
		}
		else if (chosenoption == "Microtubules-old" ) {
			for (int i = 0; i < distanceno; i++){
				distancearr[i] = di;
				di += db;
				db ++;
			}
		}
		this.realdistancearr = new double[distanceno];
		for (int i = 0; i < distanceno; i++ ) {
			realdistancearr[i] = icalibration.getX((double) distancearr[i]);
		}
		for (int i = 0; i < angleno; i++){
			radanglearr[i] = (Math.PI/(angleno)) * i;
			degranglearr[i] = radanglearr[i] * (180 / Math.PI);
		}
		this.xoffs = calculatexoffs(radanglearr, distancearr);
		this.yoffs = calculateyoffs(radanglearr, distancearr);
		this.bins = histogrambins(this.ip, roi.width, roi.height);
		this.pixarr = pixbin(bins, this.ip, roi.width, roi.height);
		tdarr = tdarrwbins();
		ImagePlus pixarrdisplay = NewImage.createShortImage("pixarr", roi.width, roi.height, 1, NewImage.FILL_BLACK);
		ImageProcessor pixarrprocessor = pixarrdisplay.getProcessor();
		short[] pixarrpixels = (short[]) pixarrprocessor.getPixels();
		for (int x = 0; x<roi.width; x++){
			for (int y = 0; y<roi.height; y++){
				pixarrpixels[y*roi.width + x] = (short) pixarr[y * maskwidth + x];
			}
		}
		pixarrdisplay.show();
		pixarrdisplay.updateAndDraw();
		if (chosenoption == "Sarcomeres") {
			shownplot = showsarccor(tdarr);
			PlotCanvas sarcplotcanvas = (PlotCanvas) shownplot.getCanvas();
			sarcplotcanvas.addMouseListener(this);
			sarcplotcanvas.addKeyListener(this);
			shownhdplot = shownplot.getPlot();
		}
		else if (chosenoption == "Microtubules") {
			shownplot = shownewtubcor(tdarr);
			PlotCanvas newtubplotcanvas = (PlotCanvas) shownplot.getCanvas();
			newtubplotcanvas.addMouseListener(this);
			newtubplotcanvas.addKeyListener(this);
			shownhdplot = shownplot.getPlot();
		}
		else if (chosenoption == "Microtubules-old") {
			this.tthetaarr = ttheta(tdarr);
			hdarr = hd(tthetaarr);
			shownplot = showtubcor(hdarr);
			PlotCanvas tubplotcanvas = (PlotCanvas) shownplot.getCanvas();
			tubplotcanvas.addMouseListener(this);
			tubplotcanvas.addKeyListener(this);
			shownhdplot = shownplot.getPlot();
		}
	}
//This code forms the math of the program
	public int[][] calculatexoffs(double[]angles, int[] distances) {
		int[][] retval = new int[angleno][distanceno];
		for (int anglei=0; anglei<angleno; anglei++){
			for (int distancei=0; distancei<distanceno; distancei++){
				retval[anglei][distancei] = (int) Math.round(Math.cos(radanglearr[anglei]) * distancearr[distancei]);
			}
		}
		return retval;
	}
	public int[][] calculateyoffs(double[] angles, int[] distances) {
		int[][] retval = new int[angleno][distanceno];
		for (int anglei=0; anglei<angleno; anglei++){
			for (int distancei=0; distancei<distanceno; distancei++){
				retval[anglei][distancei] = (int) Math.round(Math.sin(radanglearr[anglei]) * distancearr[distancei]);
			}
		}
		return retval;
	}
	public double[][][][] glcm(){
		glcmarr = new double[angleno][distanceno][binno][binno];
		//calculate distances and angles to produce glcm for
		int xoff, yoff;
		for (int anglei=0; anglei<angleno; anglei++){
			for (int distancei=0; distancei<distanceno; distancei++){
				xoff = xoffs[anglei][distancei];
				yoff = yoffs[anglei][distancei];
				int[][] glcmsubarr = glcmsub(pixarr, roixels, xoff, yoff, maskwidth,
						maskwidth, maskheight, 0, 0);
				int pairno = 0;
				for (int gx = 0; gx < glcmsubarr.length; gx++) {
					for (int gy = 0; gy < glcmsubarr[0].length; gy++) {
						int entry = glcmsubarr[gx][gy];
						glcmarr[anglei][distancei][gx][gy] = entry;
						pairno += entry;
					}
				}
				if (pairno == 0){
					validarr[anglei][distancei] = 0;
				}
				else {
					validarr[anglei][distancei] = 1;
					continue;
				}
				for(int vi = 0; vi<binno; vi++ ) {
					for (int vj = 0; vj<binno; vj++) {
						glcmarr[anglei][distancei][vi][vj] = glcmarr[anglei][distancei][vi][vj] / pairno;
					}
				}
			}
		}
	return glcmarr;
	}
	public double[][] tdarrnobins(){
		double[][] retval = new double[angleno][distanceno];
		for (int anglei = 0; anglei < angleno; anglei++) {
			for (int distancei = 0; distancei < distanceno; distancei++) {
				int xoff = xoffs[anglei][distancei];
				int yoff = yoffs[anglei][distancei];
				retval[anglei][distancei] = glcmsubfull(imarr, roixels, xoff, yoff, maskwidth,
						maskheight,anglei, distancei);
			}
		}
		return retval;
	}
	public double[][] tdarrwbins(){
		double[][] retval = new double[angleno][distanceno];
		short[] lpixarr = new short[roi.width * roi.height];
		for (int i = 0; i < roi.width*roi.height; i++) {
			lpixarr[i] = (short) pixarr[i];
		}
		for (int anglei = 0; anglei < angleno; anglei++) {
			for (int distancei = 0; distancei < distanceno; distancei++) {
				int xoff = xoffs[anglei][distancei];
				int yoff = yoffs[anglei][distancei];
				retval[anglei][distancei] = glcmsubfull(lpixarr, roixels, xoff, yoff, maskwidth,
						maskheight, anglei, distancei);
			}
		}
		return retval;
	}
	private int argmax(double[] arr) {
		int len = arr.length;
		double compval = arr[0];
		int retval = 0;
		for (int i = 0; i < len; i++) {
			if (arr[i] > compval) {
				compval = arr[i];
				retval = i;
			}
		}
		return retval;
	}
	public double glcmsubfull(short[] image, byte[] mask, int xoff, int yoff,
			int maskwidth, int maskheight, int anglei, int distancei) {
		ArrayList<Integer> xarray = new ArrayList<Integer>();
		ArrayList<Integer> yarray = new ArrayList<Integer>();
		for (int x = 0; x < maskwidth; x++) {
			int xoff1 = x + xoff;
			// int xoff2 = x - xoff;
			for (int y = 0; y < maskheight; y++) {

				int maski = y * maskwidth + x;
				if (mask[maski] != -1) {continue;}
				int imi = (y) * maskwidth + x;
				int xval = (int) image[imi];


				int yoff1 = y + yoff;
				if (yoff1 < maskheight && yoff1 >= 0 && xoff1 < maskwidth && xoff1 >= 0) {
				    int maski1 = yoff1 * maskwidth + xoff1;
				    if (mask[maski1] == -1) {
						int imi1 = (yoff1) * maskwidth + xoff1;
						int yval = (int) image[imi1];
						xarray.add(xval);
						yarray.add(yval);
				    }
				}
			}
		}
		double[] xarrayy = new double[xarray.size()];
		for (int i = 0; i < xarray.size(); i++) {
			xarrayy[i] = (double) xarray.get(i);
		}
		double[] yarrayy = new double[yarray.size()];
		for (int i = 0; i < yarray.size(); i++) {
			yarrayy[i] = (double) yarray.get(i);
		}
		double retval;
		try {
			retval = new PearsonsCorrelation().correlation(xarrayy, yarrayy);
			validarr[anglei][distancei] = 1;
		}
		catch (Exception e) {
			validarr[anglei][distancei] = 0;
			retval = 0;
		}
		if (Double.isNaN(retval)) {
			retval = 0;
			validarr[anglei][distancei] = 0;
		}
		if (chosenoption == "Sarcomeres")
			return (double) retval;
		else if (chosenoption == "Microtubules")
			return Math.abs(retval);
		return 0;

	}
	private int[] histogrambins(ImageProcessor ip, int maskwidth, int maskheight){
		int val;
		int min = 65535;
		int max = 0;
		//find minima and maxima
		for (int x = 0; x < maskwidth; x++){ for (int y = 0; y < maskheight; y++){
			if (roimask.get(x, y) == 255){
			val = ip.get(x, y);
			pixelno++;
				if (val < min){
					min = val;
				}
				if (val > max){
					max = val;
				}}
		}}
		int greylevels = max-min+1;
		int[] histogram = new int[greylevels];
		for (int x = 0; x < maskwidth; x++){ for (int y = 0; y < maskheight; y++){{
			if (roimask.get(x, y) == 255){
			val = ip.get(x, y);
			histogram[val-min]++;
		  }
		}}}
		int binthresh = pixelno / binno;
		int histsum = 0;
		int[] bins = new int [binno];
		int countbins = 0;
		for (int i=0; i<greylevels; i++){
			histsum += histogram[i];
			if (histsum > binthresh * countbins){
				bins[countbins] = i+min;
				countbins++;
				if (countbins >= binno){
					break;
				}
			}
		}
		return bins;
	}
	private int[] pixbin(int[] bins, ImageProcessor ip, int maskwidth, int maskheight){
	    int[] pixarr = new int[roi.width * roi.height];
		int val;
		for (int x = 0; x < maskwidth; x++){ for (int y = 0; y < maskheight; y++){{
			val =  ip.get(x, y);
			for (int i= (binno-1); i>=0; i--){
				if (val >= bins[i]){
					pixarr[x + (y) * maskwidth] = i;
					break;
				}
			}
		}}}
	return pixarr;}
	private double[] px(double[][] arr) {
		double[] retval = new double[binno];
		for (int vi = 0; vi < binno; vi++) {
			for (int vj = 0; vj < binno; vj ++) {
				retval[vi] += vi * arr[vi][vj];
			}
		}
		return retval;
	}
	private double[] py(double[][] arr) {
		double[] retval = new double[binno];
		for (int vi = 0; vi < binno; vi++) {
			for (int vj = 0; vj < binno; vj ++) {
				retval[vi] += vj * arr[vi][vj];
			}
		}
		return retval;
	}
	private double mu(double[] arr) {
		return new Mean().evaluate(arr);
	}
	private double sd(double[] arr) {
		return new StandardDeviation().evaluate(arr);
	}

	public double[][] td(double[][][][] arr){
		double[][] retval = new double[angleno][distanceno];
		for(int anglei = 0; anglei < angleno; anglei++) {
			for(int distancei = 0; distancei < distanceno; distancei++) {
				if (validarr[anglei][distancei] == 1) {
					double [][] subarr = arr[anglei][distancei];
					double[] pxarr = px(subarr);
					double[] pyarr = py(subarr);
					double mux = mu(pxarr);
					double muy = mu(pyarr);
					double sdx = sd(pxarr);
					double sdy = sd(pyarr);
					for (int vi = 0; vi < binno; vi ++) {
						for (int vj = 0; vj < binno; vj ++) {
							double vij = subarr[vi][vj];
							retval[anglei][distancei] += vij * (vi-mux) * (vj - muy) / (sdx * sdy);
						}
					}
				}
			}
		}
		return retval;
	}
	public double[] ttheta(double[][] tdval){
		System.out.println("validarr: " + Arrays.deepToString(validarr));
		double[] retval = new double[angleno];
		for (int anglei = 0; anglei < angleno; anglei++){
			double[] subarr = tdval[anglei];
			int[] subvalidarr = validarr[anglei];
			double subretval = 0;
			int subvalidno = 0;
			for (int distancei = 0; distancei < distanceno; distancei++){
				if(subvalidarr[distancei] == 1) {
					subretval += subarr[distancei];
					subvalidno ++;
				}
			}
			retval[anglei] = (double) subretval / (double) subvalidno;
		}
		return retval;
	}
	public double[] hd(double[] tthetaarr) {
		double[] retval = new double[angleno];
		double tthetasum = 0;
		for(int anglei = 0; anglei < angleno; anglei++) {
			tthetasum += tthetaarr[anglei];
		}
		for(int anglei = 0; anglei < angleno; anglei++) {
			retval[anglei] = tthetaarr[anglei] / tthetasum;
		}
		return retval;
	}
	public double max(double[] arr) {
		double retval = arr[0];
		for (int i = 0; i < arr.length; i++){
			if (arr[i] > retval) {
				retval = arr[i];
			}
		}
		return retval;
	}
	public double min(double[] arr) {
		double retval = arr[0];
		for (int i = 0; i < arr.length; i++){
			if (arr[i] < retval) {
				retval = arr[i];
			}
		}
		return retval;
	}
	public PlotWindow showsarccor(double[][] tdarr){
		Plot valplot = new Plot("Angle-Distance Correlation Values",  "Distance (" + units + ")", "Correlation" );
		double jump = (double) 1 / (angleno);
		for (int anglei = 0; anglei < angleno; anglei++) {
			double degreedistanceanglecor[] = new double[distanceno];
			double[] tdsubarr = tdarr[anglei];
			Color thiscolor = Color.getHSBColor((float) (jump*anglei), 1.0f, 1.0f);
			for (int distancei = 0; distancei < distanceno; distancei++) {
				degreedistanceanglecor[distancei] = tdsubarr[distancei];
			}
			valplot.setColor(thiscolor);
			valplot.add("line", (double[]) realdistancearr, tdsubarr);
			valplot.show();
			valplot.update();
		}
		return valplot.show();
	}
	public PlotWindow shownewtubcor(double[][]tdarr) {
		Plot valplot = new Plot("Angle-Distance Correlation Values", "Angle",  "Distance (" + units + ")");
		double jump = (double) 1 / (corrno);
		for (int corri = 0; corri < corrno; corri++) {
			double corrval = 0.1*corri;
			double[] tdsubcutarr = new double[angleno];
			Arrays.fill(tdsubcutarr, realdistance * 3);
			tdsubcutarr[0] = 0;
			for (int anglei = 0; anglei < angleno; anglei++) {
				for (int distancei = 0; distancei < distanceno; distancei ++) {
					if (tdarr[anglei][distancei] < corrval) {
						tdsubcutarr[anglei] = distancearr[distancei];
						break;
					}
				}
			}
			System.out.println(Arrays.toString(tdsubcutarr));
			Color thiscolor = Color.getHSBColor((float) (jump*corri), 1.0f, 1.0f);
			valplot.setColor(thiscolor);
			valplot.add("line", (double[]) degranglearr, tdsubcutarr);
			valplot.redrawGrid();
			valplot.show();
			valplot.update();
		}
		return valplot.show();
	}
	public PlotWindow showtubcor(double[] hdarr) {
		Plot valplot = new Plot("Average Angle Correlation Values", "Angle", "Correlation");
		valplot.add("line", degranglearr, hdarr);
		PlotWindow valplotwindow = valplot.show();
		return valplotwindow;
	}
	public void onclick() {

	}
	public void mousePressed(MouseEvent e) {
		int x = e.getX();
		int y = e.getY();
		Rectangle wheretheplotis = shownhdplot.getDrawingFrame();
		x = x - wheretheplotis.x;
		y = y - wheretheplotis.y;
		int plotwidth = wheretheplotis.width;
		int plotheight = wheretheplotis.height;
		double[] limits = shownhdplot.getLimits();
		double xmin = limits[0], xmax = limits[1], ymin = limits[2], ymax = limits[3];
		double xlen = xmax - xmin;
		double ylen = ymax - ymin;
		double truexcoord = ((double) x / plotwidth) * xlen + xmin;
		double trueycoord = ymax - ((double) (y) / plotheight) * ylen;
		System.out.println("x,y" + Double.toString(truexcoord) + " , " + Double.toString(trueycoord));
		clickxarray.add(truexcoord);
		clickyarray.add(trueycoord);
		if (clickxarray.size() % 2 == 1) {
			shownhdplot.addText("Peak", truexcoord, trueycoord);
		}
		else {
			shownhdplot.addText("Valley", truexcoord, trueycoord);
		}
	}
	public void keyPressed(KeyEvent e) {
		if (e.getKeyCode() == KeyEvent.VK_ENTER) {
			double[] xarrayy = new double[clickxarray.size()], yarrayy = new double[clickyarray.size()];
			for (int i = 0; i < clickxarray.size(); i++) {
				xarrayy[i] =  clickxarray.get(i);
				yarrayy[i] =  clickyarray.get(i);
			}
			if (chosenoption == "Microtubules") {
				int[] peakvalleyarr = new int[xarrayy.length];
				for (int peakvalleyi = 0; peakvalleyi < xarrayy.length; peakvalleyi++) {
					double peakvalleytheta = xarrayy[peakvalleyi];
					for (int anglei = 0; anglei < angleno; anglei++) {
						if (radanglearr[circleget(anglei, angleno)] >= peakvalleytheta && radanglearr[circleget(anglei + 1,angleno)] < peakvalleytheta) {
							peakvalleyarr[peakvalleyi] = anglei;
						}
					}
				}
				this.peakvalleyarr = peakvalleyarr;
				double dscore = dscorelet(hdarr, radanglearr, peakvalleyarr);
				iresultstable.incrementCounter();
				iresultstable.addValue("Microtubule Organization Score", (double) dscore);
				iresultstable.show("Haralick Measurements");
			}
			else if (chosenoption == "Sarcomeres") {
				double score = yarrayy[0] - yarrayy[1];
				double sarcomerelength = xarrayy[0];
				iresultstable.incrementCounter();
				iresultstable.addValue("Sarcomere Organization Score", (double) score);
				iresultstable.addValue("Sarcomere Length (" + units + ")", (double) sarcomerelength);
				iresultstable.show("Haralick Measurements");
			}
		}
	}
	public int circleget(int index, int length) {
		return Math.floorMod(index,  length);
	}
	public void keyReleased(KeyEvent e) {}
	public void keyTyped(KeyEvent e) {}

	public double interpolate(double xr, double x1, double x2, double y1, double y2) {
		return y1 + ((xr - x1) / (x2 - x1)) * (y2 - y1);
	}

	public int coordtransform(double x, double xdomain, int xsize) {
		return (int) Math.round((x/xdomain)*xsize);
	}
	public double dscorelet(double[] hdarr, double[] anglearr, int[] peakvalleyarr) {
		double retval = 0;
		double gammaval = 0;
		for(int peaki = 0; peaki < peakvalleyarr.length; peaki += 2) {
			int valley1i = peakvalleyarr[circleget(peaki - 1, peakvalleyarr.length)];
			if (valley1i > peaki) {
				valley1i = - (angleno - valley1i);
			}
			int valley2i = peakvalleyarr[circleget(peaki + 1, peakvalleyarr.length)];
			if (valley2i < peaki) {
				valley2i = - (angleno + valley2i);
			}
			double peakrad = anglearr[peaki];
			double uppeakrad = peakrad + Math.PI;
			double downpeakrad = peakrad - Math.PI;
			for (int anglei = valley1i; anglei < angleno; anglei ++) {
				double thisangle = anglearr[circleget(anglei, angleno)];
				double thishd = hdarr[circleget(anglei, angleno)];
				double d1 = Math.abs(peakrad - thisangle);
				double d2 = Math.abs(uppeakrad - thisangle);
				double d3 = Math.abs(downpeakrad - thisangle);
				double finald = min(new double[] {d1, d2, d3});
				gammaval += finald*finald;
				retval += thishd*finald*finald;
			}
		}
		gammaval = angleno / gammaval;
		return (gammaval * retval);
	}

	private void test() {
		double epsilon = 0.00001;
		double mutest = mu(new double[] {1, 1.2, 1.4, 1.6});
		double mucorrect = 1.4;
		double sdtest = sd(new double[] {1, 1.2, 1.4, 1.6});
		double sdcorrect = 0.3162278;
		assert Math.abs(mutest - 1.3) < epsilon: "mu function incorrect" + Double.toString(mutest);
		assert Math.abs(sdtest - 1.3) < epsilon: "sd function incorrect" + Double.toString(sdtest);
		int[] testglcmroi = {0, 0, 4, 5, 7,
							  1, 2, 4, 6, 0,
							  3, 4, 6, 0, 1,
							  7, 4, 0, 1, 4};
		byte[] testglcmmask = {-1, -1, -1, -1, -1,
				  			   -1, -1, -1, -1, -1,
				  			   -1, -1, -1, -1, -1,
				  			   -1, -1, -1, -1, -1};

		byte[] testglcmmask2 = {-1, -1, -1, -1, -1,
	  			   			   	-1, -1, -1, -1, -1,
	  			   			   	-1, -1,  0, -1, -1,
	  			   			   	-1, -1, -1, -1, -1};

		int[][] testglcmcorrect = {{1, 2, 0, 0, 1, 0, 0, 0},
								   {0, 0, 1, 0, 1, 0, 0, 0},
								   {0, 0, 0, 0, 1, 0, 0, 0},
								   {0, 0, 0, 0, 1, 0, 0, 0},
								   {1, 0, 0, 0, 0, 1, 2, 0},
								   {0, 0, 0, 0, 0, 0, 0, 1},
								   {2, 0, 0, 0, 0, 0, 0, 0},
								   {0, 0, 0, 0, 1, 0, 0, 0}};

		int[][] testglcmcorrect2 = {{1, 2, 0, 0, 1, 0, 0, 0},
								   {0, 0, 1, 0, 1, 0, 0, 0},
								   {0, 0, 0, 0, 1, 0, 0, 0},
								   {0, 0, 0, 0, 1, 0, 0, 0},
								   {1, 0, 0, 0, 0, 1, 1, 0},
								   {0, 0, 0, 0, 0, 0, 0, 1},
								   {1, 0, 0, 0, 0, 0, 0, 0},
								   {0, 0, 0, 0, 1, 0, 0, 0}};
		int testimwidth = 5;
		int testmaskwidth = 5;
		int testmaskheight = 4;
		int testxoff = 1;
		int testyoff = 0, testxstart = 0, testystart = 0;
		int[][] testglcmtest = glcmsub(testglcmroi, testglcmmask, testxoff, testyoff, testimwidth,
				testmaskwidth, testmaskheight, testxstart, testystart);
		int[][] testglcmtest2 = glcmsub(testglcmroi, testglcmmask2, testxoff, testyoff, testimwidth,
				testmaskwidth, testmaskheight, testxstart, testystart);
		System.out.println(Arrays.deepToString(testglcmtest));
		System.out.println(Arrays.deepToString(testglcmtest2));
		assert Arrays.equals(testglcmcorrect, testglcmtest) : "glcmsub function incorrect";
		assert Arrays.equals(testglcmcorrect2, testglcmtest2) : "glcmsub function incorrect, check masking";
	}

}
