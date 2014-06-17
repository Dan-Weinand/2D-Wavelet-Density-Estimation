/**
 *  This class represents the density model.
 *  @author Gedeon Nyengele & Daniel Weinand.
 *  
 */

import org.sf.surfaceplot.ISurfacePlotModel;



public class DensityModel implements ISurfacePlotModel{
	
	// 2-dimensional array to store the density function.
	public static double[][] density;
	
	// Discretization value used in density calculation.
	private float discretization;
	
	// Domain of the density function.
	private double[] densityDomain;
	
	// Maximum value of the density function.
	public static float maximumDensity = 0.0f;
	
	// Axis Labels.
	private String xLabel, yLabel, zLabel;
	
	/**
	 * Constructor.
	 * @param discretization : discretization value used in density calculation.
	 * @param domain : domain of the density function.
	 */
	public DensityModel(double discretization, double[] domain){
		this.discretization = (float) discretization;
		this.densityDomain  = domain;
		this.xlabel( "X" );
		this.ylabel( "Y" );
		this.zlabel( "Z" );
	} // end constructor.
	
	/**
	 * Retrieves the value of the density at a given pt.
	 * @param x : position on the x-axis.
	 * @param y : position on the y-axis.
	 */
	public float calculateZ(float x, float y){
		
		// Return meaningless value if the point does not lie within
		// the domain of the density function.
		if(!isWithinDensityDomain(x, y)) 
		{
			return Float.MIN_VALUE;
		}
		
		// Get the index positions of the point in the density table.
		int xPos = (int) ((x - densityDomain[0]) / discretization );
		int yPos = (int) ((y - densityDomain[0]) / discretization );
		
		// Retrieve the value of the density at the point.
		float densityHere = (float)  density[xPos][yPos];
		
		return densityHere;
	} // end method calculateZ(float x, float y)
	
	public int getPlotMode(){
		return ISurfacePlotModel.PLOT_MODE_SPECTRUM;
	} // end method getPlotMode()
	
	public boolean isBoxed(){
		return true;
	} // end method isBoxed()
	
	public boolean isMesh(){
		return true;
	} // end method isMesh()
	
	public boolean isScaleBox(){
		return false;
	} // end method isScaleBox()
	
	public boolean isDisplayXY(){
		return true;
	} // end method isDisplayXY()
	
	public boolean isDisplayZ(){
		return true;
	} // end method isDisplayZ()
	
	public boolean isDisplayGrids(){
		return true;
	} // end method isDisplayGrids()
	
	public int getCalcDivisions(){
		return (int) ( (densityDomain[1] - densityDomain[0]) / discretization + 1);
		//return 50;
	} // end getCalcDivisions()
	
	public int getDispDivisions(){
		return (int) ( (densityDomain[1] - densityDomain[0]) / discretization + 1);
	} // end method getDispDivisions()
	
	public float getXMin(){
		return (float) densityDomain[0];
	} // end method getXMin()
	
	public float getXMax(){
		return (float) densityDomain[1];
	} // end method getXMax()
	
	public float getYMin(){
		return (float) densityDomain[0];
	} // end method getYMin()
	
	public float getYMax(){
		return (float) densityDomain[1];
	} // end method getYMax()
	
	public float getZMin(){
		return 0.0f;
	} // end method getZMin()
	
	public float getZMax(){
		return maximumDensity;
	} // end method getZMax()
	
	public String getXAxisLabel(){
		return xLabel;
	} // end method getXAxisLabel()
	
	public String getYAxisLabel(){
		return yLabel;
	} // end method getYAxisLabel()
	
	public String getZAxisLabel(){
		return zLabel;
	} // end method getZAxisLabel()
	
	public void xlabel( String label ){
		xLabel = label;
	} // end method xlabel( String label )
	
	public void ylabel( String label ){
		yLabel = label;
	} // end method ylabel( String label )
	
	public void zlabel( String label ){
		zLabel = label;
	} // end method zlabel( String label )
	/**
	 * Resets the old density function with the new one.
	 * @param newDensity : new density function approximation.
	 * @param maxDensity : maximum value of the density function;
	 */
	public static void updateDensity(double[][] newDensity, float maxDensity ){
		density        = newDensity;
		maximumDensity = maxDensity;
	} // end method updateDensity(float[][] newDensity )
	
	/**
	 * Determines whether the given point lies within the domain
	 * of the density or not.
	 * @param x : x-coordinate of the point.
	 * @param y : y-coordinate of the point.
	 * @return
	 */
	private boolean isWithinDensityDomain( float x , float y){
	
		if(x < (float)densityDomain[0] || x > (float)densityDomain[1]) return false;
		if(y < (float)densityDomain[0] || y > (float)densityDomain[1]) return false;
		
		return true;
	} // end method isWithinDensityDomain( float x , float y)
	
	public static double[][] getDensity(){
		return density;
	}
}
