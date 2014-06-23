import java.util.ArrayList;
/**
 * Contains the coefficients and translates of both the scaling and wavelet functions.
 * 
 * @author Gedeon Nyengele & Daniel Weinand
 * 
 */
public class Transform {
	
	// Scaling function's coefficients.
	public static double[][] scalingCoefficients;
	
	// wavelet function's coefficients.
	public static ArrayList<double[][]> psiPsiCoefficients;
	public static ArrayList<double[][]> psiPhiCoefficients;
	public static ArrayList<double[][]> phiPsiCoefficients;
	
	// Scaling function's translates.
	public static ArrayList<Double> scalingTranslates;
		
	// wavelet function's translates.
	public static ArrayList<ArrayList<Double>> waveletTranslates;
	
	
	
	
	

} // end class Tranform.
