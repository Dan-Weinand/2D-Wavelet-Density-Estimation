/**
 * Encapsulates all the parameter variables used in the density estimator.
 * 
 * @author Gedeon Nyengele & Daniel Weinand.
 */
public class Settings {
	
	// Start level for both the scaling and wavelet functions.
	public static int startLevel       = 1;
	
	// Stop level for both the scaling and wavelet functions.
	public static int stopLevel        = 1;
	
	// Type of wavelet used.
	public static String waveletType   = "db6";
	
	// Types of wavelets supported.
	public static String[] waveletTypes =  {"coif1", "coif2", "coif3", "coif4", "coif5", "db2", "db3", "db4", "db5", "db6", "db7", "db8", "db9", "db10", "sym4", "sym5", "sym6", "sym7", "sym8", "sym9", "sym10", "dmey"};
	
	// Flag that determines which aging mechanism to use.
	//      0: no aging.
	//      1: Caudle aging method.
	//      2: Window method.
	public static int agingFlag        = 2;
	public static final int noAge      = 0;
	public static final int caudleAge  = 1;
	public static final int windowAge  = 2;
	
	// The distance between points in the density estimation plot
	public static double discretization = .01; 
	
	// How many samples to wait in between plot updates
	public static int updateFrequency = 100;
	
	// Size of the window for data aging using the window method.
	public static int windowSize       = 1400;
	
	// Caudle and Wegman's aging theta.
	public static double agingTheta    = .995;
	
	// Minimum and Maximum values on the domain of the density function.
	public static double[] densityRange = {-3.5, 3.5};
	
	public static double getMinimumRange() { return (densityRange[0]); }
	public static double getMaximumRange() { return (densityRange[1]); }
	
	// Minimum and Maximum values on the domain of the sample vector.
	public static double[] sampleRange = {0, 0};
	
	// Flag that determines whether or not the wavelet function should also be used.
	//      false: wavelet is OFF -> density approximation done with scaling function only.
	//      true : wavelet is ON  -> density approximation done with both scaling and wavelet functions.
	public static boolean waveletFlag      = false;
	
	// Path to the folder containing the look-up table for the wavelets.
	public static String waveletDataFolder = "";
	
	public static String dataFile = "../WaveletFiles/skewUni2Claw.csv";

}
