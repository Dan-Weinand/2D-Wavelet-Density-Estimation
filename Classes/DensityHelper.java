/**
 * Helper functions for the density estimation algorithm.
 * 
 * @author Daniel Weinand & Gedeon Nyengele
 * 
 */
import java.util.ArrayList;
import java.util.Collections;

import de.erichseifert.gral.data.DataTable;


public class DensityHelper {
	
	private static double[][] oldSamples;		// The old samples in the window
	private static int N;					// How many samples have been read in
	
	/**
	 * Checks that the sample point X is within the domain of the density function.
	 * @param X : the data point to check
	 * @return  : whether or not the point is in the domain
	 */
	public static boolean inRange (double X) {
		return (X >= Settings.getMinimumRange() && X <= Settings.getMaximumRange());
	} // end inRange.
	
	/**
	 * Updates the function coefficients based on the incoming data point and
	 * the data point leaving the sliding window.
	 * 
	 * Post: the coefficients are updated as needed
	 * 
	 * @param Xnew : the new data point to update the coefficients based on
	 */
	public static void updateCoefficients(double[] Xnew){
		
		if (Settings.waveletFlag) {
		//	updateWaveletCoefficients(Xnew);
		}
		
		// The normalizing constant for the scaling basis functions
		double scaleNormalizer = Math.pow(2, Settings.startLevel/2.0);
		double ageNorm = 1;
		if (Settings.agingFlag == Settings.windowAge) {
			ageNorm /= Settings.windowSize;
		}
		else if (Settings.agingFlag == Settings.caudleAge) {
			ageNorm *= (1 - Settings.agingTheta);
		}
		else if (Settings.agingFlag == Settings.noAge) {
			ageNorm /= (N+1)*1.0;
		}
		
/*		// Scale coefficients if Caudle aging is being used
		if (Settings.agingFlag == Settings.caudleAge) {
			for (int x1Index = 0; 
					x1Index < Transform.scalingCoefficients.length;
					x1Index++) {
				for (int x2Index = 0; 
						x2Index < Transform.scalingCoefficients.length;
						x2Index++) {
					Transform.scalingCoefficients[x1Index][x2Index] *= Settings.agingTheta;
				}
			}
				
		}*/
		
		// Recursively compute coefficients if no aging is used
		else if (Settings.agingFlag == Settings.noAge){
			for (int x1Index = 0; 
					x1Index < Transform.scalingCoefficients.length;
					x1Index++) {
				for (int x2Index = 0; 
						x2Index < Transform.scalingCoefficients.length;
						x2Index++) {
					Transform.scalingCoefficients[x1Index][x2Index] *= (N)/(N+1.0);
				}
			}
		}
		
		// Subtract old samples effect if window aging is used
		else if (Settings.agingFlag == Settings.windowAge) {
			
			// Only remove a sample if there have been more than window size samples
			if (N > Settings.windowSize){ 
				double[] Xold = oldSamples[N % Settings.windowSize];
				double x1 = Xold[0];
				double x2 = Xold[1];
			
			
				//Loop through the x1 translations for the scaling basis functions
				int x1Ind = 0;
				for (double k1 : Transform.scalingTranslates) {
				
					// Get the translated & scaled data point
					double x1Scaled = Math.pow(2, Settings.startLevel) * x1 - k1;
				
					// If the x1 wavelet supports the data point, check through x2 supports
					if (Wavelet.inSupport(x1Scaled)) {
						double phi1Here = scaleNormalizer*Wavelet.getPhiAt(x1Scaled);
						int x2Ind = 0;
						
						//Loop through the x2 translations for the scaling basis functions
						for (double k2 : Transform.scalingTranslates) {
							
							// Get the translated & scaled data point
							double x2Scaled = Math.pow(2, Settings.startLevel) * x2 - k2;
							
							// Update appropriate coefficient iff in x2 support
							if (Wavelet.inSupport(x2Scaled)) {
								double phi2Here = scaleNormalizer*Wavelet.getPhiAt(x1Scaled);
								double coeffSub = phi1Here*phi2Here/Settings.windowSize;
								
								Transform.scalingCoefficients[x1Ind][x2Ind] -= coeffSub;
							} // end updating
							x2Ind++;
						} // end x2 looping										
					} // end in x1 support
					x1Ind++;
				} // end x1 looping
			}
			
			oldSamples[N % Settings.windowSize] = Xnew;
		}
		
		// The x1 and x2 coordinates of the incoming point
		double x1 = Xnew[0];
		double x2 = Xnew[1];
	
	
		//Loop through the x1 translations for the scaling basis functions
		int x1Ind = 0;
		for (double k1 : Transform.scalingTranslates) {
		
			// Get the translated & scaled data point
			double x1Scaled = Math.pow(2, Settings.startLevel) * x1 - k1;
		
			// If the x1 wavelet supports the data point, check through x2 supports
			if (Wavelet.inSupport(x1Scaled)) {
				double phi1Here = scaleNormalizer*Wavelet.getPhiAt(x1Scaled);
				int x2Ind = 0;
				
				//Loop through the x2 translations for the scaling basis functions
				for (double k2 : Transform.scalingTranslates) {
					
					// Get the translated & scaled data point
					double x2Scaled = Math.pow(2, Settings.startLevel) * x2 - k2;
					
					// Update appropriate coefficient iff in x2 support
					if (Wavelet.inSupport(x2Scaled)) {
						double phi2Here = scaleNormalizer*Wavelet.getPhiAt(x1Scaled);
						double coeffAdd = phi1Here*phi2Here*ageNorm;
						
						Transform.scalingCoefficients[x1Ind][x2Ind] += coeffAdd;
					} // end updating
					x2Ind++;
				} // end x2 looping										
			} // end in x1 support
			x1Ind++;
		} // end x1 looping
		
		N++;
	} // end updateCoefficients
	
	/**
	 * Updates the wavelet coefficients based on the incoming data point and
	 * the data point leaving the sliding window.
	 * 
	 * Post: the coefficients are updated as needed
	 * 
	 * @param Xnew : the new data point to update the coefficients based on
	 */
/*	private static void updateWaveletCoefficients(double[] Xnew) {
		
		// Short hand for start level
		int j0 = Settings.startLevel;
		
		// Loop through for each resolution level
		for (int j = Settings.startLevel; j <= Settings.stopLevel; j++) {
			
			// The normalizing constant for the wavelet basis functions
			double waveNormalizer = Math.pow(2, j/2.0);
			if (Settings.agingFlag == Settings.windowAge) {
				waveNormalizer /= Settings.windowSize;
			}
			else if (Settings.agingFlag == Settings.caudleAge) {
				waveNormalizer *= (1 - Settings.agingTheta);
			}
			else if (Settings.agingFlag == Settings.noAge) {
				waveNormalizer /= (N+1)*1.0;
			}
			
			// Scale coefficients if Caudle aging is being used
			if (Settings.agingFlag == Settings.caudleAge) {
				for (int wavIndex = 0; 
						wavIndex < Transform.waveletCoefficients.get(j - j0).size();
						wavIndex++) {
					double newCoef = Settings.agingTheta*
							Transform.waveletCoefficients.get(j - j0).get(wavIndex);
					Transform.waveletCoefficients.get(j - j0).set(wavIndex,  newCoef);
				}
					
			}
			
			// Recursively compute coefficients if no aging is used
			else if (Settings.agingFlag == Settings.noAge){
				for (int wavIndex = 0; 
						wavIndex < Transform.waveletCoefficients.get(j - j0).size();
						wavIndex++) {
					double newCoef = (N)/(N+1.0)*Transform.waveletCoefficients.get(j - j0).get(wavIndex);
					Transform.waveletCoefficients.get(j - j0).set(wavIndex,  newCoef);
				}
			}
			
			// Subtract old samples effect if window aging is used
			else if (Settings.agingFlag == Settings.windowAge) {
				
				// Only remove a sample if there have been more than window size samples
				if (N > Settings.windowSize){ 
					double Xold = oldSamples[N % Settings.windowSize];
				
				
					//Loop through the translations for the wavelet basis functions
					int waveInd = 0;
					for (double k : Transform.waveletTranslates.get(j - j0)) {
					
						// Get the translated & scaled data point
						double xScaled = Math.pow(2, j) * Xold - k;
					
						// If the wavelet supports the data point, update the coefficient
						if (Wavelet.inSupport(xScaled)) {
							double scaleNew = Transform.waveletCoefficients.get(j - j0)
									.get(waveInd) - waveNormalizer*Wavelet.getPsiAt(xScaled);
							Transform.waveletCoefficients.get(j - j0).set(waveInd, scaleNew);
						
						}
						waveInd++;
					}
				}
			}
			
			//Loop through the translations for the wavelet basis functions
			int waveInd = 0;
			for (double k : Transform.waveletTranslates.get(j - j0)) {
				
				// Get the translated & scaled data point
				double xScaled = Math.pow(2, j) * Xnew - k;
				
				// If the wavelet supports the data point, update the coefficient
				if (Wavelet.inSupport(xScaled)) {
					double scaleNew = Transform.waveletCoefficients.get(j - j0)
							.get(waveInd) + waveNormalizer*Wavelet.getPsiAt(xScaled);
					Transform.waveletCoefficients.get(j - j0).set(waveInd, scaleNew);
					
				}
				waveInd++;
			}
		}
		
		
	} // end updateWaveletCoefficients
*/
	/**
	 * Find the maximum and minimum translation indices which
	 * support the incoming data point.
	 * 
	 * @param X : the data point
	 * @param j : the resolution level
	 * @return : An array containing the minimum and maximum
	 *           translation indices which support the data point
	 */
	private static double[] findRelevantKIndices(double X, int j) {
		//NOT YET VETTED, potentially could be used in future optimization
		
		// Get the max & min values for the wavelet's support
		double[] waveletMinMax = Wavelet.getSupport();
		
		double kMin = Math.ceil(Math.pow(2,j*X) - waveletMinMax[0]);
		double kMax = Math.floor(Math.pow(2,j*X) - waveletMinMax[1]);
		double[] kMinMax = new double[2];
		kMinMax[0] = kMin;
		kMinMax[1] = kMax;
		return (kMinMax);
		
	} //end findRelevantKIndices
	
	/**
	 * Initializes the translates for the scaling basis functions,
	 * based off of the maximum/minimum values supported and
	 * the starting resolution level.
	 * 
	 * Post: the translates in Transform are set appropriately
	 * @return the translates for the scaling basis functions
	 */
	public static void initializeTranslates() {
		
		// Initialize the scaling translates
		Transform.scalingTranslates = new ArrayList<Double> ();
		int startTranslate = (int) Math.floor((Math.pow(2,Settings.startLevel)*Settings.getMinimumRange())-Wavelet.getSupport()[1]);
		int stopTranslate = (int) Math.ceil((Math.pow(2,Settings.startLevel)*Settings.getMaximumRange())-Wavelet.getSupport()[0]);
		for (double k = startTranslate; k <= stopTranslate; k++){
			Transform.scalingTranslates.add(k);
		}
		
		
		// Initialize the wavelet translates if wavelets are being used
		if (Settings.waveletFlag) {
			Transform.waveletTranslates = new ArrayList<ArrayList<Double>> ();
		
			// Loop through resolutions
			for (int j = Settings.startLevel; j <= Settings.stopLevel; j++){
				
				ArrayList<Double> jTranslates = new ArrayList<Double> ();
				int startWTranslate = (int) Math.floor((Math.pow(2,j)*Settings.getMinimumRange())-Wavelet.getSupport()[1]);
				int stopWTranslate = (int) Math.ceil((Math.pow(2,j)*Settings.getMaximumRange())-Wavelet.getSupport()[0]);
				for (double k = startWTranslate; k <= stopWTranslate; k++){
					jTranslates.add(k);
				}
				Transform.waveletTranslates.add(jTranslates);
			}
		}
	} //end initializeTranslates
	
	/**
	 * Initializes the arrays to hold the scaling basis function
	 * and wavelet basis function coefficients based off
	 * of the resolution levels.
	 * 
	 * Post: the coefficients arrays are of the appropriate size.
	 */
	public static void initializeCoefficients() {
		N = 0;
		
		// Create window to store old samples
		if (Settings.agingFlag == Settings.windowAge) {
			oldSamples = new double[Settings.windowSize][2];
		}
		
		// Set all scaling coefficients to 0
		Transform.scalingCoefficients = new double[Transform.scalingTranslates.size()][Transform.scalingTranslates.size()];
		
/*		if (Settings.waveletFlag) {
			Transform.waveletCoefficients = new ArrayList<ArrayList<Double>> ();
			
			// Loop through resolutions
			for (int j = Settings.startLevel; j <= Settings.stopLevel; j++){
				
				// Set all wavelet at this resolution coefficients to 0
				ArrayList<Double> jCoefficients = new ArrayList<Double> (Collections.nCopies
													(Transform.waveletTranslates.get
															(j - Settings.startLevel).size(), 0.0));
				Transform.waveletCoefficients.add(jCoefficients);
			}
		}*/
	} //end initializeCoefficients
	
	/**
	 * Calculates the density at each discrete point in the
	 * range pre-specified.
	 * 
	 * Pre: the coefficient matrices have been created and
	 *      are the proper size
	 * @return the normalized density estimate
	 */
	private static double[][] getDensity() {
		
		ArrayList<Double> density = new ArrayList<Double> ();
		double scaleNormalizer = Math.pow(2, Settings.startLevel/2.0);
		
		// Calculate un-normalized density for each point in domain
		for (double i = Settings.getMinimumRange(); 
				i < Settings.getMaximumRange(); i += Settings.discretization) {
			
			// Density at point i
			double iDense = 0.0;
			
			// Cycle through translates for each point
			int scalIndex = 0;
			for (double k : Transform.scalingTranslates) {
				double Xi = Math.pow(2, Settings.startLevel)*i - k;
				
				// Only update if the point is supported
				if (Wavelet.inSupport(Xi)) {
					iDense += Transform.scalingCoefficients.get(scalIndex) 
							  * Wavelet.getPhiAt(Xi) * scaleNormalizer;
				}
				scalIndex++;
			}
			density.add(iDense);
		}
		
		if (Settings.waveletFlag) {
			density = addWaveDensity(density);
		}
		
		//Normalize density
		density = normalizeDensity(density);
		
		return(density);
	}
	
	/**
	 * Adds the contribution from the wavelet basis functions
	 * @param density The density from the scaling basis functions
	 * @return The complete unnormalized density
	 */
	private static double[][] addWaveDensity(ArrayList<Double> density) {
		
		// Short hand for start leve
		int j0 = Settings.startLevel;
			
		// Loop through for each resolution level
		for (int j = Settings.startLevel; j <= Settings.stopLevel; j++) {
			
			double waveNormalizer = Math.pow(2, j/2.0);
			
			// Calculate un-normalized density for each point in domain
			int domainIndex = 0;
			for (double i = Settings.getMinimumRange(); 
					i < Settings.getMaximumRange(); i += Settings.discretization) {
				
				// Density at point i
				double iDense = 0.0;
				
				// Cycle through translates for each point
				int wavIndex = 0;
				for (double k : Transform.waveletTranslates.get(j - j0)) {
					double Xi = Math.pow(2, j)*i - k;
					
					// Only update if the point is supported
					if (Wavelet.inSupport(Xi)) {
						iDense += Transform.waveletCoefficients.get(j - j0).get(wavIndex) 
								  * Wavelet.getPsiAt(Xi) * waveNormalizer;
					}
					wavIndex++;
				}
				density.set(domainIndex, iDense + density.get(domainIndex));
				domainIndex++;
			}
		}
		
		
		return density;
	}

	/**
	 * Takes in an un-normalized density estimate and returns the
	 * normalized version, using the normalization procedure from
	 * Gajek (1986) 'On improving density estimators which are
	 * not Bona Fide functions'.
	 * 
	 * @param unNormDensity : the un-normalized density estimate
	 *                        over the domain range.
	 * @return the normalized density.
	 */
	private static ArrayList<Double> normalizeDensity(ArrayList<Double> unNormDensity){
		
		ArrayList<Double> normDens = unNormDensity;
		int iter = 0;
		double threshold = Math.pow(10, -8);
		double densityDomainSize = Settings.getMaximumRange() - Settings.getMinimumRange();
		
		while (iter < 1000) {
			
			// Zero negative points
			for (int i = 0; i < unNormDensity.size(); i++) {
				if (normDens.get(i) < 0.0) {
					normDens.set(i,0.0);
				}
			}
			
			// Sum over probability density over interval
			double integralSum = 0.0;
			for (int i = 0; i < unNormDensity.size(); i++) {
				integralSum += normDens.get(i)*Settings.discretization;
			}
			
			// Return if error is under threshold
			if (Math.abs(integralSum - 1) < threshold) {
				return (normDens);
			}
			
			// Modify density so that it integrates to 1
			double normalizeConstant = (integralSum - 1) / densityDomainSize;
			for (int i = 0; i < unNormDensity.size(); i++) {
				
				normDens.set(i, normDens.get(i) - normalizeConstant);
			}
			
			iter++;
		}
		
		// Settle for the current approximation
		return normDens;
	} //end normalizeDensity
	
	/**
	 * Updates the plotted density based on current coefficients
	 */
	public static void updateDensity() {
		
		double[][] normDensity = getDensity();
		
		// UPDATE PLOT DENSITY

	}
	
}

