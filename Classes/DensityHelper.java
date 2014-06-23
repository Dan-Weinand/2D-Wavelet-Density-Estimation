/**
 * Helper functions for the density estimation algorithm.
 * 
 * @author Daniel Weinand & Gedeon Nyengele
 * 
 */
import java.util.ArrayList;
import java.util.Collections;


public class DensityHelper {
	
	private static double[][] oldSamples;	// The old samples in the window
	private static int N;					// How many samples have been read in
	private static double postProb = 0.0;   // Posterior probability score
	
	// Whether or not to calculate and print a cross validation score
	private static boolean postProbOn = false; 
	
	// The phi values over the density domain gridlines
	private static ArrayList<ArrayList<Double>> scalePhisHere;
	private static ArrayList<ArrayList<ArrayList<Double>>> wavePhisHere;
	private static ArrayList<ArrayList<ArrayList<Double>>> wavePsisHere;
	
	
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
		
		// Calculate Posterior probability of data point
		if (postProbOn && inRange(Xnew[0]) && inRange(Xnew[1])) {
			double[][] density = calculateDensity();
			int X1newInd = (int) Math.floor((Xnew[0] - Settings.getMinimumRange())
					/ Settings.discretization);
			int X2newInd = (int) Math.floor((Xnew[1] - Settings.getMinimumRange())
					/ Settings.discretization);
			double XnewDensity = density[X1newInd][X2newInd];
			if (XnewDensity <= .0001) {
				postProb -= 9.2;
			}
			else {
				postProb += Math.log(XnewDensity);
			}
			System.out.println(postProb);
		}
		else {
			postProb -= 9.2;
		}
		// End calculating posterior probability
		
		if (Settings.waveletFlag) {
			updateWaveletCoefficients(Xnew);
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
		
		// Scale coefficients if Caudle aging is being used
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
				
		}
		
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
								double phi2Here = scaleNormalizer*Wavelet.getPhiAt(x2Scaled);
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
						double phi2Here = scaleNormalizer*Wavelet.getPhiAt(x2Scaled);
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
	private static void updateWaveletCoefficients(double[] Xnew) {
		// Loop through resolutions
		for (double j = Settings.startLevel; j <= Settings.stopLevel; j++){
			
			int waveInd = (int) Math.floor(j - Settings.startLevel);
			
			// The normalizing constant for this resolution of wavelet functions
			double waveNormalizer = Math.pow(2, j/2.0);
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
			
			// Scale coefficients if Caudle aging is being used
			if (Settings.agingFlag == Settings.caudleAge) {
				for (int x1Index = 0; 
						x1Index < Transform.psiPsiCoefficients.get(waveInd).length;
						x1Index++) {
					for (int x2Index = 0; 
							x2Index < Transform.psiPsiCoefficients.get(waveInd).length;
							x2Index++) {
						Transform.psiPsiCoefficients.get(waveInd)[x1Index][x2Index] *= Settings.agingTheta;
						Transform.psiPhiCoefficients.get(waveInd)[x1Index][x2Index] *= Settings.agingTheta;
						Transform.phiPsiCoefficients.get(waveInd)[x1Index][x2Index] *= Settings.agingTheta;
					}
				}
					
			}
			
			// Recursively compute coefficients if no aging is used
			else if (Settings.agingFlag == Settings.noAge){
				for (int x1Index = 0; 
						x1Index < Transform.psiPsiCoefficients.get(waveInd).length;
						x1Index++) {
					for (int x2Index = 0; 
							x2Index < Transform.psiPsiCoefficients.get(waveInd).length;
							x2Index++) {
						Transform.psiPsiCoefficients.get(waveInd)[x1Index][x2Index] *= N/(N+1.0);
						Transform.psiPhiCoefficients.get(waveInd)[x1Index][x2Index] *= N/(N+1.0);
						Transform.phiPsiCoefficients.get(waveInd)[x1Index][x2Index] *= N/(N+1.0);
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
				
				
					//Loop through the x1 translations for the wavelet basis functions
					int x1Ind = 0;
					for (double k1 : Transform.waveletTranslates.get(waveInd)) {
					
						// Get the translated & scaled data point
						double x1Scaled = Math.pow(2, j) * x1 - k1;
					
						// If the x1 wavelet supports the data point, check through x2 supports
						if (Wavelet.inSupport(x1Scaled)) {
							double phi1Here = waveNormalizer*Wavelet.getPhiAt(x1Scaled);
							double psi1Here = waveNormalizer*Wavelet.getPsiAt(x1Scaled);
							int x2Ind = 0;
							
							//Loop through the x2 translations for the wavelet basis functions
							for (double k2 : Transform.waveletTranslates.get(waveInd)) {
								
								// Get the translated & scaled data point
								double x2Scaled = Math.pow(2, j) * x2 - k2;
								
								// Update appropriate coefficient iff in x2 support
								if (Wavelet.inSupport(x2Scaled)) {
									double phi2Here = waveNormalizer*Wavelet.getPhiAt(x2Scaled);
									double psi2Here = waveNormalizer*Wavelet.getPsiAt(x2Scaled);
									double psiPsiSub = psi1Here*psi2Here/Settings.windowSize;
									double psiPhiSub = psi1Here*phi2Here/Settings.windowSize;
									double phiPsiSub = phi1Here*psi2Here/Settings.windowSize;
									
									Transform.psiPsiCoefficients.get(waveInd)[x1Ind][x2Ind] -= psiPsiSub;
									Transform.psiPhiCoefficients.get(waveInd)[x1Ind][x2Ind] -= psiPhiSub;
									Transform.phiPsiCoefficients.get(waveInd)[x1Ind][x2Ind] -= phiPsiSub;
								} // end updating
								x2Ind++;
							} // end x2 looping										
						} // end in x1 support
						x1Ind++;
					} // end x1 looping
				}
			}
			
			// The x1 and x2 coordinates of the incoming point
			double x1 = Xnew[0];
			double x2 = Xnew[1];
		
		
			//Loop through the x1 translations for the scaling basis functions
			int x1Ind = 0;
			for (double k1 : Transform.waveletTranslates.get(waveInd)) {
			
				// Get the translated & scaled data point
				double x1Scaled = Math.pow(2, j) * x1 - k1;
			
				// If the x1 wavelet supports the data point, check through x2 supports
				if (Wavelet.inSupport(x1Scaled)) {
					double phi1Here = waveNormalizer*Wavelet.getPhiAt(x1Scaled);
					double psi1Here = waveNormalizer*Wavelet.getPsiAt(x1Scaled);
					int x2Ind = 0;
					
					//Loop through the x2 translations for the wavelet basis functions
					for (double k2 : Transform.waveletTranslates.get(waveInd)) {
						
						// Get the translated & scaled data point
						double x2Scaled = Math.pow(2, j) * x2 - k2;
						
						// Update appropriate coefficient iff in x2 support
						if (Wavelet.inSupport(x2Scaled)) {
							double phi2Here = waveNormalizer*Wavelet.getPhiAt(x2Scaled);
							double psi2Here = waveNormalizer*Wavelet.getPsiAt(x2Scaled);
							
							double psiPsiAdd = psi1Here*psi2Here*ageNorm;
							double psiPhiAdd = psi1Here*phi2Here*ageNorm;
							double phiPsiAdd = phi1Here*psi2Here*ageNorm;
							
							Transform.psiPsiCoefficients.get(waveInd)[x1Ind][x2Ind] += psiPsiAdd;
							Transform.psiPhiCoefficients.get(waveInd)[x1Ind][x2Ind] += psiPhiAdd;
							Transform.phiPsiCoefficients.get(waveInd)[x1Ind][x2Ind] += phiPsiAdd;
						} // end updating
						x2Ind++;
					} // end x2 looping										
				} // end in x1 support
				x1Ind++;
			} // end x1 looping
			
		} // End resolution looping
	} // End update wavelet coefficients


	/**
	 * Find the maximum and minimum translation indices which
	 * support the incoming data point.
	 * 
	 * @param X : the data point
	 * @param j : the resolution level
	 * @return : An array containing the minimum and maximum
	 *           translation indices which support the data point
	 */
	private static int[] findRelevantKIndices(double X, double j) {
		
		// Get the max & min values for the wavelet's support
		double[] waveletMinMax = Wavelet.getSupport();
		
		int jInd = (int) Math.floor(j - Settings.startLevel);
		int minimumAllowedTranslate;
		if (Settings.waveletFlag) {
			minimumAllowedTranslate = (int) Math.floor(Transform.waveletTranslates.get(jInd).get(0));
		}
		else {
			minimumAllowedTranslate = (int) Math.floor(Transform.scalingTranslates.get(0));
		}
		
		
		int kMax = (int) Math.ceil(Math.pow(2,j) * X - waveletMinMax[0]);
		int kMin = (int) Math.floor(Math.pow(2,j) * X - waveletMinMax[1]);
		int[] kMinMax = new int[2];
		kMinMax[0] = kMin - minimumAllowedTranslate;
		kMinMax[1] = kMax - minimumAllowedTranslate;
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
			for (double j = Settings.startLevel; j <= Settings.stopLevel; j++){
				
				ArrayList <Double> jTranslates = new ArrayList<Double> ();
				startTranslate = (int) Math.floor((Math.pow(2, j)*Settings.getMinimumRange())-Wavelet.getSupport()[1]);
				stopTranslate = (int) Math.ceil((Math.pow(2, j)*Settings.getMaximumRange())-Wavelet.getSupport()[0]);
				for (double k = startTranslate; k <= stopTranslate; k++){
					jTranslates.add(k);
				}
				Transform.waveletTranslates.add(jTranslates);
			}
		}
		
		initializePhiGrid();
	} //end initializeTranslates
	
	/**
	 * Initializes the arrays which hold the phi and psi values over the
	 * density domain gridlines
	 */
	private static void initializePhiGrid() {
		
		scalePhisHere = new ArrayList<ArrayList<Double>> ();
		int numGridLines = getNumGridlines();
		double scaleNormalizer = Math.pow(2, Settings.startLevel/2.0);
		double i1 = Settings.getMinimumRange();
		
		for (int x1Ind = 0; x1Ind < numGridLines; x1Ind++)
			{		
			
			ArrayList<Double> thesePhis = new ArrayList<Double> ();
			i1 += Settings.discretization;
			int[] i1RelevantIndices = findRelevantKIndices(i1, Settings.startLevel);
			int k1Max = i1RelevantIndices[1];
			int k1Min = i1RelevantIndices[0];
		
		
			// Cycle through relevant X1 translates for the line
			for (int k1Ind = k1Min; k1Ind < k1Max; k1Ind++) {
			
				double k1 = Transform.scalingTranslates.get(k1Ind);
				double Xi1 = Math.pow(2, Settings.startLevel) * i1 - k1;
				double phi1Here = Wavelet.getPhiAt(Xi1) * scaleNormalizer;
				thesePhis.add(phi1Here);
			
			}
			scalePhisHere.add(thesePhis);
		}
		
		if (Settings.waveletFlag) {
			
			wavePhisHere = new ArrayList<ArrayList<ArrayList<Double>>> ();
			wavePsisHere = new ArrayList<ArrayList<ArrayList<Double>>> ();
			
			// Loop through resolutions
			for (double j = Settings.startLevel; j <= Settings.stopLevel; j++){
				int waveInd = (int) Math.floor(j - Settings.startLevel);
				double waveNormalizer = Math.pow(2, j/2.0);
				
				ArrayList<ArrayList<Double>> jPhis = new ArrayList<ArrayList<Double>> ();
				ArrayList<ArrayList<Double>> jPsis = new ArrayList<ArrayList<Double>> ();
				
				// Loop through grid
				i1 = Settings.getMinimumRange();
				for (int x1Ind = 0; x1Ind < numGridLines; x1Ind++)
				{		
					ArrayList<Double> thesePhis = new ArrayList<Double> ();
					ArrayList<Double> thesePsis = new ArrayList<Double> ();
					i1 += Settings.discretization;
					int[] i1RelevantIndices = findRelevantKIndices(i1, j);
					int k1Max = i1RelevantIndices[1];
					int k1Min = i1RelevantIndices[0];
			
			
					// Cycle through relevant X1 translates for the line
					for (int k1Ind = k1Min; k1Ind < k1Max; k1Ind++) {
				
						double k1 = Transform.waveletTranslates.get(waveInd).get(k1Ind);
						double Xi1 = Math.pow(2, j) * i1 - k1;
						
						double phi1Here = Wavelet.getPhiAt(Xi1) * waveNormalizer;
						thesePhis.add(phi1Here);
						
						double psi1Here = Wavelet.getPsiAt(Xi1) * waveNormalizer;
						thesePsis.add(psi1Here);
				
					} // End cycling through relevant x1 translates
					jPhis.add(thesePhis);
					jPsis.add(thesePsis);
				} // End cycling through x1 grid
				wavePhisHere.add(jPhis);
				wavePsisHere.add(jPsis);
				
			} // End looping across resolutions
		} // End defining wavelet phis and psis

	}
	
	/**
	 * Initializes the arrays to hold the scaling basis function
	 * and wavelet basis function coefficients based off
	 * of the resolution levels.
	 * 
	 * Post: the coefficients arrays are of the appropriate size.
	 */
	public static void initializeCoefficients() {
		N = 0;
		postProb = 0.0;
		
		// Create window to store old samples
		if (Settings.agingFlag == Settings.windowAge) {
			oldSamples = new double[Settings.windowSize][2];
		}
		
		// Set all scaling coefficients to 0
		Transform.scalingCoefficients = new double[Transform.scalingTranslates.size()][Transform.scalingTranslates.size()];
		
		if (Settings.waveletFlag) {
			Transform.psiPsiCoefficients = new ArrayList<double[][]> ();
			Transform.psiPhiCoefficients = new ArrayList<double[][]> ();
			Transform.phiPsiCoefficients = new ArrayList<double[][]> ();
			
			// Loop through resolutions
			for (double j = Settings.startLevel; j <= Settings.stopLevel; j++){
				
				int waveInd = (int) Math.floor(j - Settings.startLevel);
				int jSize = Transform.waveletTranslates.get(waveInd).size();
				
				// Set all wavelet at this resolution coefficients to 0
				double[][] thisPsiPsi = new double[jSize][jSize];
				double[][] thisPsiPhi = new double[jSize][jSize];
				double[][] thisPhiPsi = new double[jSize][jSize];
				Transform.psiPsiCoefficients.add(thisPsiPsi);
				Transform.psiPhiCoefficients.add(thisPsiPhi);
				Transform.phiPsiCoefficients.add(thisPhiPsi);
			}
		}
	} //end initializeCoefficients
	
	/**
	 * Updates the density plot
	 */
	public static void newDensity() {
		
		double[][] density = calculateDensity();
		int numGridLines = getNumGridlines();
		
		// Find the maximum normalized density
		float maxDens = (float) 0.0;
		for (int x1Ind = 0; x1Ind < numGridLines; x1Ind++) {
			for (int x2Ind = 0; x2Ind < numGridLines; x2Ind++) {
				if (density[x1Ind][x2Ind] > maxDens) maxDens = (float) density[x1Ind][x2Ind];
			}
		}
		
		// Update the plot's density
		DensityModel.updateDensity(density, maxDens);
	}
	
	/**
	 * Calculates the density at each discrete point in the
	 * range pre-specified.
	 * 
	 * Pre: the coefficient matrices have been created and
	 *      are the proper size
	 * @return The normalized density estimate
	 */
	public static double[][] calculateDensity() {
		
		int numGridLines = getNumGridlines();
		double[][] density = new double[numGridLines][numGridLines];
		
		double i1 = Settings.getMinimumRange();
		// Calculate un-normalized density for each point in domain, looping across X1
		for (int x1Ind = 0; x1Ind < numGridLines; x1Ind++)
			{
			i1 += Settings.discretization;
			int[] i1RelevantIndices = findRelevantKIndices(i1, Settings.startLevel);
			int k1Max = i1RelevantIndices[1];
			int k1Min = i1RelevantIndices[0];					
			
			// Cycle through relevant X1 translates for the line
			for (int k1Ind = k1Min; k1Ind < k1Max; k1Ind++) {
				
				double phi1Here = scalePhisHere.get(x1Ind).get(k1Ind - k1Min);
				
				// Loop across X2 dimension
				double i2 = Settings.getMinimumRange();
				for (int x2Ind = 0; x2Ind < numGridLines; x2Ind++)
				{
					i2 += Settings.discretization;
					int[] i2RelevantIndices = findRelevantKIndices(i2, Settings.startLevel);
					int k2Max = i2RelevantIndices[1];
					int k2Min = i2RelevantIndices[0];
										
					// Cycle through relevant X2 translates
					for (int k2Ind = k2Min; k2Ind < k2Max; k2Ind++) {
						
						double phi2Here = scalePhisHere.get(x2Ind).get(k2Ind - k2Min);
						
						density[x1Ind][x2Ind] += Transform.scalingCoefficients[k1Ind][k2Ind] *
								phi1Here*phi2Here;
					} // End cycling relevant X2 translates
				} // End cycling across X2 gridlines
			} // End cycling relevant X1 translates
		} // End cycling across X1 gridlines
		
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
	private static double[][] addWaveDensity(double[][] density) {
		
		int numGridLines = getNumGridlines();
			
		// Loop through for each resolution level
		for (double j = Settings.startLevel; j <= Settings.stopLevel; j++) {
			
			int waveInd = (int) Math.floor(j - Settings.startLevel);
			ArrayList<ArrayList<Double>> jPhis = wavePhisHere.get(waveInd);
			ArrayList<ArrayList<Double>> jPsis = wavePsisHere.get(waveInd);
			
			double i1 = Settings.getMinimumRange();
			// Calculate un-normalized density for each point in domain, looping across X1
			for (int x1Ind = 0; x1Ind < numGridLines; x1Ind++)
				{
				i1 += Settings.discretization;
				int[] i1RelevantIndices = findRelevantKIndices(i1, j);
				int k1Max = i1RelevantIndices[1];
				int k1Min = i1RelevantIndices[0];					
				
				// Cycle through relevant X1 translates for the line
				for (int k1Ind = k1Min; k1Ind < k1Max; k1Ind++) {
					
					double phi1Here = jPhis.get(x1Ind).get(k1Ind - k1Min);
					double psi1Here = jPsis.get(x1Ind).get(k1Ind - k1Min);
					
					// Loop across X2 dimension
					double i2 = Settings.getMinimumRange();
					for (int x2Ind = 0; x2Ind < numGridLines; x2Ind++)
					{
						i2 += Settings.discretization;
						int[] i2RelevantIndices = findRelevantKIndices(i2, j);
						int k2Max = i2RelevantIndices[1];
						int k2Min = i2RelevantIndices[0];
											
						// Cycle through relevant X2 translates
						for (int k2Ind = k2Min; k2Ind < k2Max; k2Ind++) {
							
							double phi2Here = jPhis.get(x2Ind).get(k2Ind - k2Min);
							double psi2Here = jPsis.get(x2Ind).get(k2Ind - k2Min);
							
							density[x1Ind][x2Ind] += Transform.psiPsiCoefficients.get(waveInd)[k1Ind][k2Ind] *
									psi1Here*psi2Here;
							density[x1Ind][x2Ind] += Transform.psiPhiCoefficients.get(waveInd)[k1Ind][k2Ind] *
									psi1Here*phi2Here;
							density[x1Ind][x2Ind] += Transform.phiPsiCoefficients.get(waveInd)[k1Ind][k2Ind] *
									phi1Here*psi2Here;
						} // End cycling relevant X2 translates
					} // End cycling across X2 gridlines
				} // End cycling relevant X1 translates
			} // End cycling across X1 gridlines
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
	private static double[][] normalizeDensity(double[][] unNormDensity){
		
		// The number of gridlines across the density domain
		int numGridLines = getNumGridlines();
		
		int iter = 0;                        // How many normalization cycles have been performed
		double threshold = Math.pow(10, -8); // The maximum acceptable error
		double densityDomainSize = Settings.getMaximumRange() - Settings.getMinimumRange();
		double[][] normDens = unNormDensity;
		
		double sumall = 0.0;
		for (int i1 = 0; i1 < numGridLines; i1++) {
			
			for (int i2 = 0; i2 < numGridLines; i2++) {
				
				sumall += unNormDensity[i1][i2];
			}
		}
		
		while (iter < 1000) {
			
			// Zero negative points
			for (int i1 = 0; i1 < numGridLines; i1++) {
				
				for (int i2 = 0; i2 < numGridLines; i2++) {
					
					if (normDens[i1][i2] < 0.0) {
						normDens[i1][i2] = 0.0;
					}
				}
			}
			
			// Sum over probability density over interval
			double integralSum = 0.0;
			for (int i1 = 0; i1 < numGridLines; i1++) {
				for (int i2 = 0; i2 < numGridLines; i2++) {
					integralSum += normDens[i1][i2] * Math.pow(Settings.discretization, 2);
				}
			}
			
			// Return if error is under threshold
			if (Math.abs(integralSum - 1) < threshold) {
				return (normDens);
			}
			
			// Modify density so that it integrates to 1
			double normalizeConstant = (integralSum - 1) / Math.pow(densityDomainSize, 2);
			for (int i1 = 0; i1 < numGridLines; i1++) {
				for (int i2 = 0; i2 < numGridLines; i2++) {
					normDens[i1][i2] -= normalizeConstant;
				}
			}
			
			iter++;
		}
		
		// Settle for the current approximation
		return normDens;
	} //end normalizeDensity
	
	
	/**
	 * Gives the number of gridlines in the density domain
	 * @return the number of x1/x2 gridlines in the density domain
	 */
	private static int getNumGridlines() {
		return (int) Math.ceil((Settings.getMaximumRange() - Settings.getMinimumRange())
				/Settings.discretization);
	}
	
}

