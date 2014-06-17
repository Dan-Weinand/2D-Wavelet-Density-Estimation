/**
 * Thread to perform all calculations and update plots.
 * @author Gedeon Nyengele & Daniel Weinand.
 */

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.List;

import javax.swing.JButton;
import javax.swing.JTextField;
import javax.swing.SwingWorker;

import org.sf.surfaceplot.SurfaceCanvas;



public class DensityRunner extends SwingWorker<Object, Integer>{
	private JTextField sampLabel;
	private boolean paused;
	private boolean terminated;
	private JButton startButton, stopButton, settingsButton;
	private SurfaceCanvas canvas;
	private DensityModel dataModel;
	
	/**
	 * Constructor.
	 * @param smpLabel         	: JLabel component that shows the current sample number.
	 * @param startButton		: reference to the the applet's start button.
	 * @param stopButton		: reference to the the applet's stop button.
	 * @param settingsButton	: reference to the the applet's settings button.
	 * @param plotCanvas		: reference to the the plot.
	 * @param dtModel			: reference to the density data model.
	 */
	public DensityRunner ( JTextField smpLabel, 
			               JButton startButton, 
			               JButton stopButton, 
			               JButton settingsButton, 
			               SurfaceCanvas plotCanvas, 
			               DensityModel dtModel  )
	{
		sampLabel            = smpLabel;
		paused               = false;
		terminated           = false;
		this.startButton     = startButton;
		this.stopButton      = stopButton;
		this.settingsButton  = settingsButton;
		
		canvas    = plotCanvas;
		dataModel = dtModel;
		
		
	}
	
	/**
	 * Performs calculations necessary for the iterative plotting.
	 * This method is an implementation of the SwingWorker's abstract 
	 * method protected <Type> doInBackground.
	 */
	protected Object doInBackground(){		
		
		// Initialize the algorithm variables
		try { Wavelet.init( Settings.waveletType ); } 
		catch ( IOException e ) {}
    	DensityHelper.initializeTranslates();
    	DensityHelper.initializeCoefficients();
    	
    	// Intialize the current sample index to 1.
    	int sampInd = 1;
    	
    	// Buffer reader used to read the user-provided sample data file.
		BufferedReader dataReader;
		
		try
		{
			// Read the user file.
			dataReader = new BufferedReader( new FileReader( Settings.dataFile ) ); 
			
			// Perform next calculations only if there is a data sample left and
			// the density runner is not killed.
			while( dataReader.ready() && !isCancelled() )
			{
				// Wait if user so requested.
				while( paused )
				{
					// While pausing, close the buffer reader and stop if the user 
					// terminates the execution of the Density Runner.
					if( terminated ) 
					{
						dataReader.close();  
						return null;
					} // end if( terminated ).
					
					try
					{
						// Sleep while the Runner is still on pause.
						Thread.sleep( 500 );
						
					} // end try Thread.sleep( 500 ).
					catch(InterruptedException ex){ex.printStackTrace();}
				} // end while( paused ).
				
				// Get the next sample.
				double[] sample = new double[2];
				String lineContent = dataReader.readLine();
				String[] parts = lineContent.split(",");
				sample[0] = Double.parseDouble(parts[0]);
				sample[1] = Double.parseDouble(parts[1]);

				
				// Update the wavelet's coefficients using the new sample.
				DensityHelper.updateCoefficients( sample );
				
				// Display the current density at the user-specified frequency.
				if( sampInd % Settings.updateFrequency == 0 )
				{
					// Update the density using the current coefficients.
					DensityHelper.newDensity();
// ----------------------------------------------------------------------					

					//
					
					// Send the current sample index to the process method
					// for updating the sample index label in the applet.
					// !! the pubish method is defined in the SwingWorker class !!.
					publish( sampInd );
					
				} // end if( sampInd % Settings.updateFrequency == 0 ).

				try
				{
					// Give some time to the event queue
					// to handle the plot update.
					Thread.sleep( 1 );
					
				} // end try{ Thread.sleep ( 1 ) }
				catch( InterruptedException ex){ex.printStackTrace();}
				
				// Increment the sample index.
				sampInd++;
			} // end while( dataReader.ready() && !isCancelled() )
			
			// Close the buffer reader when execution terminates.
			dataReader.close();
			
		} // end the main try{}
		catch(Exception ex){ ex.printStackTrace();}
		
		// doInBackground has to return an object.
		// We return null as there is no meaningful object returned by this method.
		return null;
	} // end method public Object doInBackground().
	
	/**
	 * Resets the applet's UI components to their original states.
	 * This method overrides the SwingWorker class's protected void done().
	 */
	protected void done(){
		try
		{
			startButton.setEnabled( true );    // Enable the start button in the applet.
			settingsButton.setEnabled( true ); // Enable the settings button in the applet.
			stopButton.setEnabled( false );    // Disable the stop button in the applet.
			if( stopButton.getText() == "Resume" ) stopButton.setText( "Stop" ); // Assign a proper label to the stop button.
		} // end try{}
		catch(Exception ex){}
	} // end method protected void done().
	
	/**
	 * Updates the sample index label in the applet and repaints the plot in the applet.
	 * This method overrides the protected void process() method in SwingWorker class.
	 * @param published : List of published sample indexes from doInBackground method.
	 */
	protected void process( List<Integer> published ){
		
		// Get the current sample index from the List.
		Integer sampIndex = published.get( published.size() - 1 );
		
		// Update the sample index label in the applet.
		sampLabel.setText( sampIndex.toString() );
		canvas.setModel(dataModel);
		
		// Repaint the plot.
		canvas.repaint();
	} // end method protected void process().
	
	/**
	 * Pauses the Density Runner.
	 */
	public void pause()
	{
		paused = true;
	} // end method pause.
	
	/**
	 * Resumes the Density Runner.
	 */
	public synchronized void resume(){
		paused = false;
		notify();
	} // end method resume.
	
	/**
	 * Terminates the Density Runner.
	 */
	public void terminate()
	{
		terminated = true;
	}// end method terminate().
	
} // end class DensityRunner.
