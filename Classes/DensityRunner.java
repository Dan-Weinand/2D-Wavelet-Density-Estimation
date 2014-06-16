/**
 * Thread to perform all calculations and update plots.
 * @author Gedeon Nyengele & Daniel Weinand.
 */

import java.awt.Color;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.List;

import javax.swing.JButton;
import javax.swing.JTextField;
import javax.swing.SwingWorker;

import de.erichseifert.gral.data.DataTable;
import de.erichseifert.gral.plots.XYPlot;
import de.erichseifert.gral.plots.areas.AreaRenderer;
import de.erichseifert.gral.plots.areas.DefaultAreaRenderer2D;
import de.erichseifert.gral.plots.axes.Axis;
import de.erichseifert.gral.plots.lines.DefaultLineRenderer2D;
import de.erichseifert.gral.plots.lines.LineRenderer;
import de.erichseifert.gral.ui.InteractivePanel;
import de.erichseifert.gral.util.Insets2D;


public class DensityRunner extends SwingWorker<Object, Integer>{
	private JTextField sampLabel;
	private boolean paused;
	private boolean terminated;
	private JButton startButton, stopButton, settingsButton;
	private  final LineRenderer lines = new DefaultLineRenderer2D();
	private  final AreaRenderer area  = new DefaultAreaRenderer2D();
	private  final Color invis = new Color(0.0f, 0.3f, 1.0f, 0);
	private  final Color lineColor = new Color(0.0f, 0.3f, 1.0f);
	private  final Color areaColor = new Color(0.0f, 0.3f, 1.0f, 0.3f);
	private  DataTable densityTable = null;
	private  XYPlot dataPlot = null;
	private  InteractivePanel dataPanel = null;
	private  double maxHeight = 0.0;
	private  int MAX_WIDTH = 0;
	private  int MAX_HEIGHT = 0;
	
	/**
	 * Constructor.
	 * @param smpLabel         	: JLabel component that shows the current sample number.
	 * @param plotWidth			: width of the plotting area.		
	 * @param plotHeight		: height of the plotting area.
	 * @param startButton		: reference to the the applet's start button.
	 * @param stopButton		: reference to the the applet's stop button.
	 * @param settingsButton	: reference to the the applet's settings button.
	 * @param dtPlot			: reference to the the plot.
	 * @param dtPanel			: reference to the the plot's container.
	 */
	public DensityRunner ( JTextField smpLabel, int plotWidth, int plotHeight, JButton startButton, JButton stopButton, JButton settingsButton, XYPlot dtPlot, InteractivePanel dtPanel  )
	{
		sampLabel            = smpLabel;
		paused               = false;
		terminated           = false;
		this.startButton     = startButton;
		this.stopButton      = stopButton;
		this.settingsButton  = settingsButton;
		dataPlot             = dtPlot;
		dataPanel			 = dtPanel;
		MAX_WIDTH            = plotWidth;
		MAX_HEIGHT           = plotHeight;
		
		// Initialize the density table to zeros.
		densityTable = new DataTable( 2, Double.class );
		for( double x = Settings.getMinimumRange(); x <= Settings.getMaximumRange(); x += Settings.discretization )
		{
			double y = 0.0;
			densityTable.add( x, y );
		}
		
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
					catch(InterruptedException ex){}
				} // end while( paused ).
				
				// Get the next sample.
				double xNew = Double.parseDouble( dataReader.readLine() );
				
				// Update the wavelet's coefficients using the new sample.
				DensityHelper.updateCoefficients( xNew );
				
				// Display the current density at the user-specified frequency.
				if( sampInd % Settings.updateFrequency == 0 )
				{
					// Update the density using the current coefficients.
					DensityHelper.updateDensity( densityTable );
					
					// Perform routines to prepare for plotting.
					dataPlot.clear();  // Clear the content of the current plot.
					dataPlot.add( densityTable );  // Add the new density data to the plot.
					dataPlot.setLineRenderer( densityTable, lines );  // Style the plot's lines.
			        dataPlot.setAreaRenderer( densityTable, area );  // Style the area below the curve.
			        
			        // Create some padding between the plot and the window.
			        // This is necessary as otherwise the tick marks on the x-axis will not show.
			        dataPlot.setInsets( new Insets2D.Double( 5, 50, 40, 40 ) );  
			        dataPlot.getAxisRenderer( XYPlot.AXIS_Y ).setIntersection( -Double.MAX_VALUE );  // Push the y-axis all the way to the left of the plot.
			        
			        // Apply colors to the plot's lines, points, and area below the curve.
			        dataPlot.getPointRenderer( densityTable ).setColor( invis );
			        dataPlot.getLineRenderer( densityTable ).setColor( lineColor );
			        dataPlot.getAreaRenderer( densityTable ).setColor( areaColor );
			        
					// Fix the maximum height on the y-axis.
					Axis yAx = dataPlot.getAxis( XYPlot.AXIS_Y );  // get a handle of the y-axis.
					double curYMax = ( Double ) yAx.getMax(); // get the maximum y value.
					
					// Adjust the maximum y value on the plot.
					if ( curYMax > maxHeight ) {
						maxHeight = curYMax + 0.8; // Adds some space between the the curve and the plot's frame.
					}
					
					// Reset the y-axis new max value.
					Axis nYAx = new Axis(); // Create a brand new axis.
					nYAx.setMax( maxHeight ); // Set the maximum value on this new axis ( = maxHeight ).
					nYAx.setMin( yAx.getMin() ); // Set the minimum value on this new axis equal the minimum value on the old axis.
					dataPlot.setAxis( "y",  nYAx ); // Replace the old y axis with the new axis in the plot.
					
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
				catch( InterruptedException ex){}
				
				// Increment the sample index.
				sampInd++;
			} // end while( dataReader.ready() && !isCancelled() )
			
			// Close the buffer reader when execution terminates.
			dataReader.close();
			
		} // end the main try{}
		catch(Exception ex){}
		
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
		
		
		// Repaint the plot.
		dataPanel.paintImmediately( 0, 0, MAX_WIDTH, MAX_HEIGHT );
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
