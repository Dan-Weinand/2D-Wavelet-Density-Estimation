/**
 * Class to both display the current density estimate,
 * and allow the user to select various parameters.
 * 
 * @author Daniel Weinand & Gedeon Nyengele
 * 
 */

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JApplet;
import javax.swing.JButton;
import javax.swing.JPanel;
import javax.swing.JTextField;

import de.erichseifert.gral.data.DataTable;
import de.erichseifert.gral.plots.XYPlot;
import de.erichseifert.gral.plots.areas.AreaRenderer;
import de.erichseifert.gral.plots.areas.DefaultAreaRenderer2D;
import de.erichseifert.gral.plots.lines.DefaultLineRenderer2D;
import de.erichseifert.gral.plots.lines.LineRenderer;
import de.erichseifert.gral.ui.InteractivePanel;
import de.erichseifert.gral.util.Insets2D;

public class EstimatorGUI extends JApplet implements ActionListener {

	private static final long serialVersionUID = 1L;
	
	// The option buttons
	private JButton startButton;
	private JButton stopButton;
	private JButton settingsButton;
	private JButton resetButton;
	private JPanel optionsPanel;
	private JTextField sampleLabel;
	
	// The settings panel
	SettingsUI SettingsFrame;	
	
	private XYPlot dataPlot;						 // The plot of the density
	private InteractivePanel dataPanel;				 // The panel storing the density plot
	private DataTable densityTable;                  // The density information
	
	// The window size
	private static final int WINDOW_WIDTH  = 600;
	private static final int WINDOW_HEIGHT = 600;	                   
	
	// Thread to perform calculations and update plot.
	private DensityRunner runner;
	
	
	
	/**
	 * Initializes the applet's GUI and Density Estimator variables.
	 */
	public void init() {
		
		// Set the size of the applet.
		setSize( WINDOW_WIDTH, WINDOW_HEIGHT );
	
                	
    	// Initialize algorithm variables.
    	try { Wavelet.init( Settings.waveletType ); } 
		catch ( Exception e ) { }
    	DensityHelper.initializeTranslates();
    	DensityHelper.initializeCoefficients();
    	
    	// Initialize the applet's GUI.
    	initializeGUI();
         
    }
	
	/**
	 * Creates the basic GUI components,
	 * the buttons and the data plot.
	 */
	public void initializeGUI() {
		
    	JPanel GUI = new JPanel();
    	GUI.setLayout( new BorderLayout() );
    	
    	// Set up the options panel
    	optionsPanel = new JPanel();
    	optionsPanel.setAlignmentX( 0 );
    	optionsPanel.setAlignmentY( 0 );
    	
    	// Create and add the buttons to the panel
    	Dimension btnSize = new Dimension( 85, 25 );
    	startButton = new JButton ( "Start" );
    	startButton.addActionListener( this );
    	startButton.setPreferredSize( btnSize );
    	stopButton = new JButton( "Stop" );
    	stopButton.setEnabled( false );
    	stopButton.setPreferredSize( btnSize );
    	stopButton.addActionListener( this );    	
    	settingsButton = new JButton ( "Settings" );
    	settingsButton.setPreferredSize( btnSize );
    	settingsButton.addActionListener( this );
    	resetButton = new JButton( "Reset" );
    	resetButton.addActionListener( this );
    	resetButton.setPreferredSize( btnSize );
    	sampleLabel = new JTextField();
    	sampleLabel.setText( "Sample index " );
    	sampleLabel.setEditable( false );
    	sampleLabel.setPreferredSize( new Dimension( 150, 25 ) );
    	sampleLabel.setHorizontalAlignment( JTextField.CENTER );
    	optionsPanel.add( sampleLabel );
    	optionsPanel.add( startButton );
    	optionsPanel.add( stopButton );
    	optionsPanel.add( resetButton );
    	optionsPanel.add( settingsButton );
    	GUI.add( optionsPanel, BorderLayout.NORTH );
        
        // Create the plot for the data
        dataPlot = new XYPlot();
        initializePlot();
        dataPanel = new InteractivePanel( dataPlot );
        GUI.add( dataPanel, BorderLayout.CENTER );
        
        // Create the settings UI
        SettingsFrame = new SettingsUI();
        SettingsFrame.setVisible( false );
        
        // Add the frame to the display
        add( GUI );
	}

	/**
	 * Responds to action events initiated by the user.
	 */
	public void actionPerformed( ActionEvent e ) {
		
		// User presses start button
		if ( e.getSource() == startButton ) {
			
			// Update the buttons' labels to reflect proper interaction.
			startButton.setEnabled( false );
			settingsButton.setEnabled( false );
			stopButton.setEnabled( true );
			
			// Initiate the density estimation process.
			startDensityEstimation();
	        
		} // if (e.getSource() == startButton).
		
		// User presses stop/resume button
		else if( e.getSource() == stopButton ){
			
			// Pause the Density Runner if the user
			// presses the Stop button.
			if( e.getActionCommand() == "Stop" )
			{
				// Update the Stop/Resume button's label appropriately.
				stopButton.setText( "Resume" );
				runner.pause(); // Pause the Density Runner.
				
			} // end if( e.getActionCommand() == "Stop" ).
			
			// Resume the Density Runner if the user
			// presses the Resume button.
			else
			{
				// Update the Stop/Resume button's label appropriately.
				stopButton.setText( "Stop" );
				runner.resume(); // Resume the Density Runner.
			} // end else -> if( e.getActionCommand() == "Stop" ).
			
		} // end else if( e.getSource() == stopButton ).
		
		// User presses the Settings button.
		else if ( e.getSource() == settingsButton ) {
			SettingsFrame.setVisible( true ); // show the Settings Panel.
		} // end else if(e.getSource() == settingsButton).
		
		// User presses the Reset button.
		else if( e.getSource() == resetButton )
		{
			runner.terminate(); // Terminate the Density Runner.
		} // end else if(e.getSource() == resetButton).
		
	} // end method public void actionPerformed(ActionEvent e).
	

	/**
	 * Initializes the plot of the probability density function to be estimated.
	 */
	private void initializePlot() {
		
		// Create the density function data.
    	densityTable = new DataTable( 2, Double.class );
        for ( double x = Settings.getMinimumRange(); 
        		x <= Settings.getMaximumRange() + Settings.discretization;
        		x += Settings.discretization ) {
            double y = 0.0;
            densityTable.add( x, y );
        }
        
        // Initialize the density function to a zero function.
        DensityHelper.updateDensity( densityTable );
        
        // Add the density function to the plot.
        dataPlot.add( densityTable );
        
        LineRenderer lines = new DefaultLineRenderer2D();
        AreaRenderer area = new DefaultAreaRenderer2D();
        dataPlot.setLineRenderer( densityTable, lines );  // Set line rendering to the plot.
        dataPlot.setAreaRenderer( densityTable, area );   // Set area-below-the-curve rendering to the plot.
        dataPlot.setInsets( new Insets2D.Double( 5, 50, 40, 40 ) ); // Add padding between the plot's frame and applet window.
        dataPlot.getAxisRenderer( XYPlot.AXIS_Y ).setIntersection( -Double.MAX_VALUE ); // Push the y-axis all the way to the left of the plot.
        Color invis = new Color( 0.0f, 0.3f, 1.0f, 0 ); 
        Color lineColor = new Color( 0.0f, 0.3f, 1.0f );
        Color areaColor = new Color( 0.0f, 0.3f, 1.0f, 0.3f );
        
        // Apply colors to the plot's lines, points, and area below the curve.
        dataPlot.getPointRenderer( densityTable ).setColor( invis ); 
        dataPlot.getLineRenderer( densityTable ).setColor( lineColor );
        dataPlot.getAreaRenderer( densityTable ).setColor( areaColor );
		
	} // end method private void initializePlot().
	

	/**
	 * Calls the DensityRunner.
	 */
	private void startDensityEstimation(){
		
		runner = new DensityRunner( sampleLabel, this.getWidth(), this.getHeight(), startButton, stopButton, settingsButton, dataPlot, dataPanel );
		runner.execute();
	} // end method private void startDensityEstimation().
	
} // end class EstimatorGUI.
