package edu.fit.estimator2D;

/**
 * Class to both display the current density estimate,
 * and allow the user to select various parameters.
 * 
 * @author Daniel Weinand & Gedeon Nyengele
 * 
 */

import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JApplet;
import javax.swing.JButton;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextField;

import org.sf.surfaceplot.SurfaceCanvas;

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
	
	// The window size
	private static final int WINDOW_WIDTH  = 900;
	private static final int WINDOW_HEIGHT = 600;	
	
	// Plot containers.
	private SurfaceCanvas canvas = new SurfaceCanvas();
	private JPanel plotPanel;
	private JPanel commandPanel;
	private JLabel rotationLabel;
	private JLabel zoomLabel;
	private JLabel moveLabel;
	
	// Thread to perform calculations and update plot.
	private DensityRunner runner;
	
	// Data model.
	DensityModel dataModel;
	
	
	
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
    	
    	// Create the data model.
    	dataModel = new DensityModel( Settings.discretization, Settings.densityRange );
    	
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
    	resetButton.setEnabled( false );
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
        plotPanel     = new JPanel();
        commandPanel  = new JPanel();
        plotPanel.setLayout( new BorderLayout() );
        commandPanel.setLayout(new FlowLayout(FlowLayout.CENTER, 50, 5));
        
        rotationLabel = new JLabel( "Rotate: Mouse Click & Drag" );
        zoomLabel = new JLabel( "Zoom: Shift Key + Mouse Click & Drag" );
        moveLabel = new JLabel( "Move: Control Key + Mouse Click & Drag" );
        
        commandPanel.add(rotationLabel);
        commandPanel.add(zoomLabel);
        commandPanel.add(moveLabel);
        
        plotPanel.add(canvas, BorderLayout.CENTER);
        initializePlot();
        canvas.setModel(dataModel);
        
        plotPanel.add(commandPanel, BorderLayout.SOUTH);       
        GUI.add( plotPanel, BorderLayout.CENTER );
        
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
			resetButton.setEnabled( true );
			
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
		
		float incX = (float)(dataModel.getXMax() / dataModel.getCalcDivisions());
		float incY = (float)(dataModel.getYMax() / dataModel.getCalcDivisions());
		int iMax = (int)(dataModel.getXMax() / incX);
		int jMax = (int)(dataModel.getYMax() / incY);
		double[][] zvals = new double[iMax+1][jMax+1];
		for(int i = 0; i <= iMax; i++)
		{
			for(int j = 0; j <= jMax; j++)
			{
				zvals[i][j] = 0;
			}
		}
		
		DensityModel.updateDensity(zvals, 1.5f);
		
	} // end method private void initializePlot().
	

	/**
	 * Calls the DensityRunner.
	 */
	private void startDensityEstimation(){
		
		runner = new DensityRunner( sampleLabel, startButton, stopButton, settingsButton, canvas, dataModel );
		runner.execute();

	} // end method private void startDensityEstimation().
	
} // end class EstimatorGUI.
