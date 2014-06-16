import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.IOException;

import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JSeparator;
import javax.swing.JTextField;
import javax.swing.border.Border;

/**
 * Allows the user to set/select settings.
 * @author Gedeon Nyengele & Daniel Weinand
 *
 */
public class SettingsUI extends JFrame{

	private static final long serialVersionUID = 1L;
	private JLabel startLevelLabel, stopLevelLabel, agingFlagLabel, discretizationLabel, windowSizeLabel, agingThetaLabel;
	private JLabel waveletFlagLabel, fileNameLabel, updateFrequencyLabel;
	private JLabel densityFrom, densityTo, waveletTypeLabel;
	private JTextField startLevel, stopLevel, discretization, windowSize, agingTheta, updateFrequency;
	private JTextField densityRangeFrom, densityRangeTo;
	private JComboBox<String> agingFlag, waveletType;
	private JPanel waveletPanel, plotPanel, dataPanel, densityRangePanel;
	private JCheckBox waveletFlag;
	private JFileChooser fileChooser;
	private JButton btnOpenFile, btnSaveSettings;
	private JPanel content;
	
	
	/**
	 * Overloaded constructor.
	 * @param width  : width of the applet.
	 * @param height : height of the applet.
	 */
	public SettingsUI(int width, int height){
		super("Estimator Settings");		
		this.setSize(width, height);
		content = (JPanel) this.getContentPane();
		init();
		
		
	}
	
	/**
	 * Default constructor.
	 */
	public SettingsUI()
	{
		super("Estimator Settings");
		setSize(300, 480);
		content = (JPanel) this.getContentPane();
		init();
	}
	
	/**
	 * Initializes all the GUI components and
	 * add action listener on the save button.
	 */
	public void init()
	{
		this.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
		this.setLayout(new FlowLayout());
		this.setResizable(false);
		this.setLocationRelativeTo(null);
		
		// Initialize labels.
		startLevelLabel      = new JLabel("Resolution Start Level: ");
		stopLevelLabel       = new JLabel("Resolution Stop Level: ");
		agingFlagLabel       = new JLabel("Aging Flag: ");
		discretizationLabel  = new JLabel("Discretization: ");
		windowSizeLabel      = new JLabel("Window Size: ");
		agingThetaLabel      = new JLabel("Aging Theta: ");
		waveletFlagLabel     = new JLabel("Wavelet Flag: ");
		updateFrequencyLabel = new JLabel("Plot Update Frequency: ");
		fileNameLabel        = new JLabel("Sample Data File: ");
		densityFrom          = new JLabel("From: ");
		densityTo            = new JLabel("To: ");
		waveletTypeLabel     = new JLabel("Wavelet Type: ");
		
		// Initialize textFields.
		startLevel           = new JTextField(Settings.startLevel + "" , 3);
		stopLevel			 = new JTextField(Settings.stopLevel + "" , 3);
		discretization       = new JTextField(Settings.discretization + "" , 5);
		windowSize			 = new JTextField(Settings.windowSize + "" , 10);
		agingTheta           = new JTextField(Settings.agingTheta + "", 5);
		updateFrequency      = new JTextField(Settings.updateFrequency + "", 10);
		densityRangeFrom     = new JTextField(Settings.densityRange[0] + "", 5);
		densityRangeTo       = new JTextField(Settings.densityRange[1] + "", 5);
		
		// Initialize combo boxes.
		String[] agingFlags  = {"No Aging", "Caudle", "Window"};
		agingFlag            = new JComboBox<String>( agingFlags );
		agingFlag.setSelectedIndex(Settings.windowAge);  // sets the default to window aging. windowAging is constant 2 in Settings class.
		waveletType          = new JComboBox<String>( Settings.waveletTypes );
		waveletType.setSelectedIndex(9);
		
		// Initialize check Boxes.
		waveletFlag          = new JCheckBox("Enable wavelet.");
		
		// Initialize buttons.
		btnOpenFile          = new JButton("Browse...");
		btnSaveSettings      = new JButton("Save settings");
		btnSaveSettings.addActionListener(new ActionListener(){
			public void actionPerformed(ActionEvent e){
				saveSettings();
			}
		});
		btnOpenFile.addActionListener(new ActionListener(){
			public void actionPerformed(ActionEvent e){
				fileChooser = new JFileChooser();
				fileChooser.setCurrentDirectory( new java.io.File(".") );
				fileChooser.setDialogTitle("Estimator Sample Data File");
				fileChooser.setFileSelectionMode(JFileChooser.FILES_ONLY);
				fileChooser.setAcceptAllFileFilterUsed(false);
				
				if(fileChooser.showOpenDialog(null) == JFileChooser.APPROVE_OPTION)
				{
					java.io.File file = fileChooser.getSelectedFile();
					try {
						Settings.dataFile = file.getCanonicalPath();
					} catch (IOException e1) {}
					
				}
			}
		});
		
		
		// Initialize panels.
		waveletPanel         = new JPanel( new FlowLayout() );		
		plotPanel            = new JPanel( new FlowLayout() );
		dataPanel            = new JPanel( new FlowLayout() );
		densityRangePanel    = new JPanel( new FlowLayout() );
		Border panelBorder   = BorderFactory.createLineBorder(Color.BLACK, 2);
		waveletPanel.setBorder( BorderFactory.createTitledBorder(panelBorder, "Wavelet Settings" ));
		plotPanel.setBorder(BorderFactory.createTitledBorder(panelBorder, "Plot Settings" ));
		dataPanel.setBorder(BorderFactory.createTitledBorder(panelBorder, "Data Settings" ));
		densityRangePanel.setBorder(BorderFactory.createTitledBorder(panelBorder, "Density Range" ));
		
		// Add waveletPanel components.
		GridLayout wavPanelLayout = new GridLayout(0, 2);
		waveletPanel.setLayout(wavPanelLayout);
		waveletPanel.add(startLevelLabel);        waveletPanel.add(startLevel);
		waveletPanel.add(stopLevelLabel);         waveletPanel.add(stopLevel);
		waveletPanel.add(agingFlagLabel);         waveletPanel.add(agingFlag);
		waveletPanel.add(discretizationLabel);    waveletPanel.add(discretization);
		waveletPanel.add(windowSizeLabel);        waveletPanel.add(windowSize);
		waveletPanel.add(agingThetaLabel);        waveletPanel.add(agingTheta);   
		waveletPanel.add(waveletTypeLabel);       waveletPanel.add(waveletType);
		waveletPanel.add(waveletFlagLabel);       waveletPanel.add(waveletFlag);
		
		
		// Add plotPanel components.
		plotPanel.setPreferredSize( new Dimension(270, 55) );
		plotPanel.add(updateFrequencyLabel);      plotPanel.add(updateFrequency);
		
		// Add dataPanel components.
		dataPanel.setPreferredSize( new Dimension(270, 60) );
		dataPanel.add(fileNameLabel);
		dataPanel.add(btnOpenFile);
		
		// Add densityRangePanel components.
		densityRangePanel.setPreferredSize( new Dimension(270, 55) );
		densityRangePanel.add(densityFrom);	
		densityRangePanel.add(densityRangeFrom);
		densityRangePanel.add(densityTo);	
		densityRangePanel.add(densityRangeTo);
		
		// Add the Panels to the frame.
		content.add(waveletPanel, BorderLayout.LINE_START);
		content.add(plotPanel, BorderLayout.LINE_START);
		content.add(dataPanel, BorderLayout.LINE_START);
		content.add(densityRangePanel, BorderLayout.LINE_START);
		content.add(new JSeparator(JSeparator.HORIZONTAL));
		content.add(btnSaveSettings);
		
		
	}
	
	private void saveSettings(){
		
		
		// Process all integer-valued settings.
		Settings.startLevel        = Integer.parseInt(startLevel.getText());
		Settings.stopLevel         = Integer.parseInt(stopLevel.getText());
		Settings.windowSize        = Integer.parseInt(windowSize.getText());
		Settings.updateFrequency   = Integer.parseInt(updateFrequency.getText());
		
		// Process all real-valued settings.
		Settings.discretization    = Double.parseDouble(discretization.getText());
		Settings.agingTheta        = Double.parseDouble(agingTheta.getText());
		Settings.densityRange[0]   = Double.parseDouble(densityRangeFrom.getText());
		Settings.densityRange[1]   = Double.parseDouble(densityRangeTo.getText());
		
		// Process checkboxes.
		Settings.waveletFlag       = waveletFlag.isSelected();
		
		// Process comboboxes.
		Settings.agingFlag         = agingFlag.getSelectedIndex();
		
		// Process string-valued settings.
		Settings.waveletType       = Settings.waveletTypes[waveletType.getSelectedIndex()];
		
		// Close the Settings window.
		setVisible(false);
		dispose();
	}
	
	
	

}
