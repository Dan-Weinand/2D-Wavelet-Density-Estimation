1-D Streaming Density Estimation
===================

###Authors: Dan Weinand and Gedeon Nyengele
Funding for writing this applet was generously provided the Florida Institute of Technology's AMALTHEA REU program.  The authors would like to thank their graduate mentor Mark Moyou and adviser Dr. Adrian Peter. 

Based off of a Matlab implementation by Mark Moyou and Eddy Ihou. The algorithms used are from García-Treviño and Barria's "Online wavelet-based density estimation for non-stationary streaming data" and Wegman and Caudle's "Density estimation from streaming data using wavelets".

## Use
This applet is intended to perform density estimation on a 1-dimensional data set using wavelets. Users may upload their own csv data files containing a single column of values (and no header) using the settings menu. The main advantage of this applet is that it displays a changing density distribution for non-stationary density distributions.

The default data set (which is included in this repository) is a concatenation of a strongly skewed unimodal distribution for 3000 samples followed by 5000 samples from a 5-peaked "claw" distribution. This is included to demonstrate the ability of the algorithms to respond to changes in the distribution functions which are being approximated.

### General functionality
To upload your own data, select settings. Next, from the settings window in data settings select 'browse'. Then select your own csv file. You may then select from the various user parameters as desired. When you are done with the settings, click 'Save settings'. Then in the main screen, click 'start' to begin the density estimation process. The sample index will display how many samples have been read in when the plot was updated.

### Recommendations and settings

#### Data sets

The current default settings work well for normalized data but are not  necessarily limited to such data sets. However, if non-normalized data is to be used, we strongly recommend that the density domain setting be edited in order to properly accommodate the data being used.

#### Selecting a resolution level

The resolution to select is arguably the most important parameter. Higher resolutions give a density distribution which is more responsive to local  fluctuations but which is also more susceptible to noise. See Jansen, Malfait, and Bultheel's 1995 "Generalized Cross validation for wavelet thresholding" for a technique applicable to wavelet resolution thresholding in the case of stationary density distributions. Rigorously theoretically grounded methods for optimizing the resolution given a non-stationary distribution have not been proposed to the extent of our knowledge.


#### Deciding how to age the data

Many data sets are stationary with respect to their distribution function. It is typically a good idea to consider whether or not it makes sense to look at the data as non-stationary before trying to optimize aging parameters. Assuming the data is non-stationary, an important consideration is how to implement aging. If exponential aging like that proposed by Caudle and Wegman 2009 is used, the optimal value of theta needs to be found. As discussed by García-Treviño and Barria 2012, selecting the optimal theta can be an exceedingly difficult process rather opaque to the end user. In general, increasing theta reduces the model's sensitivity to noise and local changes in the underlying density function while decreasing the parameter has the reverse effect. While a value of theta near .999 has seen some success empirically, it is difficult to understand how much to modify the parameter in order fine tune the estimate's sensitivity to both noise and local changes in the nature of the distribution function. In contrast, using the sliding window approach generally feels much more intuitive. Since the technique simply uses a sliding average of W samples, the user may pick how many samples to average at a time and the trade off between noise insensitivity versus rapid adjustment to changes in the distribution function is fairly straightforward. Based on this, we recommend the sliding window approach unless there is a clear-cut reason for using exponential discounting.

#### Other considerations

Enabling wavelets at a stopping resolution level of J is exactly equivalent to disabling wavelets at a starting resolution of (J + 1) in terms of the density distribution function. The wavelet type and order to be selected are relatively unimportant in actually finding the density estimate, although a general rule is that higher orders of a given wavelet family will give smoother estimates. If wavelets are to be used, make sure to select a wavelet type which is dyadic in nature (the Daubechies, Symlet or Coiflet wavelets are all dyadic).

There are several factors which influence the speed of the density estimation process. The discretization level determines how well the plot shown demonstrates the actual calculated density distribution. Reducing the level will increase the speed of the algorithm but may result in a plot which appears to be somewhat piece-wise linear and jagged. The update frequency influences the algorithm's speed both by determining how often the plot is visually updating (repainting the plot can easily begin to dominate computation for very frequent updates) and how often the density must be updated. Broadly speaking, if speed is a concern one should update as rarely as feasible. Finally, higher order wavelets have larger tables defined and thus computing phi and psi will take longer to interpolate. For this reason mid-range Daubechies wavelets are recommended as generally a good balance between smoothness and speed.

## Background information
There is a wealth of information available to help the interested user familiarize themselves with density estimators, dealing with data streams, and the theory of wavelets.  Some useful resources are given below.

### Wavelets
[A Gentle Introduction to Wavelets](http://web.media.mit.edu/~rehmi/wavelet/wavelet.html)

[A Practical Guide to Wavelet Analysis](http://paos.colorado.edu/research/wavelets/)

[Surfing the Wavelets](http://www.wavelet.org/tutorial/)

[An introduction to wavelets](http://www.eecis.udel.edu/~amer/CISC651/IEEEwavelet.pdf)

### Density Estimation
[Density Estimation for statistics and data analysis](http://ned.ipac.caltech.edu/level5/March02/Silverman/paper.pdf)

[Nonparametric density function estimation](http://igpphome.ucsd.edu/~cathy/Classes/SIO223A/sio223a.chap9.pdf)

[Nonparametric density estimation using wavelets](http://stat.duke.edu/sites/default/files/papers/1993-26.pdf)

### Streaming data
[Density estimation over data stream](http://alumni.cs.ucr.edu/~wli/publications/deosd.pdf)

[Online wavelet-based density estimation for non-stationary streaming data](http://www.sciencedirect.com/science/article/pii/S0167947311003082)

