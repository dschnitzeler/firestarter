# firestarter
--- Program description in a nutshell ---

Firestarter is a QU-fitting based algorithm for extracting information from polarization data. How the algorithm works is described in Schnitzeler, 2018, MNRAS, or ‘DS18’. 

In short, the program maximizes the likelihood to find the best-fitting parameters and their uncertainties, based on your measurements and the (1-sigma) uncertainties on these measurements. The program assumes that each source intrinsically emits a power-law flux density spectrum in Stokes Q and U (intrinsically = before this emission is depolarized), then it fits for the spectral index of this power law, and for the values of Stokes Q and U associated with this power law. The program assumes that the noise in Q and U, and the frequency channels, are statistically independent. Noise variances in Q and U do not have to be the same in each channel, and are allowed to vary across the frequency band. 

Typically, the user provides the core program ‘Fitit’ in fs.pro with these measurements, a list of source types that can be fitted, and the maximum number of source components that should be fitted. The fitting program then creates a list of all the models that should be fitted automatically, and calculates metrics like the Akaike Information Criterion and the Bayesian Information Criterion for each fit to the data. When all models have been fitted, they are ranked automatically for each of these metrics. The significance of a detection is quantified using the log likelihood ratio, and the signal-to-noise ratio by using the amplitude of the polarized signal and the covariance matrix of the noise in Stokes Q and U.
 


--- Installing the software ---

The software package consists of four programs and a library. Please install these all in the same subdirectory, or add the programs in ‘lib’ to your IDL or GDL path. In addition, you should install the CMSVLIB and MPFIT libraries (both written by Craig Markwardt). These can be downloaded from 

  http://cow.physics.wisc.edu/~craigm/idl/cmsave.html

  http://cow.physics.wisc.edu/~craigm/idl/fitting.html



--- Quick tour ---

The are four core programs:

fs.pro > Fitit: the main fitting program.

  Example: Fitit, meas_arr=meas_arr, cmp_arr=1, n_cmp_max=5

readit.pro > Readit: use this program to interpret the output from fs.pro, and to plot the data and the best-fitting model fit, as selected by the Bayesian Information Criterion.

  Example: Readit, meas_arr=meas_arr

batch_fitit.pro > Batch_fitit: calls ‘Fitit’ repeatedly, for example, when running a Monte Carlo simulation.

  Example: Batch_fitit, test=101

batch_readit.pro > Batch_readit: calculates model weights, and model-weighted parameters. Useful when running Monte Carlo simulations.

  Example: Batch_readit, test=101

I would start with running Fitit and Readit. Use Batch_fitit and Batch_readit if you want to generate mock observations and apply Fitit to each of these simulations. Information on how to set up your own Monte Carlo simulation can be found in the section ‘How to set up a Monte Carlo simulation’ at the beginning of batch_fitit.pro.

More information on how to call the programs and their output can be found at the start of each program. Alternatively, you can type ‘fitit’, ‘batch_fitit’, ‘readit’, or ‘batch_readit’ after compiling each of the core programs.



--- Request ---

If you find a bug, or even fix a big, please let me know so that I can update this repository. You can contact me at d{insert my surname}(at){name of that e-mail service provided by Google}.com.
