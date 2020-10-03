Tethered particle motion (TPM) is a single molecule particle tracking technique where a bead is attached to a DNA tether and acts as a reporter of DNA dynamics (https://en.wikipedia.org/wiki/Tethered_particle_motion). Here, I have uploaded the code that is used to create a calibration plot from the raw data for different DNA tether sizes.

selectbeads63X.m reads the raw tethered particle motion data. It corrects the drift using using immobile beads and selects good tethered beads which are tethered to single DNA tethers and hence are symmetric. The outpput for each file is saved as data structure c in a .mat file. The analysis is repeated for all the raw data files (1 file for each filed of view).

collectALLdata.m combines all the outputs from all the field of views for a single DNA tether length in a single file. The output mat file has the time series of the x-y location of each tether bead. It also has the location of the raw data stored in the variable SOURCE.

exp_8s_63X_short.m runs on the output of collectALLdata.m to calculate all the numbers and plots histograms. Root mean square distance and mean square distance as reporter of the size of the DNA tether. 

rho_dist.m runs on analyzed8s.mat and creates a normalised histogram for all the tethers which are further used for clustering analysis. The output data was saved in the same file and a new excel file is also created.

rho_dist_selected.m is used to create a distribution plot of selected beads

rho_select_hist.m is used to create mean rho and mean rho square distributions of selected beads. The mean and standard deviation of rho average and rho square average of the combined traces are calculated and stored.

The above steps are repeated for experiments on different DNA tether lengths and the mean and standard deviation of rho average and rho square average from each experiments are used to plot a calibration curve.

