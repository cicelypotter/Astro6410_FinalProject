# Astro6410_FinalProject

## Redshift vs. g-r Color
The rate at which the universe is expanding is called the Hubble Flow. It can be measured by measuring the distances to galaxies and their receding velocities, as those further from us should be moving away more quickly. Currently, there are two measured values which are in disagreement by over 4σ. This code estimates the redshift value of a galaxy based on color magnitude measurements via Kth Nearest Neighbor Density Estimation. The results are then compared with the accepted redshift values.

## Method
First, ~2000 random galaxies are selected from [nsa_v0_1_2.fits](http://nsatlas.org/data), using np.random.seed(0) so that they are the same galaxies each time. The g and r magnitudes and the redshifts of these galaxies are stacked in an array, and the array is sigma clipped at 3σ. Then, the KNeighborsDensity function from astroML.density_estimation is fit to/trained on these galaxies, and evaluated on a 2000x2000 grid where the minimum and maximum for x and y are the minimum and maximum of the color magnitudes (g-r) and redshifts of these same galaxies, respectively.

The KNeighborsDensity function is set to the 'bayesian' method, and the number of neighbors to select (k) is optimized by minimizing the difference between the redshift estimation and accepted redshift value of 1000 other galaxies, whose redshifts are estimated the same way as the final 'science' set. This was done manually, but could probably be done with a loop.

Redshifts are estimated by selecting the closest color magnitude on the grid to that of each galaxy, finding the location of highest density on the grid corresponding to that color (thus the most likely place for a data point to be found), and in turn finding which redshift this location implies.

## Results
Unfortunately, my computer could not handle looping to optimize k or running the code with k>50 for the ~2000 galaxy training set. I suspect that the optimal number of neighbors is much higher; maybe around ~400, since I've had success with k=10 for a group of 50 galaxies. I intended to use the CHPC, and did spend a while watching the tutorial videos and reading the documentation. I was able to connect and get started, but not able to run python despite trying numerous different ways. I am sure I will figure it out, but didn't give myself enough time to learn it before the project was due. The value for k that I settled on was (), and the maximum difference between accepted and estimated redshift that this gave me for the 'science' galaxy set was ()
