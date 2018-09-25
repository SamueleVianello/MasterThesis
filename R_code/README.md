
# 1. Multivariate Merton Model
Files `MultivariateMertonModel.R` and  `CalibrationMVMertonModel.R` implement a library that computes density, loglikelihood and calibrates a Multivariate Merton model for the evolution of n-assets.

The efficient implementation only goes up to 4 assets at the same time, and still the calibration of the model may be quite time consuming depending on the number of observations you have (with about 2100 obs of 4 assets, it takes me around 3 minutes).

I introduced a procedure to obtain any NxN covariance matrix for the diffusion part through sequential computations of a 4-asset model.
This procedure can still be improved (so far I only compute the mean on the components that are calibrated more than once) by fixing the parameters that have already been compute and calibrating only on the remaing free coefficients.

All the computations on the data for this part is done in the file `main_calibration.R`.

# 2. Correlation Significance

Studies on the correlation between Bitcoin and other assets and its significance are performed in `main_correlation.significance.R`.
I computed the sample correlations of the log-returns and then performed three tests (Pearson's, Spearman's and a Permutation test) to check whether their values are significantly different from zero.

The results clearly show that the historical correlation between Bitcoin and the other assets taken into consideration is *not significantly different* from zero.

# 3. Rolling Correlation

In `main_rolling_correlation.R` I computed the rolling correlations and plotted the relative graphs. I included in each plot also the p-value for their significances.



