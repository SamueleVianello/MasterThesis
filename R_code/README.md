
# 1. Multivariate Merton Model
Files `MultivariateMertonModel.R` and  `CalibrationMVMertonModel.R` implement a library that computes density, loglikelihood and calibrates a Multivariate Merton model for the evolution of n-assets.

The efficient implementation only goes up to 4 assets at the same time, and still the calibration of the model may be quite time consuming.

Anyway, I introduced a procedure to obtain any NxN covariance matrix for the diffusion part since that is what I am interested in.
This procedure can be improved (so far I only compute the mean on the components that are calibrated more than once) by fixing the parameters that have already been compute and calibrating only on the remaing free coefficients.
