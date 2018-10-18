
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

# 4. Markowitz Efficient Frontier and Allocation

In `test_markowitz.R` I am developing the plots and weights for the efficient frontier with and without Bitcoin, then the graphs of the allocations. The functions needed to do so are in `MarkowitzMeanVariancePortfolio.R`. Both constrained (=no short-selling) and unconstrained frontiers are implemented.

I also computed the efficient frontier using the CVaR (alternatively the VaR) as the measure of risk in `test_empirical_VaR_optimization.R`, implementing two different approaches to such computations. 

The first is a "bootstrap" technique, where given the sample of daily returns, I built their empirical distribution and then extract 255
values from it to obtain a year evolution and thus a yearly return. Repeating this step N=10000 times I obtained N different scenarios which are then used to compute the portfolio yearly returns and finally the VaR/CVaR.

The second approach focuses on the daily return and given the last 5 years of observations (approximately 255*5=1275 vectors of daily returns)  the portfolio VaR/CVaR is computed.

Results from each approach are in the `Results` folder.
