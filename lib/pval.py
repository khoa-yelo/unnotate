import numpy as np
from scipy.stats import norm

DEFAULT_X0 = 200
DEFAULT_MU = 0.985
DEFAULT_SIGMA = 0.005

def local_gauss(x0, x, y, h):
    """
    Compute local Gaussian mean & std at x0 using bandwidth h.
    """
    w = np.exp(-0.5*((x - x0)/h)**2)
    if np.sum(w) == 0:
        return DEFAULT_MU, DEFAULT_SIGMA
    mu = np.sum(w * y) / np.sum(w)
    var = np.sum(w * (y - mu)**2) / np.sum(w)
    return mu, np.sqrt(var)

def upper_p_value(x0, y0, x, y, h=0.5):
    """
    Returns the one‐tailed p‐value that Y > y0 under
    Y | X=x0 ~ N(mu0, sigma0^2) estimated locally.
    """
    if not x0:
        x0 = DEFAULT_X0
    mu0, sigma0 = local_gauss(x0, x, y, h)
    z = (y0 - mu0) / sigma0
    return 1 - norm.cdf(z)

def is_upper_outlier(x0, y0, x, y, h=0.5, alpha=0.05):
    """
    Returns True if y0 is significantly above the
    (1-alpha) upper bound at x0.
    """
    mu0, sigma0 = local_gauss(x0, x, y, h)
    z_crit = norm.ppf(1 - alpha)
    return (y0 - mu0) > z_crit * sigma0
   