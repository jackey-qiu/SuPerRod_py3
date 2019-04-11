"""
Routines to fit a model to measured data.

Requires scipy version 0.17 or higher.
"""
# exceptions etc.
from . import XRRLibraryVersionError
from distutils.version import LooseVersion
# external libraries
import scipy
if LooseVersion(scipy.__version__) < LooseVersion("0.17"):
    raise XRRLibraryVersionError("scipy", "0.17", scipy.__version__)
from scipy.optimize import curve_fit
import numpy


def fit_model(model, xdata, ydata, weights=1):
    # make fit func
    fit_func, fit_params = model.get_fit_func()
    # prepare init values and bounds
    p_init = [p.value for p in fit_params]
    p_min = [p.fit_min for p in fit_params]
    p_max = [p.fit_max for p in fit_params]
    p_bounds = (p_min, p_max)
    print(p_init)
    # do the fitting
    p_opt, p_cov = curve_fit(fit_func, xdata, ydata, p_init, sigma=weights,
                             bounds=p_bounds,method="trf", max_nfev=100000000) #max_nfev=10000000
    p_err = numpy.sqrt(numpy.diag(p_cov))
    # set layer paramters to the obtained values
    for n in range(len(fit_params)):
        fit_params[n].value = p_opt[n]
        fit_params[n].fit_error = p_err[n]
    return p_opt, p_err, fit_func
