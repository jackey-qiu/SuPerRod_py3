"""
Alternate fitting routine with NLopt

using NLopt versio 2.5.0-cp37-win64
"""

import nlopt as nl
import numpy as np
#import pandas
import matplotlib.pyplot as plt
from time import time


def fit_model(model, xdata, ydata, weights=1, verbose=False, plot=True, tol=1e-3):
    fitfunc, fit_params = model.get_fit_func()
    p_init = [p.value for p in fit_params]
    p_min = [p.fit_min for p in fit_params]
    p_max = [p.fit_max for p in fit_params]
    
    # plot initial values
    if plot:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        plt.plot(xdata, ydata, 'ko', label='data')
        fit, = ax.plot(xdata, abs(fitfunc(xdata, *p_init)), 'r-', label='fit')
        print (*p_init)
        plt.xlabel(r'$q_\mathrm{z}\>[\>\mathrm{\AA}^{-1}]$')
        plt.ylabel(r'$R\>/\>R_\mathrm{F}$')
        plt.legend()
        plt.ion()
        plt.show()

    # objective function, grad=None for gradient free algrotithm
    def costfunc(p, grad):
        # calculate gradient numerically
        if grad.size > 0:
            for i in range(0, len(p)):
                mask = np.zeros(len(p))
                mask[i] = 1
                grad[i] = (costfunc(p + tol*mask, np.array([])) - costfunc(p, np.array([]))) / tol
        fx = abs(fitfunc(xdata, *p)) # abs: workaround for imag values in fitfunc
        r = ydata - fx # risiduals
        if plot: # update every iteration
            fit.set_ydata(fx)
            plt.draw()
            plt.pause(0.001)
        cost = sum((r / weights)**2)
        if verbose:
            print(cost)
        return cost

    # algorithm naming scheme:
    # AB_xxx
    # A = {L, G}   local, global
    # B = {D, N}   gradient-based, derivative-free
    opt = nl.opt(nl.GN_CRS2_LM, len(p_init)) #GN_DIRECT_L_RAND #GN_CRS2_LM #GN_DIRECT_L
    #nlopt_set_local_optimizer(nl.LD_MMA) 
    opt.set_lower_bounds(p_min)
    opt.set_upper_bounds(p_max)
    opt.set_min_objective(costfunc)
    opt.set_xtol_rel(tol)
    opt.set_ftol_rel(tol)
    t0 = time()
    p_opt = opt.optimize(p_init)
    p_err = np.zeros(len(p_init)) # no error without gradient
    t1 = time() - t0
    if verbose:
        print('optimization time = {:2.0f}h {:2.0f}m {:2.0f}s'.format(t1//3600, (t1%3600)//60, t1%60))
    for n in range(len(fit_params)):
        fit_params[n].value = p_opt[n]
        fit_params[n].fit_error = p_err[n]
        
    return p_opt, p_err, fitfunc
