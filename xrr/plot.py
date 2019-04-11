# -*- coding: utf-8 -*-
import matplotlib
import matplotlib.pyplot as plt
import numpy
import sympy
import scipy.special

from .layer import StackedLayer


matplotlib.rcParams["mathtext.fontset"] = "stixsans"
matplotlib.rcParams["font.size"] = 13
matplotlib.rcParams["text.usetex"] = True
matplotlib.rcParams['text.latex.preamble'] = [r"\usepackage[helvet]{sfmath}"]
matplotlib.rcParams["legend.fontsize"] = 13
matplotlib.rcParams["figure.figsize"] = (8, 5)


def plot_layer(ax, zdata, layer, fill, parent=None):
    if isinstance(layer, StackedLayer):
        for stack_layer in layer.stack:
            plot_layer(ax, zdata, stack_layer, fill, layer)
    else:    
        expr = layer.expr
        if parent is not None:
            expr = parent._substitute(expr)
        f_layer = sympy.lambdify(layer.z, expr,
                                 ("numpy", {"erf": scipy.special.erf}))
        ax.plot(zdata, f_layer(zdata), c=layer.color, lw=1, zorder=2)
        
    if parent is None and fill:
        f_layer = sympy.lambdify(layer.z, layer.expr,
                                 ("numpy", {"erf": scipy.special.erf}))
        ax.fill_between(zdata, 0, f_layer(zdata), color=layer.color, alpha=0.4,
                        zorder=0)


def plot_model(model, zmin, zmax, npoints=100, filename=None, dpi=150):
    fig = plt.figure()
    fig.patch.set_facecolor('white')
    ax = fig.gca()
    z = numpy.linspace(zmin, zmax, npoints)
    rho = model(z)
    ax.plot(z, rho, c="#424242", lw=1.5, zorder=3)
    for layer in model._layers:
        plot_layer(ax, z, layer, True)
    ax.set_xlim(zmin, zmax)
    ax.set_xlabel(r"z [\AA]")
    ax.set_ylabel(r"$\rho(z)/\rho_\infty$")
    if filename is None:
        plt.show()
    else:
        plt.savefig(filename, dpi=dpi)
