# -*- coding: utf-8 -*-
import numpy
import sympy
import scipy.special
#import math

from .layer import Layer


class ElectronDensityModel(object):
    
    def __init__(self):
        super(ElectronDensityModel, self).__init__()
        self._layers = []
        
    def __call__(self, z):
        z_sym = sympy.Symbol("z_model")
        model = 0
        for layer in self._layers:
            model += layer.expr.subs(layer.z, z_sym)
        f = sympy.lambdify(z_sym, model, ("numpy", {"erf": scipy.special.erf}))
        return f(z)
        
    def get_master(self):
        qz = sympy.Symbol("qz", real=True)
        temp = 0
        for layer in self._layers:
            temp += layer.structure_factor(qz)
        master = abs(temp)**2
        return qz, master
        
    def get_fit_func(self):
        qz = sympy.Symbol("qz", real=True)

        fit_symbols, fit_params = [], []

        temp = 0
        for layer in self._layers:
            temp += layer.structure_factor(qz, fit=True)
            fit_params += filter(lambda x: x.fit, layer.params)
        master = abs(temp)**2

        fit_symbols = list(map(lambda x: x.symbol, fit_params))
        f = sympy.lambdify([qz] + fit_symbols, master,
                           ("numpy", {"erf": scipy.special.erf}))

        return f, fit_params
        
    def get_params(self):
        fit_params = []
        for layer in self._layers:
            fit_params += layer.params
        return fit_params
            
    def add(self, layer):
        if isinstance(layer, Layer):
            self._layers.append(layer)
            layer.set_name("layer_{0}".format(len(self._layers) - 1))
        else:
            raise TypeError()
