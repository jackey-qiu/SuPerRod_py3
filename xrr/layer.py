# -*- coding: utf-8 -*-
import numpy
import scipy
import sympy
import uuid


class LayerParameter(object):
    
    def __init__(self, value, symbol, fit=False, fit_min=-numpy.inf,
                 fit_max=numpy.inf):
        super(LayerParameter, self).__init__()
        self.symbol = symbol
        self.value = value
        self._fit = fit
        self._min = fit_min
        self._max = fit_max
        self._fit_err = None
        self._symbol_name = symbol.name
        self._name = symbol.name
        self.symbol.name = str(uuid.uuid4())
        
    @property
    def fit(self):
        return self._fit
        
    @fit.setter
    def fit(self, do_fit):
        self._fit = do_fit
        
    @property
    def fit_min(self):
        return self._min
    
    @fit_min.setter
    def fit_min(self, fit_min):
        self._min = fit_min
        
    @property
    def fit_max(self):
        return self._max
    
    @fit_max.setter
    def fit_max(self, fit_max):
        self._max = fit_max
        
    @property
    def fit_error(self):
        return self._fit_err
        
    @fit_error.setter
    def fit_error(self, err):
        self._fit_err = err
        
    @property
    def name(self):
        return self._name
        
    @name.setter
    def name(self, name):
        self._name = name

            
class Layer(object):
    
    def __init__(self, density=1.0, color="#A3A3A3"):
        super(Layer, self).__init__()
        self.z = sympy.Symbol("z", real=True)
        self._s_density = sympy.Symbol("rho", real=True, positive=True)
        self.density = LayerParameter(density, self._s_density)
        self._master = None
        self._master_fit = None
        self._term = sympy.numbers.Zero()
        self.color = color
        self.name = "layer"
        
    def __call__(self, z):
        z_sym = sympy.Symbol("z_model")
        expr = self.expr.subs(self.z, z_sym)
        f = sympy.lambdify(z_sym, expr, ("numpy", {"erf": scipy.special.erf}))
        return f(z)
        
    def get_expr(self, fit=False):
        return self._substitute(self._term, fit=fit)
        
    expr = property(get_expr)
    
    def set_name(self, name):
        self.name = name
        self.set_param_names()
        
    def set_param_names(self):
        for param in self.params:
            param.name = "{0}_{1}".format(self.name, param._symbol_name)
    
    @property
    def params(self):
        # find layer parameters
        attrs = list(self.__dict__.keys())
        def f(x):
            return (not x.startswith("_") and
                    isinstance(self.__dict__[x], LayerParameter))
        p_attrs = list(filter(f, attrs))
        params = list(map(lambda x: self.__dict__[x], p_attrs))
        return params        
        
    def _substitute(self, expr, fit=False):
        # substitute parameters
        params = self.params
        if fit:
            # sort out fit parameters
            params = list(filter(lambda x: not x.fit, self.params))
        subs = list(map(lambda x: (x.symbol, x.value), params))
        return expr.subs(subs)
        
    def get_fit_parameters(self):
        params = self.params
        # get fit parameters
        params = list(filter(lambda x: x.fit, self.params))
        return params
        
    def structure_factor(self, qz, fit=False):
        if not fit and self._master is None:
            integrand = sympy.diff(self._term, self.z) * sympy.exp(1j*qz*self.z)
            master = sympy.integrate(integrand, (self.z, -sympy.oo, sympy.oo))
            self._master = qz, self._substitute(master, False)
            return self._master[1]
        elif fit and self._master_fit is None:
            integrand = sympy.diff(self._term, self.z) * sympy.exp(1j*qz*self.z)
            master = sympy.integrate(integrand, (self.z, -sympy.oo, sympy.oo))
            self._master_fit = qz, self._substitute(master, True)
            return self._master_fit[1]
        elif not fit and self._master is not None:
            old_qz, old_master = self._master
            master = old_master.subs(old_qz, qz)
            self._master = qz, master
            return self._master[1]
        elif fit and self._master_fit is not None:
            old_qz, old_master = self._master_fit
            master = old_master.subs(old_qz, qz)
            self._master_fit = qz, master
            return self._master_fit[1]
            
            
class StackedLayer(Layer):
    
    def __init__(self, **kwargs):
        super(StackedLayer, self).__init__(**kwargs)
        self.stack = []
        
    def get_expr(self, fit=False):
        term = sympy.numbers.Zero()
        for layer in self.stack:
            term += layer.get_expr(fit=fit)
        return self._substitute(term, fit=fit)
    expr = property(get_expr)
    
    def set_name(self, name):
        super(StackedLayer, self).set_name(name)
        for d in enumerate(self.stack):
            n, layer = d
            layer_name = "{0}:{1}".format(name, n)
            layer.set_name(layer_name)
        

class GaussianLayer2(StackedLayer):
      
    def __init__(self, density, mu, sigma, **kwargs):
        super(GaussianLayer2, self).__init__(**kwargs)
        
        
        self.rho_g=sympy.Symbol("r_gau", real=True, positive=True)
        self.sigma_g=sympy.Symbol("sigma_gau", real=True, positive=True)
        self.z_g=sympy.Symbol("z_gau", real=True)
        
        self.density = LayerParameter(density, self.rho_g)
        self.sigma = LayerParameter(sigma, self.sigma_g)
        self.mu = LayerParameter(mu, self.z_g)
        
        self._make_term()
    
    
    def _make_term(self):
        
        z_g= self.z_g
        sigma_g= self.sigma_g
        amp=  self.rho_g #/(self.sigma_g *  sympy.sqrt(2 * sympy.pi) )
        layer= GaussianLayer(amp, z_g, sigma_g)
        
        self.stack.append(layer)
       
        #self.layer = (self.rho_g / (self.sigma_g *  (sympy.sqrt(sympy.pi)))* sympy.exp(-(self.z - self.z_g)**2 / 
        #              (2 * self.sigma_g**2)))
        
        
    def structure_factor(self, qz, fit=False):
        g = 1j* qz * self.rho_g * sympy.exp(1j * qz * self.z_g - (qz**2 * self.sigma_g**2) / 2.0 )
        gau = self._substitute(g, fit=fit)
        return gau
    
class GaussianLayer(Layer):
    
    c_master = None
    c_density = sympy.Symbol("c_density", real=True, positive=True)
    c_mu = sympy.Symbol("c_mu", real=True)
    c_sigma = sympy.Symbol("c_sigma", real=True, positive=True)
    c_z = sympy.Symbol("c_z", real=True)
    c_qz = sympy.Symbol("c_qz", real=True)
    
    def __init__(self, density, mu, sigma, **kwargs):
        super(GaussianLayer, self).__init__(density=density, **kwargs)
        self._s_mu = sympy.Symbol("mu", real=True)
        self._s_sigma = sympy.Symbol("sigma", real=True, positive=True)
        self.mu = LayerParameter(mu, self._s_mu)
        self.sigma = LayerParameter(sigma, self._s_sigma)
        self._term = (self._s_density * sympy.exp(-(self.z - self._s_mu)**2 / 
                      (2 * self._s_sigma**2)))
        
    def structure_factor(self, qz, fit=False):
        if self.c_master is None:
            expr = (self.c_density / (self.c_sigma *(sympy.sqrt(2*sympy.pi)))*sympy.exp(-(self.c_z - self.c_mu)**2 / 
                    (2 * self.c_sigma**2)))
            integrand = (-1*sympy.diff(expr, self.c_z) *
                         sympy.exp(1j*self.c_qz*self.c_z))
            master = sympy.integrate(integrand,
                                     (self.c_z, -sympy.oo, sympy.oo))
            self.__class__.c_master = master
        subs = [(self.c_qz, qz),
                (self.c_z, self.z),
                (self.c_density, self._s_density),
                (self.c_mu, self._s_mu),
                (self.c_sigma, self._s_sigma)]
        master_expr = self.c_master.subs(subs)
        master = self._substitute(master_expr, fit=fit)
        return master    
   
        
class GaussianLayer3(Layer):
    
    c_master = None
    c_density = sympy.Symbol("c_density", real=True, positive=True)
    c_mu = sympy.Symbol("c_mu", real=True)
    c_sigma = sympy.Symbol("c_sigma", real=True, positive=True)
    c_z = sympy.Symbol("c_z", real=True)
    c_qz = sympy.Symbol("c_qz", real=True)
    
    def __init__(self, density, mu, sigma, **kwargs):
        super(GaussianLayer3, self).__init__(density=density, **kwargs)
        self._s_mu = sympy.Symbol("mu", real=True)
        self._s_sigma = sympy.Symbol("sigma", real=True, positive=True)
        self.mu = LayerParameter(mu, self._s_mu)
        self.sigma = LayerParameter(sigma, self._s_sigma)
        self._term = (self._s_density * sympy.exp(-(self.z - self._s_mu)**2 / 
                      (2 * self._s_sigma**2)))
        
    def structure_factor(self, qz, fit=False):
        if self.c_master is None:
            expr = (self.c_density / (self.c_sigma *(sympy.sqrt(sympy.pi)))*sympy.exp(-(self.c_z - self.c_mu)**2 / 
                    (2 * self.c_sigma**2)))
            integrand = (-1*sympy.diff(expr, self.c_z) *
                         sympy.exp(1j*self.c_qz*self.c_z))
            master = sympy.integrate(integrand,
                                     (self.c_z, -sympy.oo, sympy.oo))
            self.__class__.c_master = master
        subs = [(self.c_qz, qz),
                (self.c_z, self.z),
                (self.c_density, self._s_density),
                (self.c_mu, self._s_mu),
                (self.c_sigma, self._s_sigma)]
        master_expr = self.c_master.subs(subs)
        master = self._substitute(master_expr, fit=fit)
        return master    
        

class GaussianLayer_lead(Layer):
    
    c_master = None
    c_density = sympy.Symbol("c_density", real=True, positive=True)
    c_mu = sympy.Symbol("c_mu", real=True)
    c_sigma = sympy.Symbol("c_sigma", real=True, positive=True)
    c_z = sympy.Symbol("c_z", real=True)
    c_qz = sympy.Symbol("c_qz", real=True)
    
    def __init__(self, density, mu, sigma, **kwargs):
        super(GaussianLayer_lead, self).__init__(density=density, **kwargs)
        self._s_mu = sympy.Symbol("mu", real=True)
        self._s_sigma = sympy.Symbol("sigma", real=True, positive=True)
        self.mu = LayerParameter(mu, self._s_mu)
        self.sigma = LayerParameter(sigma, self._s_sigma)
        self._term = (self._s_density * sympy.exp(-(self.z - self._s_mu)**2 / 
                      (2 * self._s_sigma**2)))
        
    def structure_factor(self, qz, fit=False):
        if self.c_master is None:
            expr = (self.c_density / (self.c_sigma *(sympy.sqrt(sympy.pi)))*sympy.exp(-(self.c_z - self.c_mu)**2 / 
                    (2 * self.c_sigma**2)))
            integrand = (sympy.diff(expr, self.c_z) *
                         sympy.exp(1j*self.c_qz*self.c_z))
            Pb_Zeff = 80 - (80/82.5)**2.37
            S = (13.4118 +
                 31.0617 * sympy.exp(-0.6902*qz**2/(16*sympy.pi**2)) +
                 13.0637 * sympy.exp(-2.3576*qz**2/(16*sympy.pi**2)) +
                 18.442 * sympy.exp(-8.618*qz**2/(16*sympy.pi**2)) +
                 5.9696 * sympy.exp(-47.2579*qz**2/(16*sympy.pi**2)))
            f = (S + 0.077) / (Pb_Zeff + 0.077)
            master = f *  sympy.integrate(integrand,
                                     (self.c_z, -sympy.oo, sympy.oo))
            self.__class__.c_master = master
        subs = [(self.c_qz, qz),
                (self.c_z, self.z),
                (self.c_density, self._s_density),
                (self.c_mu, self._s_mu),
                (self.c_sigma, self._s_sigma)]
        master_expr = self.c_master.subs(subs)
        master = self._substitute(master_expr, fit=fit)
        return master
        




class ErrorLayer(Layer):
    
    def __init__(self, density, mu, sigma, **kwargs):
        super(ErrorLayer, self).__init__(**kwargs)
        self.density.value = density
        self._s_mu = sympy.Symbol("mu", real=True)
        self._s_sigma = sympy.Symbol("sigma", real=True, positive=True)
        self.mu = LayerParameter(mu, self._s_mu)
        self.sigma = LayerParameter(sigma, self._s_sigma)
        self._term = (0.5 * self._s_density * (sympy.erf((self.z - self._s_mu) /
                      (self._s_sigma*sympy.sqrt(2))) + 1))
    