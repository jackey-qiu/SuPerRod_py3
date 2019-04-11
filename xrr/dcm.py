import math
import sympy
import scipy

from .layer import StackedLayer, LayerParameter, GaussianLayer,ErrorLayer




class DCMLayer(StackedLayer):
    
    display_layers = 25
    
    def __init__(self, d, si, sb,  r_beach, z_beach, sigma_beach, density_w, mu_w, sigma_w, Z, position=0, formfactor=lambda x: 1.0, **kwargs):
        super(DCMLayer, self).__init__(**kwargs)
        self._s_d = sympy.Symbol("d", real=True, positive=True)
        self._s_position = sympy.Symbol("position", real=True, positive=True)
        self._s_si = sympy.Symbol("si", real=True, positive=True)
        self._s_sb = sympy.Symbol("sb", real=True, positive=True)
        
        self.d = LayerParameter(d, self._s_d)
        self.position = LayerParameter(position, self._s_position)
        self.sigma_i = LayerParameter(si, self._s_si)
        self.sigma_b = LayerParameter(sb, self._s_sb)
        self.Z = Z
        self.form_factor = formfactor
        
        self.rho_ad=sympy.Symbol("r_beach", real=True, positive=True)
        self.sigma_ad=sympy.Symbol("sigma_beach", real=True, positive=True)
        self.z_ad=sympy.Symbol("z_beach", real=True)
        
        self.r_beach = LayerParameter(r_beach, self.rho_ad)
        self.sigma_beach = LayerParameter(sigma_beach, self.sigma_ad)
        self.z_beach = LayerParameter(z_beach, self.z_ad)
        
        
        self._s_density_w= sympy.Symbol("density_w", real=True, positive=True)
        self._s_mu_w = sympy.Symbol("mu_w", real=True)
        self._s_sigma_w = sympy.Symbol("sigma_w", real=True, positive=True)
        self.mu_w = LayerParameter(mu_w, self._s_mu_w)
        self.sigma_w = LayerParameter(sigma_w, self._s_sigma_w)
        self.density_w = LayerParameter(density_w, self._s_density_w)
        

        
        self._make_term()


   
    
    def _make_term(self):
        
        mu_ad= self.z_ad
        sigma_ad= self.sigma_ad
        amp_ad=  self.rho_ad /(self.sigma_ad *  math.sqrt(2 * math.pi) )
        layer= GaussianLayer(amp_ad, mu_ad, sigma_ad)
        
        self.stack.append(layer)
    
        d_w=self._s_density_w
        mu_w=self._s_mu_w
        s_w=self._s_sigma_w 
        layer= ErrorLayer(d_w, mu_w, s_w)  
        self.stack.append(layer)
        
        for n in range(self.display_layers):
            mu = self._s_d * n - self._s_position
            sigma = sympy.sqrt(n * self._s_sb**2 + self._s_si**2)
            amp = (self._s_d * self._s_density /
                  (sigma * math.sqrt(2 * math.pi)))
            layer = GaussianLayer(amp, mu, sigma)
            self.stack.append(layer)
      
                   
    
            
        
    def structure_factor(self, qz, fit=False):
        f = (self.form_factor(qz) + 0.2) / (self.Z + 0.2)
#(self._s_density / (self._s_density - self._s_density_w  ))
        dcm =   1j * qz * self._s_d * self._s_density * (sympy.exp(1j * qz * self._s_position - qz**2 * self._s_si**2 / 2.0) /
                (1 - sympy.exp(1j*qz * self._s_d - qz**2 * self._s_sb**2 / 2.0)))
        beach = 1j* qz * self.rho_ad * sympy.exp(1j * qz * self.z_ad - (qz**2 * self.sigma_ad**2) / 2.0 ) 
        electrolyte= -self._s_density_w * sympy.exp(-1j *qz * self._s_mu_w - qz**2 * self._s_sigma_w**2 /2 ) 
        
        temp=  1/ (self._s_density + self._s_density_w ) * ( f *(dcm+ beach) + electrolyte)
        fac = self._substitute(temp, fit=fit)
        

        return fac
        #return sympy.simplify(fac)


class DCMLayerModified(DCMLayer):
    
    display_layers = 25
    
    def __init__(self, rho_first, z_first, d, si, sb, density_w, mu_w, sigma_w, Z, position=0, formfactor=lambda x: 1.0, **kwargs):
        self._s_density_first = sympy.Symbol("rho_first", real=True, positive=True)
        self._s_z_first = sympy.Symbol("z_first", real=True, positive=True)
        self.density_first = LayerParameter(rho_first, self._s_density_first)
        self.z_first = LayerParameter(z_first, self._s_z_first)
        super(DCMLayerModified, self).__init__(d, si, sb, 0, 0, 0, density_w, mu_w, sigma_w, Z, position=0, formfactor=lambda x: 1.0, **kwargs)
        
    def _make_term(self):
        
        d_w=self._s_density_w
        mu_w=self._s_mu_w
        s_w=self._s_sigma_w 
        layer= ErrorLayer(d_w, mu_w, s_w)  
        self.stack.append(layer)
        
        # first layer modification
        amp = self._s_d * self._s_density_first / (self._s_si * math.sqrt(2 * math.pi))
        mu = self._s_z_first
        sigma = self._s_si
        layer = GaussianLayer(amp, mu, sigma)
        self.stack.append(layer)
        # normal DCM layers
        for n in range(1, self.display_layers):
            mu = self._s_d * n - self._s_position
            sigma = sympy.sqrt(n * self._s_sb**2 + self._s_si**2)
            amp = (self._s_d * self._s_density /
                  (sigma * math.sqrt(2 * math.pi)))
            layer = GaussianLayer(amp, mu, sigma)
            self.stack.append(layer)
        
    def structure_factor(self, qz, fit=False):
        f = (self.form_factor(qz) + 0.2) / (self.Z + 0.2) #(1 / ( 0.015)) *
        
        
        
        # modification by Andrea 03/2019
      
        dcm= (1j * qz* self._s_d * self._s_density * sympy.exp(- (qz**2 * self._s_si**2 )/2.0 ) * ( sympy.exp( 1j * qz * self._s_d - (qz**2 * self._s_sb**2 )/2.0 ) / (1 - sympy.exp( 1j * qz * self._s_d - (qz**2 * self._s_sb**2 )/2.0 ))  ))
        mod= ( 1j * qz* self._s_d *  self._s_density_first  *  sympy.exp(( 1j * qz* self._s_z_first - (qz**2 * self._s_si**2 )/2.0 )  ))
        #mod= ( 1j * qz* self._s_d * ( self._s_density_first -self._s_density )  *  sympy.exp(( 1j * qz* self._s_z_first - (qz**2 * self._s_si**2 )/2.0 )  ))
        electrolyte = -self._s_density_w * sympy.exp(-1j *qz * self._s_mu_w - qz**2 * self._s_sigma_w**2 /2 )
        temp = 1/ (0.015 + self._s_density_w ) * ( f* (mod + dcm) + electrolyte) #65*
        
        fac = self._substitute(temp, fit=fit)
        return sympy.simplify(fac)
