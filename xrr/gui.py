# -*- coding: utf-8 -*-
from gi.repository import Gtk, GObject
from matplotlib.figure import Figure
from matplotlib.backends.backend_gtk3cairo import FigureCanvasGTK3Cairo
#from matplotlib.backends.backend_gtk3agg import FigureCanvasGTK3Agg
from matplotlib.backends.backend_gtk3 import NavigationToolbar2GTK3
import matplotlib.pyplot as plt
import numpy

def fresnel(qc, qz, roughness=0):
  return numpy.exp(-qz**2 * roughness**2 / 2) * \
                  numpy.abs((qz - numpy.sqrt((qz**2 - qc**2) + 0j)) / \
                  (qz+numpy.sqrt((qz**2 - qc**2)+0j)))**2


class ModelPlot(Gtk.VBox):
    
    def __init__(self, model):
        super(ModelPlot, self).__init__()
        self._model = model
        self._prepare_canvas()
        self._prepare_tools()
        self.update()
        
    def _prepare_canvas(self):
        self.fig = Figure()
        self.fig.patch.set_visible(False)
        self.canvas = FigureCanvasGTK3Cairo(self.fig)
        self.pack_start(self.canvas, True, True, 0)
        self.canvas.draw()
        
    def _prepare_tools(self):
        toolbar = NavigationToolbar2GTK3(self.canvas, self)
        self.pack_start(toolbar, False, True, 0)
        self.reorder_child(toolbar, 0)
        
    def update(self):
        pass


class DensityPlot(ModelPlot):
    
    def __init__(self, model, zlims=(-20, 20)):
        self._zmin, self._zmax = zlims
        super(DensityPlot, self).__init__(model)
        
    def update(self):
        z = numpy.linspace(self._zmin, self._zmax, 200)
        rho = self._model(z)
        ax = self.fig.gca()
        ax.clear()
        ax.plot(z, rho)
        self.canvas.draw()
        
        
class ReflectivityPlot(ModelPlot):
    
    def __init__(self, model, data, qzlims=(0, 3)):
        self._qzmin, self._qzmax = qzlims
        self._data = data
        self._func, self._params = model.get_fit_func()
        self._line_model = None
        self._first_run = True
        self._ax = None
        super(ReflectivityPlot, self).__init__(model)
        
    def update(self):
        if self._ax is None:
            self._ax = self.fig.gca()
        ax = self._ax
        self._func, params = self._model.get_fit_func()
        #params = self._model.get_params()
        p = [param.value for param in params]
        
        if self._first_run:
            ax.set_yscale("log")
            ax.plot(self._data[0], self._data[1], ls="none", marker="o",
                    mec="#424242", mew=1.5, mfc="white")      
            self._first_run = False
        if self._line_model is not None:
            line = [line for line in ax.lines if line.get_label()=="model"][0]
            ax.lines.remove(line)
            
        qz = numpy.linspace(self._qzmin, self._qzmax, 100)
        self._line_model = ax.plot(qz, self._func(qz, *p) * fresnel(0.0667, qz),
                                   c="#cc0000", label="model")
        
        self.canvas.draw()
        
        
class ParameterList(Gtk.Table):
    
    __gsignals__ = {"parameter-changed": (GObject.SIGNAL_RUN_FIRST, None,
                               (GObject.TYPE_PYOBJECT,))}
    
    def __init__(self, model):
        super(ParameterList, self).__init__(2, 10)
        self._model = model
        self._prepare_list()
        
    def _prepare_list(self):
        f, params = self._model.get_fit_func()
        n = 0
        for param in params:
            l = Gtk.Label(param.name)
            self.attach(l, 0, 1, n, n+1, yoptions=Gtk.AttachOptions.SHRINK)
            cntrl = Gtk.SpinButton()
            cntrl.set_digits(2)
            adjustment = Gtk.Adjustment(param.value, -1000, 1000, 0.1, 1, 1)
            cntrl.set_adjustment(adjustment)
            cntrl.connect("value-changed", self._cb_param_changed, param)
            self.attach(cntrl, 1, 2, n, n+1, yoptions=Gtk.AttachOptions.SHRINK)
            n += 1
            
    def _cb_param_changed(self, widget, param):
        param.value = widget.get_value()
        self.emit("parameter-changed", param)


class ModelGUIMainWindow(Gtk.Window):
    
    def __init__(self, model, data):
        super(ModelGUIMainWindow, self).__init__()
        self.set_default_size(800, 600)
        
        self._p_density = DensityPlot(model)
        self._p_reflectivity = ReflectivityPlot(model, data)
        self._l_params = ParameterList(model)
        self._l_params.connect("parameter-changed", self._cb_param_changed)
        
        vbox = Gtk.VBox()
        plot_box = Gtk.HBox()
        plot_box.pack_start(self._p_density, True, True, 0)
        plot_box.pack_start(self._p_reflectivity, True, True, 0)
        vbox.pack_start(plot_box, True, True, 0)
        vbox.pack_start(self._l_params, True, False, 0)
        self.add(vbox)
        
    def _cb_param_changed(self, plist, param):
        self._p_density.update()
        self._p_reflectivity.update()
        

class ModelGUI(object):
    
    def __init__(self, model, xdata, ydata):
        super(ModelGUI, self).__init__()
        self._model = model
        self._data = xdata, ydata
        self._window = ModelGUIMainWindow(model, self._data)
        self._window.connect("destroy", self._cb_close_window)
    
    def show(self):
        self._window.show_all()
        Gtk.main()
        
    def _cb_close_window(self, *args):
        Gtk.main_quit()
