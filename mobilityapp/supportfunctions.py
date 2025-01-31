import numpy as np
from lmfit import Model, Parameters, minimize

from tkinter import *
import tkinter as tk

class CreateToolTip(object):
    """
    create a tooltip for a given widget
    """
    def __init__(self, widget, text='widget info'):
        self.waittime = 500     #miliseconds
        self.wraplength = 180   #pixels
        self.widget = widget
        self.text = text
        self.widget.bind("<Enter>", self.enter)
        self.widget.bind("<Leave>", self.leave)
        self.widget.bind("<ButtonPress>", self.leave)
        self.id = None
        self.tw = None

    def enter(self, event=None):
        self.schedule()

    def leave(self, event=None):
        self.unschedule()
        self.hidetip()

    def schedule(self):
        self.unschedule()
        self.id = self.widget.after(self.waittime, self.showtip)

    def unschedule(self):
        id = self.id
        self.id = None
        if id:
            self.widget.after_cancel(id)

    def showtip(self, event=None):
        x = y = 0
        x, y, cx, cy = self.widget.bbox("insert")
        x += self.widget.winfo_rootx() + 25
        y += self.widget.winfo_rooty() + 20
        # creates a toplevel window
        self.tw = tk.Toplevel(self.widget)
        # Leaves only the label and removes the app window
        self.tw.wm_overrideredirect(True)
        self.tw.wm_geometry("+%d+%d" % (x, y))
        label = tk.Label(self.tw, text=self.text, justify='left',
                       background="#ffffff", relief='solid', borderwidth=1,
                       wraplength = self.wraplength)
        label.pack(ipadx=1)

    def hidetip(self):
        tw = self.tw
        self.tw= None
        if tw:
            tw.destroy()
            
def perform_deriv_fit(Vg,G,dGdVg,Gsmooth,smoothing,Vmin,Vmax):
    
    if Vmin != 'Min':
        Vmin=float(Vmin)
        min_ind=(np.abs(Vg - Vmin)).argmin()
    else:
        min_ind=1 #To avoid ends
        Vmin=Vg[min_ind]
    if Vmax != 'Max':
        Vmax=float(Vmax)
        max_ind=(np.abs(Vg - Vmax)).argmin()
    else:
        max_ind=int(Vg.shape[0]-2)
        Vmax=Vg[max_ind]

    initial_a = -1
    initial_x0 = (Vmax-Vmin)/2
    initial_u0 = (Vmax-Vmin)/10
        
    def asym_lorentzian(A, x0, a, u0, x):
        u = 2 * u0 / (1 + np.exp(a * (x - x0)))
        return 2 * A / (np.pi * u) * (1 + 4 * ((x - x0) / u)**2)**-1

    # Define the lmfit model for the asymmetric Lorentzian
    def model(params, x):
        A = params['A'].value
        x0 = params['x0'].value
        a = params['a'].value
        u0 = params['u0'].value
        return asym_lorentzian(A, x0, a, u0, x)

    # Create a lmfit Parameters object with initial guesses for the asymmetric lorentzian
    params = Parameters()
    params.add('A', value=np.max(dGdVg),min=0)  # Use max of data as the initial guess for A
    params.add('a', value=initial_a,max=0)
    params.add('x0', value=initial_x0)
    params.add('u0', value=initial_u0,min=0)

        # Define the residual function to fit the gradient
    def residual(params, x, data):
        model_values = model(params, x)
        return data - model_values

    # Perform the fit using lmfit
    result = minimize(residual, params, args=(Vg[min_ind:max_ind], dGdVg[min_ind:max_ind]))

    # Generate two versions of the best-fit curve; one for plotting, one for extracting values more precisely
    Vgfine=np.linspace(Vg[min_ind],Vg[max_ind],10001) #Makes the method more robust to sparse datasets.
    deriv_fit = model(result.params, Vg)
    deriv_fitfine = model(result.params, Vgfine)

    infl_slope=deriv_fitfine.max()    #Slope in G(Vg) at the inflection point is simply the maximum of the derivative
    infl_ind=deriv_fitfine.argmax()   #Array index where inflection point occurs
    Vg_infl=Vgfine[infl_ind]             #Value of Vg at the inflection point

    if smoothing !=0:
        G_infl=np.interp(Vg_infl, Vg, Gsmooth)#Gsmooth[infl_ind]
    else:
        G_infl=np.interp(Vg_infl, Vg, Gsmooth)#G[infl_ind]  #Value of G at the inflection point

    G_intercept = G_infl-infl_slope*Vg_infl  #Finding 'threshold' voltage
    Vth=-G_intercept/infl_slope
    
    V_Rs = Vg_infl+2*(Vg_infl-Vth)   #Vg above which we will use to calculate series resistance
    
    thresholdline=infl_slope*Vg+G_intercept      #Draw a line tangential with the inflection point. Vg-intercept is Vth

    V0 = Vth-2*(Vg_infl-Vth)         #Vg for which density extrapolates to zero.

    return V0,Vth,Vg_infl,V_Rs,thresholdline,deriv_fit

def drude(x, Rs,mu,Vth,L,c):    # drude fit 
    return 1/(Rs + L**2/(mu*c*(x-Vth)))
model_drude = Model(drude)

def perform_Rs_fit(Vg,G,V0,V_Rs,initial_Rs,initial_mu,L,c):
    params_Rs = Parameters()
    params_Rs.add('Rs',value=initial_Rs)
    params_Rs.add('mu',value=initial_mu)
    params_Rs.add('Vth',value=V0,vary=False) #Do not vary this parameter, since we know it now.
    params_Rs.add('L',value=L,vary=False)
    params_Rs.add('c',value=c,vary=False)

    # Perform the Drude fit over limited range
    V_Rs_ind = (np.abs(Vg - V_Rs)).argmin()
    result_drudeRs = model_drude.fit(G[V_Rs_ind:], params_Rs, x=Vg[V_Rs_ind:])
    Rs = result_drudeRs.params['Rs'].value             ## This is THE value of Rs
    mu_Rs = result_drudeRs.params['mu'].value          # Meaningless value

    Rs_fit = drude(Vg[V_Rs_ind:],Rs,mu_Rs,V0,L,c)
    
    return Rs,Rs_fit,result_drudeRs,V_Rs_ind

def perform_drude_fit(Vg,G,Vth,initial_Rs,initial_mu,L,c):
    params_drude = Parameters()
    params_drude.add('Rs',value=initial_Rs) #Allow to var
    params_drude.add('mu',value=initial_mu)
    params_drude.add('Vth',value=Vth) #Can/should vary now.
    params_drude.add('L',value=L,vary=False)
    params_drude.add('c',value=c,vary=False)

    # Perform the Drude fit over 'full' range. Which is still Vg>Vth
    Vth_ind=(np.abs(Vg - Vth)).argmin()
    result_drude = model_drude.fit(G[Vth_ind:], params_drude, x=Vg[Vth_ind:])
    Rs_drude = result_drude.params['Rs'].value
    mu_drude = result_drude.params['mu'].value
    Vth_drude=result_drude.params['Vth'].value
    
    drude_fit = drude(Vg[Vth_ind:],Rs_drude,mu_drude,Vth_drude,L,c)
    
    return mu_drude,drude_fit,Rs_drude,result_drude,Vth_ind