import numpy as np
from lmfit import Model, Parameters, minimize
from scipy.signal import savgol_filter
import matplotlib.pyplot as plt

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
            
def perform_deriv_fit(Vg,G,dGdVg,Gsmooth,smoothing,Vmin,Vmax,holes=False):
    
    if holes:           #Flipping the data as if it was electrons means we can stick with good initial guesses no matter what
        dGdVg=-dGdVg[::-1]

    minoptions=['Min','','min']
    if Vmin not in minoptions:
        Vmin=float(Vmin)
        min_ind=(np.abs(Vg - Vmin)).argmin()
    else:
        min_ind=1 #To avoid ends
        Vmin=Vg[min_ind]
    maxoptions=['Max','','max']
    if Vmax not in maxoptions:
        Vmax=float(Vmax)
        max_ind=(np.abs(Vg - Vmax)).argmin()
    else:
        max_ind=int(Vg.shape[0]-2)
        Vmax=Vg[max_ind]

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
    params.add('x0', value=(Vmax+Vmin)/2)  # Use the middle of the data as the initial guess for x0
    params.add('a', -1,max=0)
    params.add('u0', value=(Vmax-Vmin)/10,min=0)
    
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
    
    if holes: #Since we flipped the derivative data earlier, we need to reverse everything.
        infl_slope=-infl_slope
        deriv_fit=-deriv_fit[::-1]
        Vg_infl=Vgfine[-infl_ind]             #Value of Vg at the inflection point
    else:
        Vg_infl=Vgfine[infl_ind]             #Value of Vg at the inflection point

    if smoothing !=0:
        G_infl=np.interp(Vg_infl, Vg, Gsmooth)#Gsmooth[infl_ind]
    else:
        G_infl=np.interp(Vg_infl, Vg, G)#G[infl_ind]  #Value of G at the inflection point

    G_intercept = G_infl-infl_slope*Vg_infl  #Finding 'threshold' voltage
    Vth=-G_intercept/infl_slope
    
    V_Rs = Vg_infl+2*(Vg_infl-Vth)   #Vg above which we will use to calculate series resistance
    
    inflectionline=infl_slope*Vg+G_intercept      #Draw a line tangential with the inflection point. Vg-intercept is Vth

    V0 = Vth-2*(Vg_infl-Vth)         #Vg for which density extrapolates to zero.

    return V0,Vth,Vg_infl,V_Rs,inflectionline,deriv_fit,result

def manual_inflection(Vg,G,Gsmooth,smoothing,Vg_infl,infl_slope):
    if smoothing !=0:
        G_infl=np.interp(Vg_infl, Vg, Gsmooth)#Gsmooth[infl_ind]
    else:
        G_infl=np.interp(Vg_infl, Vg, G)#G[infl_ind]  #Value of G at the inflection point

    G_intercept = G_infl-infl_slope*Vg_infl  #Finding 'threshold' voltage
    Vth=-G_intercept/infl_slope
    
    V_Rs = Vg_infl+2*(Vg_infl-Vth)   #Vg above which we will use to calculate series resistance
    
    inflectionline=infl_slope*Vg+G_intercept      #Draw a line tangential with the inflection point. Vg-intercept is Vth

    V0 = Vth-2*(Vg_infl-Vth)         #Vg for which density extrapolates to zero.

    return V0,Vth,V_Rs,inflectionline

def drude(x, Rs,mu,Vth,L,c,holes=False):    # drude fit 
    if holes:
        return 1/(Rs + L**2/(mu*c*(Vth-x)))
    else:
        return 1/(Rs + L**2/(mu*c*(x-Vth)))
model_drude = Model(drude)

def perform_Rs_fit(Vg,G,V0,V_Rs,initial_Rs,initial_mu,L,c,holes=False):
    params_Rs = Parameters()
    params_Rs.add('Rs',value=initial_Rs)
    params_Rs.add('mu',value=initial_mu)
    params_Rs.add('Vth',value=V0,vary=False) #Do not vary this parameter, since we know it now.
    params_Rs.add('L',value=L,vary=False)
    params_Rs.add('c',value=c,vary=False)
    params_Rs.add('holes',value=holes,vary=False)

    # Perform the Drude fit over limited range, where Rchannel starts to fall below Rs
    V_Rs_ind = (np.abs(Vg - V_Rs)).argmin()
    if holes:
        result_drudeRs = model_drude.fit(G[:V_Rs_ind], params_Rs, x=Vg[:V_Rs_ind])
    else:
        result_drudeRs = model_drude.fit(G[V_Rs_ind:], params_Rs, x=Vg[V_Rs_ind:])
    Rs = result_drudeRs.params['Rs'].value             ## This is THE value of Rs
    mu_Rs = result_drudeRs.params['mu'].value          # Meaningless value

    if holes:
        Rs_fit = drude(Vg[:V_Rs_ind],Rs,mu_Rs,V0,L,c,holes)
    else:
        Rs_fit = drude(Vg[V_Rs_ind:],Rs,mu_Rs,V0,L,c)
    
    return Rs,Rs_fit,V_Rs_ind,result_drudeRs

def perform_drude_fit(Vg,G,Vth,initial_Rs,initial_mu,L,c,holes=False):
    params_drude = Parameters()
    params_drude.add('Rs',value=initial_Rs) #Allow to vary for the purpose of illustration.
    params_drude.add('mu',value=initial_mu)
    params_drude.add('Vth',value=Vth) #Can/should vary now.
    params_drude.add('L',value=L,vary=False)
    params_drude.add('c',value=c,vary=False)
    params_drude.add('holes',value=holes,vary=False)

    # Perform the Drude fit over 'full' range. Which is still Vg>Vth
    Vth_ind=(np.abs(Vg - Vth)).argmin()
    if holes:
        result_drude = model_drude.fit(G[:Vth_ind], params_drude, x=Vg[:Vth_ind])
    else:
        result_drude = model_drude.fit(G[Vth_ind:], params_drude, x=Vg[Vth_ind:])
    Rs_drude = result_drude.params['Rs'].value
    mu_drude = result_drude.params['mu'].value
    Vth_drude=result_drude.params['Vth'].value
    
    if holes:
        drude_fit = drude(Vg[:Vth_ind],Rs_drude,mu_drude,Vth_drude,L,c,holes)
    else:
        drude_fit = drude(Vg[Vth_ind:],Rs_drude,mu_drude,Vth_drude,L,c)
    
    return mu_drude,drude_fit,Rs_drude,Vth_ind,result_drude

def perform_entire_prodecure(Vg,G,smoothing,Vmin,Vmax,L,C,CperA,initial_Rs,initial_mu,holes=False,plotting=True):

    datadict={}
    paramdict={}
    
    if Vg[-1]<Vg[0]:
        Vg=Vg[::-1]
        G=G[::-1]
    datadict['Vg (V)']=Vg
    datadict['G (S)']=G

    if smoothing!=0:
        Gsmooth=savgol_filter(G, int(G.shape[0]*smoothing), 5, deriv=0, delta=1.0, axis=-1, mode='interp', cval=0.0)
        dGdVg = np.gradient(Gsmooth, Vg)
    else:
        Gmooth=0 #Necessary later because I'm not a good developer, sorry
        dGdVg = np.gradient(G, Vg)
    datadict['dGdVg (S/V)']=dGdVg

    V0,Vth,Vg_infl,V_Rs,inflectionline,deriv_fit,result_deriv_fit=perform_deriv_fit(Vg,G,dGdVg,Gsmooth,smoothing,Vmin=Vmin,Vmax=Vmax,holes=holes)
    datadict['dGdVg fit (S/V)']=deriv_fit
    datadict['Inflection fit (S)']=inflectionline
    paramdict['V0 (V)']=V0
    paramdict['Vth (V)']=Vth
    paramdict['Vg_infl (V)']=Vg_infl
    paramdict['V_Rs (V)']=V_Rs
        
    Rs,Rs_fit,V_Rs_ind,result_drudeRs=perform_Rs_fit(Vg,G,V0,V_Rs,initial_Rs,initial_mu,L,C,holes)
    paramdict['Rs (Ohm)']=Rs
    if holes:
        datadict['Vg for Rs fit (V)']=Vg[:V_Rs_ind]
    else:
        datadict['Vg for Rs fit (V)']=Vg[V_Rs_ind:]
    datadict['Rs fit (S)']=Rs_fit
    
    if holes:
        density=CperA*(V0-Vg)/1.602176634e-19
        mu_eff=L**2/(C*(V0-Vg)*((1/G)-Rs))
    else:
        density=CperA*(Vg-V0)/1.602176634e-19
        mu_eff=L**2/(C*(Vg-V0)*((1/G)-Rs))
    datadict['density (1/m2)']=density
    datadict['mu_eff (m2/Vs)']=mu_eff
    
    mu_drude,drude_fit,Rs_drude,Vth_ind,result_drude=perform_drude_fit(Vg,G,Vth,initial_Rs,initial_mu,L,C,holes)
    paramdict['mu_FET (m2/Vs)']=mu_drude
    #paramdict['Rs_drude (Ohm)']=Rs_drude
    if holes:
        datadict['Vg for mu_FET fit (V)']=Vg[:Vth_ind]
    else:
        datadict['Vg for mu_FET fit (V)']=Vg[Vth_ind:]
    datadict['mu_FET fit (S)']=drude_fit
    
    if plotting==True:
        plt.plot(Vg,G,label='data')
        plt.plot(Vg,inflectionline,label='inflection')
        if holes:
            plt.plot(Vg[:V_Rs_ind],Rs_fit,label='Rs fit')
            plt.plot(Vg[:Vth_ind],drude_fit,label='drude fit')
        else:
            plt.plot(Vg[V_Rs_ind:],Rs_fit,label='Rs fit')
            plt.plot(Vg[Vth_ind:],drude_fit,label='drude fit')
        plt.xlabel('Vg (V)')
        plt.ylabel('Conductance (S)')
        plt.ylim([G.min()-(G.max()-G.min())/10,G.max()+(G.max()-G.min())/10])
        plt.legend()
        plt.show()

        plt.plot(Vg,dGdVg)
        plt.plot(Vg,deriv_fit)
        plt.xlabel('Vg (V)')
        plt.ylabel('$dG/dVg$ (S/V)')
        plt.show()
        
        mu_drude_array=np.full(Vg.shape[0],mu_drude)
        plotstart=(np.abs(Vg - (2*Vth-Vg_infl))).argmin()
        if holes:
            plt.plot(density[:plotstart]*1e-12/1e4,mu_eff[:plotstart]*1e4,label='mu_eff')
            plt.plot(density[:plotstart]*1e-12/1e4,mu_drude_array[:plotstart]*1e4,label='mu_drude')
        else:
            plt.plot(density[plotstart:]*1e-12/1e4,mu_eff[plotstart:]*1e4,label='mu_eff')
            plt.plot(density[plotstart:]*1e-12/1e4,mu_drude_array[plotstart:]*1e4,label='mu_drude')
        plt.xlabel('Carrier density x 10$^{12}$ (cm$^{-2}$)')
        plt.ylabel('Mobility (cm$^2$/(Vs))')
        plt.legend()
        plt.show()
    return datadict,paramdict