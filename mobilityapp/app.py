# Import necessary libraries
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg, NavigationToolbar2Tk)
from matplotlib.figure import Figure

from scipy.signal import savgol_filter
import numpy as np
import json
import csv

from tkinter import *
from tkinter import ttk
import tkinter as tk

from .supportfunctions import *

def run():
    root = Tk()
    root.title('Calculate mobility from two terminal data')
    ## Pick the datafile
    dataframe = ttk.Frame(root, padding='3 3 12 12')
    dataframe.grid(column=0,row=0,sticky='N')
    data={}       #Dictionaries to store data and parameters both loaded and created since tkinter cannot read function returns
    paramdict={}
    exportdatadict={}
    loadedfile=StringVar()
    loadedfile.set('No file loaded')
    #exportdata = [["Vg (V)", "dGdVg (S/V)","Inflection fit (S)","R_s fit (S)","density (1/m2)","mu_eff (m2/Vs)","mu_FET_fit (S)"]]
    exportdatavar = StringVar()
    exportparamsvar=StringVar()
    def set_exportdata():
        exportdataarray=[]
        for param in exportdatadict:
            exportdataarray.append(f'{param}   {exportdatadict[param]}')
        exportdatavar.set(exportdataarray)
    def set_exportparams():
        exportparams=[]
        for param in paramdict:
            exportparams.append(f'{param}   {paramdict[param]}')
        exportparamsvar.set(exportparams)
    def load_data():
        data.clear()  #Clear dictionaries when loading a new datafile to make sure no old data/params get stuck around
        paramdict.clear()
        exportdatadict.clear()
        exportdatavar.set('')
        exportparamsvar.set('')
        clear_plot()
        clear_deriv()
        clear_mobility()
        
        filename = tk.filedialog.askopenfilename()
        with open(filename) as f:
            loadeddata=np.loadtxt(f)
        loadedfile.set(filename)
        if loadeddata.shape[0]>loadeddata.shape[1]: #Assume that the data is going to have more measurement points than columns
            loadeddata=np.transpose(loadeddata)
        numcols=loadeddata.shape[0]
        data['numcols']=numcols
        for i in range(numcols):
            data[str(i)]=loadeddata[i]
        menu=Vgcol_dropmenu['menu']
        menu.delete(0, 'end')
        for i in range(numcols):
            menu.add_command(label=str(i),command=lambda value=i: Vgcolumn.set(value))
        menu=Gcol_dropmenu['menu']
        menu.delete(0, 'end')
        for i in range(numcols):
            menu.add_command(label=str(i),command=lambda value=i: Gcolumn.set(value))
            
            
    loadframe = ttk.Frame(dataframe, padding='3 3 10 10')
    loadframe.grid(column=0,row=0,sticky=('N,W,E,S'))
    loadlabel=Label(loadframe,text='Select a data file (.dat, .txt or .csv)\n Should contain Vg and G (or I) in separate columns.')
    loadlabel.grid(row=0)
    loadbutton = Button(loadframe, text='1) Select data file', fg='green', command=load_data)
    loadbutton.grid(row=1)
    loadtt=CreateToolTip(loadbutton,'The data is loaded by numpy.loadtxt. It assumes at least two columns, one for '
                        'Vg and one for G (or I). More columns are OK, simply choose correct columns below. '
                        'If you have a more esoteric datatype, e.g. database, please reshape it and save as csv-like.')
    loadedfileentry=Entry(loadframe,textvariable=loadedfile,width=45)
    loadedfileentry.config(state='disabled')
    loadedfileentry.grid(row=2)
    colframe = ttk.Frame(dataframe, padding='3 3 10 10')
    colframe.grid(column=0,row=1)
    Label(colframe,text='Specify the data columns').grid(row=0,columnspan=4)
    Label(colframe,text='Vg column').grid(row=1,column=0)
    Label(colframe,text='G (or I) column').grid(row=1,column=2)
    Gcolumn=StringVar()
    Gcolumn.set('1')
    Vgcolumn=StringVar()
    Vgcolumn.set('0')
    # Vgcol_entry=Entry(colframe,textvariable=Vgcolumn,width=4).grid(row=1,column=1,padx='0 10')
    # Gcol_entry=Entry(colframe,textvariable=Gcolumn,width=4).grid(row=1,column=3)
    Vgcol_dropmenu=OptionMenu(colframe,Vgcolumn,*[str(i) for i in range(2)])
    Vgcol_dropmenu.grid(row=1,column=1)
    Gcol_dropmenu=OptionMenu(colframe,Gcolumn,*[str(i) for i in range(2)])
    Gcol_dropmenu.grid(row=1,column=3)
    unitsframe= ttk.Frame(dataframe, padding='3 3 10 10')
    unitsframe.grid(column=0,row=2)
    Label(unitsframe,text='Select units for conductance data').grid(row=0)
    options=['S','2e2/h','e2/h']
    Gunits=StringVar()
    Gunits.set('S')
    Gdropmenu=OptionMenu(unitsframe,Gunits,*options)
    Gdropmenu.grid(row=1)
    Gdroptt=CreateToolTip(Gdropmenu,'Very important to make sure this is set correctly')
    GorIframe = ttk.Frame(dataframe, padding='3 3 10 10')
    GorIframe.grid(column=0,row=3)
    Label(GorIframe,text='Or, specify Vsd if supplying current').grid(row=0,columnspan=4)
    options=['G provided','I provided']
    GorI=StringVar()
    GorI.set('G provided')
    GorIdropmenu=OptionMenu(GorIframe,GorI,*options)
    GorIdropmenu.grid(row=1,column=0,sticky='E')
    Label(GorIframe,text='Vsd').grid(row=1,column=1,sticky='E')
    Vsd=StringVar()
    Vsd.set('0')
    Vsd_entry=Entry(GorIframe,textvariable=Vsd,width=6).grid(row=1,column=2)
    Label(GorIframe, text='(V)').grid(row=1,column=3,sticky='W')

        
    ## Enter parameters
    paramsframe = ttk.Frame(root, padding='3 3 12 3')
    paramsframe.grid(column=0,row=1,sticky=('N,W,E,S'))
    Label(paramsframe,text='4) Enter geometrical properties of your device',fg='green').grid(row=0,columnspan=3)
    Label(paramsframe,text='Enter intial guesses for series resistance and mobility').grid(row=4,columnspan=3)
    c=StringVar()
    c.set('5.3e-15')
    L=StringVar()
    L.set('3.6e-6')
    W=StringVar()
    W.set('280e-9')
    initial_Rs=StringVar()
    initial_Rs.set('10000')
    initial_mu=StringVar()
    initial_mu.set('4000')
    def update_params(*args):
        if Ctype.get()=='Width':
            Wunits.config(text='(F)')
        else:
            Wunits.config(text='(F/m^2)')
    Ctype=StringVar()
    Ctype.set('Width')
    Ctype.trace('w',update_params)
    Label(paramsframe, text='Capacitance').grid(row=1,sticky='E')
    Label(paramsframe, text='Length').grid(row=2,sticky='E')
    Cdropmenu=OptionMenu(paramsframe,Ctype,*['Width','Cap per area'])
    Cdropmenu.grid(row=3,sticky='E')
    #Label(paramsframe, text='Width').grid(row=3,sticky='E')
    Label(paramsframe, text='Initial R_s').grid(row=5,sticky='E')
    Label(paramsframe, text='Initial mu').grid(row=6,sticky='E')
    Label(paramsframe, text='(F)').grid(row=1,column=2,sticky='W')
    Label(paramsframe, text='(m)').grid(row=2,column=2,sticky='W')
    Wunits=Label(paramsframe, text='(m)')
    Wunits.grid(row=3,column=2,sticky='W')
    Label(paramsframe, text='(Ohm)').grid(row=5,column=2,sticky='W')
    Label(paramsframe, text='(cm^2/Vs)').grid(row=6,column=2,sticky='W')
    c_entry=Entry(paramsframe,textvariable=c)
    c_entry.grid(row=1,column=1)
    L_entry=Entry(paramsframe,textvariable=L).grid(row=2,column=1)
    W_entry=Entry(paramsframe,textvariable=W)
    W_entry.grid(row=3,column=1)
    initRs_entry=Entry(paramsframe,textvariable=initial_Rs)
    initRs_entry.grid(row=5,column=1)
    initmu_entry=Entry(paramsframe,textvariable=initial_mu)
    initmu_entry.grid(row=6,column=1)
    ctt=CreateToolTip(c_entry,'It is important to have an accurate calculation/simulation of the '
                    'capacitance to obtain the proper absolute value of mobility. However, '
                    'a wrong value will not change the trend of mobility vs density, just the absolute value.')
    CreateToolTip(initRs_entry,'Setting initial guesses sensibly can help fitting the series resistance properly.')
    CreateToolTip(initmu_entry,'Setting initial guesses sensibly can help fitting the series resistance properly.')
    CreateToolTip(Cdropmenu,'Calculating the density requires the capacitance per unit area, '
                            'which for planar devices is simply C/(length*width). '
                            'However, nanowires and other irregularly shaped devices do not '
                            'have a well-defined width. In this case, we instead need to know '
                            'the capacitance per area as well as the total capcitance. ')
    for child in paramsframe.winfo_children():
        child.grid_configure(padx=2,pady=2)
        
        
    def VgandG(convertunits=False):
        Vg=data[Vgcolumn.get()]
        G=data[Gcolumn.get()]

        if Vg[-1]<Vg[0]: #make sure Vg ascending; easier to work with
            Vg=Vg[::-1]
            G=G[::-1]

        if np.average(G[0:6])>np.average(G[-6:-1]):
            holes=True
        else:
            holes=False

        if GorI.get()=='I provided':
            convertunits=False
            G=G/float(Vsd.get())
            
        if convertunits==True:
            if Gunits.get()=='2e2/h':
                G=G*7.748091729e-5
            elif Gunits.get()=='e2/h':
                G=G*7.748091729e-5/2
                
        return Vg,G,holes
        
        
    ## Window to plot the loaded data
    dataframe = ttk.Frame(root, padding='3 3 12 3')
    dataframe.grid(column=1,row=0,sticky='N')
    datatopframe = ttk.Frame(dataframe)
    datatopframe.grid(row=0)
    databottomframe=ttk.Frame(dataframe)
    databottomframe.grid(row=1)
    ax={}
    figuresize=(3.75,2.25)
    plt.rcParams.update({'font.size': 8})
    fig1 = Figure(figsize=figuresize, dpi=100)
    canvas1 = FigureCanvasTkAgg(fig1, master=databottomframe)  # A tk.DrawingArea.
    canvas1.draw()
    canvas1.get_tk_widget().pack()#(side=TOP, fill=BOTH, expand=1)
    toolbar = NavigationToolbar2Tk(canvas1, databottomframe)
    toolbar.update()
    canvas1.get_tk_widget().pack()#(side=TOP, fill=BOTH, expand=1)
    def plot_data():
        fig1.clf()
        Vg,G,holes=VgandG(convertunits=False)
        
        exportdatadict['Vg (V)']=Vg
        exportdatadict['G (S)']=G
        set_exportdata()
            
        ax[0]=fig1.add_subplot()
        ax[0].plot(Vg,G,'k',label='data')
        ax[0].set_xlabel('Gate voltage (V)')
        if GorI.get()=='I provided':
            ax[0].set_ylabel('Conductance (S)')
        else:
            ax[0].set_ylabel(f'Conductance ({Gunits.get()})')
        ax[0].set_ylim([G.min()-(G.max()-G.min())/10,G.max()+(G.max()-G.min())/10])
        
        fig1.tight_layout()
        ax[0].legend()
        fig1.canvas.draw()
        
    def clear_plot():
        fig1.clf()
        fig1.canvas.draw()
        
    plotdatabutton = Button(datatopframe, text='2) Plot data/Refresh', command=plot_data,fg='green').grid(row=0)
    clearplotbutton = Button(datatopframe, text='Clear plot', command=clear_plot).grid(row=0,column=2)
    
    
    
    ## Window to plot the derivlorentzian fit
    derivframe = ttk.Frame(root, padding='3 3 12 3')
    derivframe.grid(column=2,row=0)
    derivtopframe = ttk.Frame(derivframe)
    derivtopframe.grid(row=0)
    derivbottomframe=ttk.Frame(derivframe)
    derivbottomframe.grid(row=1)
    derivextraframe=ttk.Frame(derivframe)
    derivextraframe.grid(row=2)

    fig2 = Figure(figsize=tuple((4.5,figuresize[1])), dpi=100)
    canvas2 = FigureCanvasTkAgg(fig2, master=derivbottomframe)  # A tk.DrawingArea.
    canvas2.draw()
    canvas2.get_tk_widget().pack()#(side=TOP, fill=BOTH, expand=1)
    toolbar = NavigationToolbar2Tk(canvas2, derivbottomframe)
    toolbar.update()
    canvas2.get_tk_widget().pack()#(side=TOP, fill=BOTH, expand=1)
    def plot_deriv():
        #First do the analysis
        Vg,G,holes=VgandG(convertunits=True)
        
        #Range over which to smooth, as a percentage of the data range
        smoothing=float(smoothingbut.get())/100
        
        if smoothing!=0:
            Gsmooth=savgol_filter(G, int(G.shape[0]*smoothing), 5, deriv=0, delta=1.0, axis=-1, mode='interp', cval=0.0)
            dGdVg = np.gradient(Gsmooth, Vg)
        else:
            Gsmooth=0
            dGdVg = np.gradient(G, Vg)

        exportdatadict['dGdVg (S/V)']=dGdVg
        set_exportdata()
            
        fig2.clf()
        ax[1]=fig2.add_subplot()
        # if holes:
        #     ax[1].plot(Vg,dGdVg*1e3,'k',label='data')
        # else:
        ax[1].plot(Vg,dGdVg*1e3,'k',label='data')
        ax[1].set_xlabel('Gate voltage (V)')
        ax[1].set_ylabel(f'$dG/dV_g$ (mS/V)')
        ax[1].legend()
        fig2.tight_layout()
        fig2.canvas.draw()
            
        # try:
        V0,Vth,Vg_infl,V_Rs,inflectionline,deriv_fit,result_deriv_fit=perform_deriv_fit(Vg,G,dGdVg,Gsmooth,smoothing,Vmin=Vmin.get(),Vmax=Vmax.get(),holes=holes)
        paramdict['V0 (V)']=V0
        paramdict['Vth (V)']=Vth
        paramdict['Vg_infl (V)']=Vg_infl
        paramdict['V_Rs (V)']=V_Rs
                                
        set_exportparams()
        
        exportdatadict['dGdVg fit (S/V)']=deriv_fit
        exportdatadict['Inflection fit (S)']=inflectionline
        set_exportdata()
        # if holes:
        #     ax[1].plot(Vg,-deriv_fit[::-1]*1e3,label='fit')
        # else:
        ax[1].plot(Vg,deriv_fit*1e3,label='fit')
        ax[1].legend()
        fig2.canvas.draw()
        if GorI.get()=='G provided':
            if Gunits.get()=='2e2/h':
                inflectionline=inflectionline/7.748091729e-5
            elif Gunits.get()=='e2/h':
                inflectionline=inflectionline/7.748091729e-5*2
        ax[0].plot(Vg,inflectionline,label='slope at inflection')
        ax[0].legend()
        fig1.canvas.draw()
        # except Exception as e:
        #     print(e)
        #     fig2.clf()
        #     canvas = FigureCanvasTkAgg(fig2, master=derivframe)  # A tk.DrawingArea.
        #     canvas.draw()
        #     canvas.get_tk_widget().grid(row=1,columnspan=11)
        
    def clear_deriv():
        fig2.clf()
        fig2.canvas.draw()
        
    plotderivbutton = Button(derivtopframe, text='3) Fit dGdVg/Refresh', command=plot_deriv,fg='green')
    plotderivbutton.grid(row=0)
    clearderivbutton = Button(derivtopframe, text='Clear plot', command=clear_deriv)
    clearderivbutton.grid(row=0,column=10)
    CreateToolTip(plotderivbutton,'Fit on the derivative of G(Vg). '
                            'The aim is to find the maximum of dG/dVg, which is the inflection point in G(Vg). '
                            'The fit can easily fail, but do not give up! Make sure units are correctly selected. '
                            'Try different amount of smoothing,  '
                            'and limiting the Vg range which the fit is performed over. '
                            'If you have to give up, you can enter the peak coordinates manually below.')
    smoothingbut=StringVar()
    smoothingbut.set('0')
    Label(derivtopframe, text='Smoothing').grid(row=0,column=1,sticky='E')
    smoothing_entry=Entry(derivtopframe,textvariable=smoothingbut,width=3).grid(row=0,column=2)
    Label(derivtopframe, text='%').grid(row=0,column=3,sticky='W')
    Vmin=StringVar()
    Vmin.set('Min')
    Label(derivtopframe, text='Vmin').grid(row=0,column=4,sticky='E')
    Vmin_entry=Entry(derivtopframe,textvariable=Vmin,width=6).grid(row=0,column=5)
    Label(derivtopframe, text='V').grid(row=0,column=6,sticky='W')
    Vmax=StringVar()
    Vmax.set('Max')
    Label(derivtopframe, text='Vmax').grid(row=0,column=7,sticky='E')
    Vmax_entry=Entry(derivtopframe,textvariable=Vmax,width=6).grid(row=0,column=8)
    Label(derivtopframe, text='V').grid(row=0,column=9,sticky='W')
    
    def plot_manual_inflection():
        Vg,G,holes=VgandG(convertunits=True)
        
        #Range over which to smooth, as a percentage of the data range
        smoothing=float(smoothingbut.get())/100
        
        if smoothing!=0:
            Gsmooth=savgol_filter(G, int(G.shape[0]*smoothing), 5, deriv=0, delta=1.0, axis=-1, mode='interp', cval=0.0)
        else:
            Gsmooth=0

        V0,Vth,V_Rs,inflectionline=manual_inflection(Vg,G,Gsmooth,smoothing,
                                                    float(Vg_inflman.get()),float(dGdVg_inflman.get())*1e-3)

        paramdict['V0 (V)']=V0
        paramdict['Vth (V)']=Vth
        paramdict['Vg_infl (V)']=float(Vg_inflman.get())
        paramdict['V_Rs (V)']=V_Rs
                                
        set_exportparams()

        exportdatadict['Inflection fit (S)']=inflectionline
        set_exportdata()

        if GorI.get()=='G provided':
            if Gunits.get()=='2e2/h':
                inflectionline=inflectionline/7.748091729e-5
            elif Gunits.get()=='e2/h':
                inflectionline=inflectionline/7.748091729e-5*2
        ax[0].plot(Vg,inflectionline,label='slope at inflection')
        ax[0].legend()
        fig1.canvas.draw()

    Label(derivextraframe,text='Enter peak position manually:').grid(column=0,row=0)
    Label(derivextraframe,text='Vg_infl').grid(column=1,row=0)
    Vg_inflman=StringVar()
    Vg_inflmanentry=Entry(derivextraframe,textvariable=Vg_inflman,width=5).grid(column=2,row=0)
    Label(derivextraframe,text='V').grid(column=3,row=0)
    Label(derivextraframe,text='dG/dVg_infl').grid(column=4,row=0)
    dGdVg_inflman=StringVar()
    dGdVg_inflmanentry=Entry(derivextraframe,textvariable=dGdVg_inflman,width=5).grid(column=5,row=0)
    Label(derivextraframe,text='(ms/V)').grid(column=6,row=0)
    manualinflbutton=Button(derivextraframe,text='Enter',command=plot_manual_inflection)
    manualinflbutton.grid(column=7,row=0)
    
    
    
    ## Window to plot final mobility vs density plot
    mobframe = ttk.Frame(root, padding='3 3 12 3')
    mobframe.grid(column=1,row=1)
    mobtopframe = ttk.Frame(mobframe)
    mobtopframe.grid(row=0)
    mobbottomframe=ttk.Frame(mobframe)
    mobbottomframe.grid(row=1)
    fig3 = Figure(figsize=figuresize, dpi=100)
    canvas3 = FigureCanvasTkAgg(fig3, master=mobbottomframe)  # A tk.DrawingArea.
    canvas3.draw()
    canvas3.get_tk_widget().pack()#(side=TOP, fill=BOTH, expand=1)
    toolbar = NavigationToolbar2Tk(canvas3, mobbottomframe)
    toolbar.update()
    canvas3.get_tk_widget().pack()#(side=TOP, fill=BOTH, expand=1)
    def plot_mobility():
        fig3.clf()
        
        Vg,G,holes=VgandG(convertunits=True)
            
        Rs,Rs_fit,V_Rs_ind,result_drudeRs=perform_Rs_fit(Vg,G,
                                                        paramdict['V0 (V)'],paramdict['V_Rs (V)'],
                                                        float(initial_Rs.get()),float(initial_mu.get())*1e-4,
                                                        float(L.get()),float(c.get()),holes)
        
        paramdict['Rs (Ohm)']=Rs             
        set_exportparams()
        
        if holes:
            exportdatadict['Vg for Rs fit (V)']=Vg[:V_Rs_ind]
        else:
            exportdatadict['Vg for Rs fit (V)']=Vg[V_Rs_ind:]
        exportdatadict['Rs fit (S)']=Rs_fit
        
        if Ctype.get()=='Width':
            Cperarea=float(c.get())/(float(L.get())*float(W.get()))
        else:
            Cperarea=float(W.get())
        #density=float(c.get())*(Vg-paramdict['V0 (V)'])/(float(L.get())*float(W.get())*1.6e-19)
        if holes:
            density=Cperarea*(paramdict['V0 (V)']-Vg)/1.602176634e-19
            mu_eff=float(L.get())**2/(float(c.get())*(paramdict['V0 (V)']-Vg)*((1/G)-Rs))
        else:
            density=Cperarea*(Vg-paramdict['V0 (V)'])/1.602176634e-19
            mu_eff=float(L.get())**2/(float(c.get())*(Vg-paramdict['V0 (V)'])*((1/G)-Rs))
        exportdatadict['density (1/m2)']=density
        exportdatadict['mu_eff (m2/Vs)']=mu_eff
        
        set_exportdata()
        plotstart=(np.abs(Vg - (2*paramdict['Vth (V)']-paramdict['Vg_infl (V)']))).argmin()
            
        ax[2]=fig3.add_subplot()
        if holes:
            ax[2].plot(density[:plotstart]*1e-12/1e4,mu_eff[:plotstart]*1e4,'k',label='mu_eff')
        else:
            ax[2].plot(density[plotstart:]*1e-12/1e4,mu_eff[plotstart:]*1e4,'k',label='mu_eff')
        ax[2].set_xlabel('Carrier density x 10$^{12}$ (cm$^{-2}$)')
        ax[2].set_ylabel('Mobility (cm$^2$/(Vs))')
        ax[2].legend()
        fig3.tight_layout()
        fig3.canvas.draw()
        
        if GorI.get()=='G provided':
            if Gunits.get()=='2e2/h':
                Rs_fit=Rs_fit/7.748091729e-5
            elif Gunits.get()=='e2/h':
                Rs_fit=Rs_fit/7.748091729e-5*2
        
        ax[0].plot(exportdatadict['Vg for Rs fit (V)'],Rs_fit,label='fit for R_s')
        ax[0].legend()
        fig1.canvas.draw()
        
    def clear_mobility():
        fig3.clf()
        fig3.canvas.draw()
    Rsbutton=Button(mobtopframe, text='5) Find Rs and plot mu_eff/Refresh', command=plot_mobility,fg='green')
    Rsbutton.grid(row=0)
    Button(mobtopframe, text='Clear plot', command=clear_mobility).grid(row=0,column=2)
    CreateToolTip(Rsbutton,'If the fit is not working well, and the fit in the first panel (orange line) seems '
                'bad, try checking the initial fit values for R_s and mu')
    def plot_drude():
        
        Vg,G,holes=VgandG(convertunits=True)
            
        mu_drude,drude_fit,Rs_drude,Vth_ind,result_drude=perform_drude_fit(Vg,G,
                                                        paramdict['Vth (V)'],
                                                        float(initial_Rs.get()),float(initial_mu.get())*1e-4,
                                                        float(L.get()),float(c.get()),holes)
        
        paramdict['mu_FET (m2/Vs)']=mu_drude
        #paramdict['Rs_drude (Ohm)']=Rs_drude           
        set_exportparams()
        
        if holes:
            exportdatadict['Vg for mu_FET fit (V)']=Vg[:Vth_ind]
        else:
            exportdatadict['Vg for mu_FET fit (V)']=Vg[Vth_ind:]
        exportdatadict['mu_FET fit (S)']=drude_fit
        set_exportdata()
        
        mu_drude_array=np.full(Vg.shape[0],mu_drude) # Make an array of the correct size for plotting
        plotstart=(np.abs(Vg - (2*paramdict['Vth (V)']-paramdict['Vg_infl (V)']))).argmin()
        if holes:
            ax[2].plot(exportdatadict['density (1/m2)'][:plotstart]*1e-12/1e4,mu_drude_array[:plotstart]*1e4,label='mu_FET')
        else:
            ax[2].plot(exportdatadict['density (1/m2)'][plotstart:]*1e-12/1e4,mu_drude_array[plotstart:]*1e4,label='mu_FET')
        ax[2].legend()
        fig3.canvas.draw()
        
        if GorI.get()=='G provided':
            if Gunits.get()=='2e2/h':
                drude_fit=drude_fit/7.748091729e-5
            elif Gunits.get()=='e2/h':
                drude_fit=drude_fit/7.748091729e-5*2
        
        ax[0].plot(exportdatadict['Vg for mu_FET fit (V)'],drude_fit,label='mu_FET fit')
        ax[0].legend()
        fig1.canvas.draw()
        
    drudebutton=Button(mobtopframe, text='Fit and plot mu_FET', command=plot_drude)
    drudebutton.grid(row=0,column=1)
    exporterrorframe=ttk.Frame(root, padding='3 3 12 3')
    exporterrorframe.grid(column=2,row=1)
    ### Export frame, including scrollable windows of the produced data and parameters
    exportframe = ttk.Frame(exporterrorframe, padding='3 3 12 3')
    exportframe.grid(column=0,row=0)
    def export_data():
        jsonexportdata={}
        for param in exportdatadict:
            jsonexportdata[param]=list(exportdatadict[param])
        filename = tk.filedialog.asksaveasfilename(title='Select file name and type.',
                                                defaultextension='.json',filetypes=[('JSON (*.json)','*.json'),('CSV (*.csv)','*.csv')])
        if '.json' in filename:
            with open(filename, 'w', encoding='utf-8') as f:
                json.dump(jsonexportdata, f, ensure_ascii=False,indent=4)
        elif '.csv' in filename:
            with open(filename, 'w', newline='') as f:
                writer = csv.writer(f)
                writer.writerow([param for param in exportdatadict])
                for i in range(len(exportdatadict['Vg (V)'])):
                    row = []
                    for param in exportdatadict:
                        try:
                            row.append(exportdatadict[param][i])
                        except IndexError:
                            row.append('')
                    writer.writerow(row)
        else:
            print('File type must be .json or .csv')
    exportdatabutton=Button(exportframe, text='6) Export data', command=export_data, fg='green')
    exportdatabutton.grid(row=0,column=0,stick='E')
    CreateToolTip(exportdatabutton,'Export data as json or csv. Reimporting a json to python with json.load() returns a dictionary')
    def copy_data():
        maxsize=1000
        for param in exportdatadict:
            if len(exportdatadict[param])>maxsize:
                maxsize=len(exportdatadict[param])
        np.set_printoptions(threshold=maxsize)
        root.clipboard_clear()
        datacopy=''
        for param in exportdatadict:
            datacopy+=param+': '+str(exportdatadict[param])+'\n'
        root.clipboard_append(datacopy)
        np.set_printoptions(threshold=1000)
    copydatabutton=Button(exportframe, text='Copy data', command=copy_data)
    copydatabutton.grid(row=0,column=1,sticky='W')
    def export_params():
        filename = tk.filedialog.asksaveasfilename(title='Select file name and type.',
                                                defaultextension='.json',filetypes=[('JSON (*.json)','*.json'),('Comma separated values (*.csv)','*.csv')])
        if '.json' in filename:
            with open(filename, 'w', encoding='utf-8') as f:
                json.dump(paramdict, f, ensure_ascii=False, indent=4)
        elif '.csv' in filename:
            with open(filename, 'w', newline='') as f:
                writer = csv.writer(f)
                for param in paramdict:
                    writer.writerow([param,paramdict[param]])
    Button(exportframe, text='7) Export parameters', command=export_params, fg='green').grid(row=0,column=3)
    def copy_params():
        root.clipboard_clear()
        paramcopy=''
        for param in paramdict:
            paramcopy+=param+': '+str(paramdict[param])+'\n'
        root.clipboard_append(str(paramcopy))
    Button(exportframe, text='Copy params', command=copy_params).grid(row=0,column=4,columnspan=2)
    databox=Listbox(exportframe,listvariable=(exportdatavar),height=11,width=35)
    databox.grid(column=0,row=1,columnspan=2,sticky=('N W E S'))
    s1 = ttk.Scrollbar(exportframe, orient=VERTICAL, command=databox.yview)
    s1.grid(column=2,row=1,sticky=('N S'))
    databox['yscrollcommand'] = s1.set
    s2 = ttk.Scrollbar(exportframe, orient=HORIZONTAL, command=databox.xview)
    s2.grid(column=0,row=2,columnspan=2,sticky=('W E'))
    databox['xscrollcommand'] = s2.set
    parambox=Listbox(exportframe,listvariable=(exportparamsvar),height=11,width=30)
    parambox.grid(column=3,row=1,columnspan=2,sticky=('N W E S'))
    s3 = ttk.Scrollbar(exportframe, orient=VERTICAL, command=parambox.yview)
    s3.grid(column=5,row=1,sticky=('N S'))
    parambox['yscrollcommand'] = s3.set
    s3 = ttk.Scrollbar(exportframe, orient=HORIZONTAL, command=parambox.xview)
    s3.grid(column=3,row=2,columnspan=2,sticky=('W E'))
    parambox['xscrollcommand'] = s3.set
    # errorframe= ttk.Frame(exporterrorframe, padding='3 3 12 12')
    # errorframe.grid(column=0,row=1)
    # errorbox=Text(errorframe,height=6,width=50)
    # errorbox.grid(column=0,row=0,sticky=('N W E S'))
    # s4=ttk.Scrollbar(errorframe,orient=VERTICAL,command=errorbox.yview)
    # s4.grid(column=1,row=0,sticky=('N S'))
    # errorbox['yscrollcommand']=s4.set
    # s5=ttk.Scrollbar(errorframe,orient=HORIZONTAL,command=errorbox.xview)
    # s5.grid(column=0,row=1,sticky=('W E'))
    # errorbox['xscrollcommand']=s5.set
    root.mainloop()