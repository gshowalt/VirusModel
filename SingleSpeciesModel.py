import numpy as np
  
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib import cm
import matplotlib.tri as tri
from matplotlib.colors import LogNorm

from collections import Counter

from functools import wraps

import csv
import sys

import itertools
from itertools import islice, cycle, chain

from scipy.interpolate import griddata
from scipy import interpolate
from scipy.integrate import odeint

import decimal as dc
from decimal import Decimal

from bokeh.io import curdoc, push_notebook, show, output_notebook
from bokeh.layouts import column, row, layout, widgetbox
from bokeh.models import ColumnDataSource, Slider, TextInput, Range1d, HoverTool, LinearColorMapper, BoxAnnotation, Toggle
from bokeh.models import BoxSelectTool, BoxZoomTool, LassoSelectTool, CustomJS, InputWidget
from bokeh.plotting import figure, show, output_file
from bokeh.palettes import viridis, plasma, brewer
from bokeh.transform import transform
from bokeh.models.widgets import Button, Toggle

import altair as alt

import csv
import sys

import seaborn as sns
import pandas as pd
import math
from ipywidgets import interact

runs = 10
time = 100
temp_list = (range(-15,-2,1))

mux = 1
betx = 1
phix = 1
gamx = 1


"""WITHOUT AN INFECTED CLASS"""
def f(s,t, beta, mu, phi, gamma):
    alpha = 1.2e-7 #4.2e-7 at +8, or 1.2e-7 at lower temps, at -5 --> mu = 0.25/day = 0.01/hr = 1e-8
    # alpha is a coefficient that we'd like to change with temperature? Or change eta?
    lamb = 0.001 #/(0.1/mu)
    #latent period/"lysis rate" in per hour
    #nutrient transfer coefficient to bacteria (ug/cell * hr)
    Q = 0.022
    #half saturation constant (ug/mL)
    #eta = 1e5
    #conversion rate between uptake and production (cells/ug)
    d = 0.002
    #constant of bacterial death (1/hr)
    m = 0.0001
    #constant of viral decay (1/hr)
    g = 0.2
    #POM transfer coefficient from bacteria (ug/cell*hr)
    n = 0.99
    #POM transfer coefficient from viral lysis ug/[burst]cell
    #gamma is a lysogeny value
    N = s[0]
    B = s[1]
    V = s[2]
    P = s[3]
    #systems of equations below
    if N < 0:
        N = 0
    if B < 1:
        B = 1
    if V < 1:
        V = 1
    gamma = (gamma)/(mu+0.1)
    dNdt = - alpha * (N / (N + Q)+ P/(Q+P)) * B + g * (alpha  * (N/(N+Q))*B) + (n * 1e-7 * (gamma) * V * B)
    if N < 0:
        N = 0
    #nutrient term
    dBdt = (mu) * (N/(Q + N)) * B - gamma * phi * V * B - d*B
    if B < 1:
        B = 1
    dVdt =  gamma*beta*phi*V* B  - phi * V * B -  m*V
    if V < 1:
        V = 1
    #virus term
    #dPdt = (g * (0.0083*1e-7))*B + (n * 1e-7 * phi * V * B*RCR) + 1e-10*m*V + 1.0e-7*d*B - (P/(P+Q))*alpha * B
    dPdt = g * alpha  * (N/ (N+Q))*B + n * 1e-7 * (1-gamma)*phi*B*V
    #POM term
    return [dNdt, dBdt, dVdt, dPdt]

def equation_run(temp_list, mu_x, beta_x, phi_x, gamma_x, time, runs):
    nuts_val = 0.1
    nuts_init = nuts_val
    t = np.linspace(1,time,1000)
    DOM = []
    NUTS = []
    BAC = []
    INF = []
    VIR = []
    DOM_End = []
    RATIO = []
    VIRendlist = []
    BACendlist =[]
    INFendlist = []
    Totalendlist = []
    DOM_Sum = []
    gamma_list2 = []
    Bac_init_lo = 1e4
    Bac_init_hi = 1e6
    Vir_init_lo = 1
    Vir_init_hi = 100
    beta_list2 = []
    phi_list2 = []
    mu_list2 = []


    for j in range(0, len(temp_list)):
        temp = temp_list[j]

        if temp < -1:
            RCR = 0.0716*temp**4 + 2.9311*temp**3 + 34.108*temp**2 + 45.826*temp + 3.5125 #Fit from Wells and Deming, 2006
            BCF = -0.0106 * temp **2 - 0.519 * temp + 0.2977
            sal = 32 * BCF
        else:
            RCR = 1
            BCF = 1
            sal = 32

        # define the changing parameters (make lists through which to iterate)
        beta_list = []
        mu_list = []
        phi_list = []
        gamma_list = []

        beta_mu = 0.0064 * temp**3 - 0.3047 * temp ** 2 + 0.7701 * temp + 93.605
        beta_std = 0.0095 * temp **3 - 0.5184 * temp**2 + 2.2456 * temp + 126.59
        if beta_std < 0:
            beta_std = 0.

        mu_mu = 2e-5*temp**3 + 0.0008 * temp **2 + 0.0091 * temp + 0.0386
        mu_std = 2e-5*temp**3 + 0.0009 * temp **2 + 0.0144 * temp + 0.0818
        if mu_std<0:
            mu_std = 0.001

        #phi as function of temperature    
        phi_mu = 6e-13 * temp **5 - 2e-11 * temp ** 4 + 1e-10 * temp ** 3 + 3e-9 * temp ** 2 - 3e-8 * temp + 5e-8
        phi_std = 1e-13 * temp ** 5 - 1e-12 * temp ** 4 - 5e-11 * temp ** 3 + 1e-10 * temp ** 2 + 5e-9 * temp + 2e-8
        #phi as function of salt
        #phi_mu = -1e-11*sal**2 +4e-9*sal - 9e-8
        #phi_std = -2e-11*sal**2 + 4e-9*sal - 9e-8
        if phi_std < 0:
            phi_std = 0
    # RUN THROUGH THE NUMBER OF RUNS
    #first make a list of the temperature-depdendant parameters so it exists in entirety
        for i in range(0, runs):
            beta = beta_x*np.random.normal(beta_mu, beta_std)
            if beta < 0:
                beta = 1
            beta_list.append(beta)
            beta_list2.append(beta)

            mu = mu_x*np.random.normal(mu_mu, mu_std)
            if mu <= 0:
                mu = 0.001
            mu_list.append(mu)
            mu_list2.append(mu)

            phi = phi_x*np.random.normal(phi_mu, phi_std)
            phi = phi * RCR
            if phi < 0:
                phi = 0
            phi_list.append(phi)
            phi_list2.append(phi/RCR)

            gamma = gamma_x*mu
            if gamma > 1:
                gamma = 1
            gamma_list.append(gamma)
            gamma_list2.append(gamma)

        for i in range(0, runs):
            #establish run-dependant variables
            nuts_val = nuts_init * BCF #seawater values * Brine concentrating factor
            B_init = np.random.randint(Bac_init_lo, Bac_init_hi) * BCF #seawater values * brine concentrating factor
            V_init = np.random.randint(Vir_init_lo, Vir_init_hi) * B_init #Bacterial concetration times number of viruses (Ratio)       
       

            #initial conditions and integration
            s0=[nuts_val, B_init, V_init,0]
            s = odeint(f,s0,t, args = (beta_list[i], mu_list[i], phi_list[i], gamma_list[i]))

            #manipulate data into form wanted
            NUTS = NUTS + list(x for x in s[:,0])
            BAC = BAC + list(x for x in s[:,1])
            VIR = VIR + list(x for x in s[:,2])
            DOM = DOM +  list(x for x in s[:,3])
            DOM_End.append(s[-1,3])
            DOM_Sum.append(sum(DOM))
            VIRend = s[-1,2]
            BACend = s[-1,1]
            VIRendlist.append(VIRend)
            BACendlist.append(BACend)
            Totalendlist.append((BACend))
    return [Totalendlist, VIRendlist, gamma_list2]

""" Run the functions"""
# --------  Run Model and Plot it -------- #
[Totalendlist, VIRendlist, gamma_list2] =  equation_run(temp_list, mux, betx, phix, gamx, time, runs)

temp_list2=list(itertools.chain.from_iterable(itertools.repeat(x, runs) for x in temp_list))
size_list = [i*500 for i in gamma_list2]
lytic = [i*100 for i in gamma_list2]
source = ColumnDataSource(data=dict(x=Totalendlist, y=VIRendlist, temp_list2 = temp_list2, size_list = size_list, lytic = lytic))

hover = HoverTool(tooltips=[("Ice Temp (C)","@temp_list2"),("Bacterial Population", "$x"),("Viral Population", "$y"), ("Lytic %", "@lytic")])
mapper = LinearColorMapper(palette=viridis(256), low=min(temp_list2), high=max(temp_list2))

p = figure(plot_height=500, plot_width=500, x_axis_type = 'log', y_axis_type = 'log',x_range = (1e4, 1e8), y_range = (1e5,1e10), tools = [hover])
p.xaxis.axis_label = 'Bacteria per mL'
p.yaxis.axis_label = 'Virus per mL'


# plot the ratio lines 
Ratiolines = np.linspace(1e4, 1e8,1000)
Lines = p.line(Ratiolines, 1*Ratiolines)
Lines2 = p.line(Ratiolines, 10*Ratiolines)
Lines3 = p.line(Ratiolines, 100*Ratiolines)
Lines4 =p.line(Ratiolines, 1000*Ratiolines)
Lines5= p.line(Ratiolines, 10000*Ratiolines)
toggle2 = Toggle(label="Ratio Lines", button_type="success", active=True)
toggle2.js_link('active', Lines , 'visible')
toggle2.js_link('active', Lines2 , 'visible')
toggle2.js_link('active', Lines3 , 'visible')
toggle2.js_link('active', Lines4 , 'visible')
toggle2.js_link('active', Lines5 , 'visible')

#plot the data
p.circle('x', 'y', fill_color= transform('temp_list2', mapper), size = 'size_list', fill_alpha=0.6, line_color=None, source = source)
layout = row(
    p, toggle2
)


show(layout)
