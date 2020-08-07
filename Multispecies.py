#import everything here
import numpy as np
import pandas as pd
from pandas import DataFrame
import math

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.tri as tri
import matplotlib.colors as mcolors
from matplotlib.colors import LogNorm
from random import paretovariate
import random

from scipy.integrate import odeint
from collections import Counter
import csv
import sys
import itertools
from itertools import islice, cycle, chain
from scipy.interpolate import griddata
from scipy import interpolate
from datetime import datetime
import os

import seaborn as sns
import altair as alt

from bokeh.io import curdoc
from bokeh.layouts import column, row, layout, widgetbox
from bokeh.models import ColumnDataSource, Slider, TextInput, Range1d, HoverTool, LinearColorMapper, BoxAnnotation, Toggle
from bokeh.models import BoxSelectTool, BoxZoomTool, LassoSelectTool, CustomJS, InputWidget
from bokeh.plotting import figure, show, output_file
from bokeh.palettes import viridis, plasma, brewer
from bokeh.transform import transform
from bokeh.models.widgets import Button, Toggle



# -------- Define Classes -------- #

class GenoB(object):
    def __init__(self, hi, pop, delt, infect, mu):
        self.hi = hi
        self.pop = pop
        self.delt = delt
        self.infect = infect
        self.mu = mu

class GenoV(object):
    def __init__(self, hi, burst, adsorp,  pop, delt,infect, decay,s, loss):
        self.hi = hi
        self.burst = burst
        self.adsorp = adsorp
        self.pop = pop
        self.delt = delt
        self.infect = infect
        self.decay = decay   
        self.s = s
        self.loss = loss

# -------- Fitess Landscape ------- #
def Fit(h):
    if h>=1:
        delt = 0
    if h <= 0:
        delt = 0
    if h>0 and h<1:
#        delt = float(np.random.normal(1,10))*h      
        delt = 0.4*h+1.2
        delt = 1
        if delt < 0:
            delt = 0
    return delt

# --------------------------------- Evolution/population model below ------------------------------------------#

def round_nearest(x, a):
    x = float(x)
    return round(round(x / a) * a, -int(math.floor(math.log10(a))))


# -------- Seed Function -------- #
def Seed(x,y,xmin, xmax, ymin, ymax, temp, mu_x, phi_x):
  
    genotypesB = []
    genotypesV = []

    if temp < -1:
        RCR = 0.0716*temp**4 + 2.9311*temp**3 + 34.108*temp**2 + 45.826*temp + 3.5125 #Fit from Wells and Deming, 2006
        BCF = -0.0106 * temp **2 - 0.519 * temp + 0.2977
        sal = 32 * BCF
    else:
        RCR = 1
        BCF = 1
        sal = 32

    beta_mu = 0.0064 * temp**3 - 0.3047 * temp ** 2 + 0.7701 * temp + 93.605
    beta_std = 0.0095 * temp **3 - 0.5184 * temp**2 + 2.2456 * temp + 126.59
    if beta_std < 0:
        beta_std = 0.

    mu_mu = 2e-5*temp**3 + 0.0008 * temp **2 + 0.0091 * temp + 0.0386
    mu_std = 2e-5*temp**3 + 0.0009 * temp **2 + 0.0144 * temp + 0.0818
    if mu_std<0:
        mu_std = 0.001

        #phi as function of temperature    
        #phi_mu = 6e-13 * temp **5 - 2e-11 * temp ** 4 + 1e-10 * temp ** 3 + 3e-9 * temp ** 2 - 3e-8 * temp + 5e-8
        #phi_std = 1e-13 * temp ** 5 - 1e-12 * temp ** 4 - 5e-11 * temp ** 3 + 1e-10 * temp ** 2 + 5e-9 * temp + 2e-8
        #phi as function of salt

    phi_mu = -1e-11*sal**2 +4e-9*sal - 9e-8
    phi_std = -2e-11*sal**2 + 4e-9*sal - 9e-8
    if phi_std < 0:
        phi_std = 0



    for i in range(0,x):
        h = float((np.random.randint(100,1000)))*1/1000
        pop = np.random.pareto(1) * BCF * ((xmax+xmin)/ 2)
        delt = Fit(h)
        mu = mu_x*np.random.normal(mu_mu, mu_std)
        if mu <= 0:
            mu = 0.001
        i = GenoB(h, pop, delt,0,mu)
        genotypesB.append(i)

    for i in range(0,y):
        h = float((np.random.randint(100,1000)))*1/1000
        pop = np.random.pareto(1) * BCF * ((ymax+ymin)/ 2)
        delt = Fit(h)
        adsorp =  np.random.normal(phi_mu, phi_std)
        adsorp = phi_x * adsorp * RCR
        burst = np.random.normal(beta_mu, beta_std)
        if burst < 0:
             burst = 1
        if adsorp < 0:
            adsorp = 0
        decay = 1e-8 #(np.random.randint(10, 2600) * (1/10000))
        loss = 0
        s = 1000
        i = GenoV(h,  burst, adsorp, pop,  delt,0, decay, s, loss)
        genotypesV.append(i)
   

    return genotypesB, genotypesV


# -------- Bin Function -------- #
def Bin_same(tol, MasterList):
#curently tolerance is set to the same for virus and bacteria, would have to update otherwise

    for x in MasterList:
        for i in range(0, len(x)):
            for j in range(i+1, len(x)):
#                print "i = ", i, "j = ", j
                if x[i].hi >= x[j].hi - tol and x[i].hi <= x[j].hi + tol :
                    x[i].pop = x[i].pop + x[j].pop
                    x[j].hi = 0
                    x[j].pop = 0
            #remove duplicate items
        for item in list(x):
            if item.hi <= 0 or item.pop <= 0 or item.delt == 0:
                x.remove(item)
#               for i in range(0, len(x)): print x[i].hi, x[i].pop
    return MasterList

#Masterlist = Bin(0.01, MyList)

# -------- Infect Function -------- #
def Infect(genotypesB, genotypesV):
#calculates the infection parameters for each genotype
    for i in range(0, len(genotypesB)):
        infect_sum = 0
        for j in range(0, len(genotypesV)):
            gamma = (genotypesB[i].mu/(genotypesB[i].mu + 0.1))
            infect_sum = infect_sum + gamma * genotypesV[j].adsorp * (np.exp(-genotypesV[j].s * ((genotypesB[i].hi - genotypesV[j].hi))**2) * genotypesB[i].pop * genotypesV[j].pop)
            if infect_sum <= 0:
                infect_sum = 0
        genotypesB[i].infect = infect_sum

    for j in range(0, len(genotypesV)):
        infect_sum = 0
        gammaless_infect = 0
        for i in range(0, len(genotypesB)):
            gamma = (genotypesB[i].mu/(genotypesB[i].mu + 0.1))
            infect_sum = infect_sum + gamma * genotypesV[j].adsorp * (np.exp(-genotypesV[j].s * ((genotypesB[i].hi - genotypesV[j].hi))**2) * genotypesB[i].pop * genotypesV[j].pop)
            gammaless_infect = gammaless_infect + genotypesV[j].adsorp * (np.exp(-genotypesV[j].s * ((genotypesB[i].hi - genotypesV[j].hi))**2) * genotypesB[i].pop * genotypesV[j].pop)
            if infect_sum <= 0:
                infect_sum  = 0
        genotypesV[j].infect = infect_sum
        genotypesV[j].loss = gammaless_infect
    
    
    return genotypesB, genotypesV

# ------- Mutation Function -------- #
def Mutate(tol,genotypesB, genotypesV):
#calculates in the number of mutants for each instance of class and assigns new genotypes 
    #MB = 10e-6 #mutation rate; start with 10e6
    #MV = 10e-6 #mutation rate; start with 10e-6
    MB = 0
    MV = 0
    mutants_B = []
    mutants_V = []

    #calculate mutant population and genotype of bacteria
    for i in range(0, len(genotypesB)):
        mutBpop = MB * genotypesB[i].mu * (N)/(Q + N) * (genotypesB[i].pop * genotypesB[i].delt)
        mub = np.round(np.random.normal(0.0,0.05,1), 2)  
        mut_hi = genotypesB[i].hi + mub
        delt = Fit(mut_hi)
        mutmu = genotypesB[i].mu
        new = GenoB(mut_hi,mutBpop, delt, 0, mutmu)
        mutants_B.append(new)

    #bin mutant bacteria with rest of population
    for i in range(0, len(genotypesB)):
        for j in range(0, len(mutants_B)):
            if genotypesB[i].hi >= mutants_B[j].hi - tol and genotypesB[j].hi <= mutants_B[j].hi + tol :
                genotypesB[i].pop = genotypesB[i].pop + mutants_B[j].pop
                mutants_B[j].hi = 0
                mutants_B[j].pop = 0

    #remove duplicate items
    for item in list(mutants_B):
        if item.hi <= 0 or item.pop <= 0 or item.delt ==0:
            mutants_B.remove(item)
        else:
            genotypesB.append(item)
            mutants_B.remove(item)
#calculate mutant population and genotype of viruses
    for i in range(0, len(genotypesV)):
        mutVpop = MV * genotypesV[i].burst * genotypesV[i].infect
        muv = np.round(np.random.normal(0.0,0.05,1), 2)
        mut_hi = genotypesV[i].hi + muv
        delt = Fit(mut_hi)
        mut_burst = np.random.randint(5,100)
        mut_adsorp = (np.random.randint(1,10) * (1e-10))
        mut_decay = (np.random.randint(10, 2600) * (1/10000))
        mut_gam = 1#np.random.randint(0,10)*0.1
        mut_s = np.random.normal(250, 10, 1)
        new = GenoV(mut_hi,mut_burst, mut_adsorp, mutVpop, delt , 0, mut_decay, mut_gam, mut_s)
        mutants_V.append(new)

    #bin mutant virus with rest of population
    for i in range(0, len(genotypesV)):
        for j in range(0, len(mutants_V)):
            if genotypesV[i].hi >= mutants_V[j].hi - tol and genotypesV[j].hi <= mutants_V[j].hi + tol :
                genotypesV[i].pop = genotypesV[i].pop + mutants_V[j].pop
                mutants_V[j].hi = 0
                mutants_V[j].pop = 0

    #remove duplicate items
    for item in list(mutants_V):
        if item.hi <= 0 or item.pop <= 0 or item.delt == 0:
            mutants_V.remove(item)
        else:
            genotypesV.append(item)
            mutants_V.remove(item)

    MyList = [genotypesB, genotypesV]
    [genotypesB, genotypesV] = Bin_same(tol, MyList)

    return genotypesB, genotypesV


# -------- Population update -------- #
def PopUpdate(t, genotypesB, genotypesV, N, P):
    for i in range(0, len(genotypesB)):
        dNdt = - alpha * (N / (N + Q)) * (genotypesB[i].pop * genotypesB[i].delt) + g1* alpha * (N/(N+Q)) * genotypesB[i].pop  + n2 * (1e-7)* genotypesB[i].infect
        dBdt = genotypesB[i].mu  * (N/(Q + N)) * (genotypesB[i].pop * genotypesB[i].delt) - genotypesB[i].infect - d * genotypesB[i].pop
        dPdt = g1* alpha * (N/(N+Q)) * genotypesB[i].pop  + n2 * (1e-7)* genotypesB[i].infect
        N = N + dNdt * t
        P = P + dPdt * t
        if N < 0:
            N = 0
        if P < 0:
            P = 0
        #print ('nutrient concentation is', N)
        #print ('recycled change', dPdt)
        #print ('recycled total', P)
        genotypesB[i].pop = genotypesB[i].pop +  dBdt * t
        if genotypesB[i].pop < 0:
            genotypesB[i].pop = 1
        if N < 0:
            N = 0
        if P < 0:
            P = 0
       


    for j in range(0, len(genotypesV)):
        dVdt = genotypesV[j].burst * genotypesV[j].infect  - genotypesV[j].loss - genotypesV[j].decay*genotypesV[j].pop
        genotypesV[j].pop = genotypesV[j].pop + dVdt * t
        if genotypesV[j].pop < 0:
            genotypesV[j].pop = 1   
#        print ("the type of genotype data index", j , "is", type(genotypesV[j].hi))
#        vpoplist.append(int(genotypesB[i].pop))
#        vlist.append(float(genotypesB[i].hi))


    return genotypesB, genotypesV, N, P


    
#-------- Write hard-coded parameters -------- #

temperature = -6

alpha = 1.2e-7*3**((temperature-23)/10) # changes with temperature following q10 law
print (alpha)
Q = 0.022 #0.022
d = 0.0000001
#m = 0.0000001
g1 = 0.1
n2 = 0.99
#s = 100
#phi = 1e-9

tolerance = 0.01
time = 2500
step = 1
runs = 10 #be careful in balancing these parameters; at -6, 2500 hours, 10 runs = ~2 hours of run time
timepoints = np.arange(0,time+step,step)



# -------- Initialize populations and bin at given tolerance level --------- #
#tolerance = 0.005

templist = np.arange(-6,-5,1)





N_list = []
P_list = []
BACendlist = []
VIRendlist = []
TEMPendlist = []
GAMendlist = []
BACinit = []
VIRinit = []
Pop_list = []
popsum = 0
mu_x = 0.1
phi_x = 1e-5
s = "750"
PopB_list = []
PopV_list = []
Shan_list = []
Simp_list = []
PopB_list_internal = np.zeros(shape = (runs,time+1))
PopV_list_internal  = np.zeros(shape = (runs,time+1))
VBR_list_internal  = np.zeros(shape = (runs,time+1))
Shannon_list  = np.zeros(shape = (runs,time+1))
P_list_internal = np.zeros(shape = (runs,time+1))
Shan_list_internal  = []
Simp_list_internal  = []


for x in range(0,runs):
    #populations = []
    genotypes = np.arange(0,1+tolerance, tolerance)
    dfB = pd.DataFrame(1, index = timepoints , columns = genotypes)
    dfV = pd.DataFrame(1, index = timepoints , columns = genotypes)
    dfB.columns.name = 'genotypes'
    dfB.index.name = 'timelist'
    dfV.columns.name = 'genotypes'
    dfV.index.name = 'timelist'

    if temperature < -1.8:
        BCF = -0.0106 * temperature **2 - 0.519 * temperature + 0.2977
    else:
        BCF = 1

    N = 0.12*BCF
    P = 0

    # Seed the population and clear lists from last temp
    [genotypesB, genotypesV]= Seed(100,100,1e1,5e3,1e1,5e4,temperature, mu_x, phi_x)
    MyList = [genotypesB, genotypesV]
    [genotypesB1, genotypesV1] = Bin_same(tolerance,MyList)
    viruspopsum = 0
    bactopopsum = 0

    t = 0


    #update populations 
    [genotypesB1, genotypesV1] = Bin_same(tolerance,MyList)
    for t in timepoints:
        popBsum = 0
        popVsum = 0
        populations = []
        shannon_inst = 0
        simpsons_inst = 0
        [genotypesB, genotypesV] = Infect(genotypesB1, genotypesV1)
        [genotypesB, genotypesV, N, P] = PopUpdate(t, genotypesB, genotypesV, N, P)
        N_list.append(N)
        P_list.append(P) 
        [genotypesB, genotypesV] = Mutate(tolerance,genotypesB, genotypesV)
        [genotypesB, genotypesV] = Bin_same(tolerance,MyList)

        # sum up all the bacterial populations for plot over time
        for i in range(0, len(genotypesB)):
            gen_cat = round_nearest(genotypesB[i].hi, tolerance)
            if gen_cat in dfB.columns:
                if genotypesB[i].pop <=1:
                    genotypesB[i].pop = 1
                dfB.loc[t, gen_cat] = genotypesB[i].pop
            if genotypesB[i].pop > 5:
                populations.append(float(genotypesB[i].pop))
                popBsum = popBsum + float(genotypesB[i].pop)


        # sum up all the virus populations for plot over time
        for i in range(0, len(genotypesV)):
            gen_cat = round_nearest(genotypesV[i].hi, tolerance)
            if gen_cat in dfV.columns:
                if genotypesV[i].pop <=1:
                    genotypesV[i].pop = 1
                dfV.loc[t, gen_cat] = genotypesV[i].pop
            if genotypesV[i].pop > 10:
                popVsum = popVsum + float(genotypesV[i].pop)

        # calculate diversity indices
        diversity_sum = sum(populations) 
        for i in range(0, len(populations)):
            fraction = populations[i]/diversity_sum
            shannon_inst += -1*(fraction) * np.log((fraction))
            #simpsons_inst = (populations[i]/diversity_sum)**2
        
        PopB_list_internal[x, t] = popBsum
        PopV_list_internal[x, t] = popVsum
        VBR_list_internal[x, t] = popVsum/(popBsum+1)
        P_list_internal[x,t] = P

        Shannon_list[x, t] = shannon_inst
        #Simp_list_internal.append(1/simpsons_inst)
        
        

    Shan_list.append(shannon_inst)
    #Simp_list.append(1/simpsons_inst)



    #sum the last column of data into a list for the VBR plots    
    # take the averages and then make the VBR
    #PopB_list = [PopB_list[i]/runs for i in range(0,len(PopB_list))]
    #PopV_list = [PopV_list[i]/runs for i in range(0,len(PopV_list))]
    
    
    #
    VBR_list = [PopV_list[i]/PopB_list[i] for i in range(0, len(PopB_list))]
    dfB = dfB.transpose()
    dfV = dfV.transpose()
    location = int(time)
    loc2 = int(0)
    Bstarttotal = dfB[loc2].sum()
    Vstarttotal = dfV[loc2].sum()
    BACinit.append(Bstarttotal)
    VIRinit.append(Vstarttotal)
    BACtotal = dfB[location].sum()
    BACendlist.append(BACtotal)
    VIRtotal = dfV[location].sum()
    VIRendlist.append(VIRtotal)
    TEMPendlist.append(temperature)
    # print(BACendlist[-1])
    print ("run",x,"out of", runs," completed")
    
# Take averages at each time point and plot w/ standard error

PopB_average = PopB_list_internal.mean(axis = 0)
PopB_stdev = PopB_list_internal.std(axis = 0)
PopV_average = PopV_list_internal.mean(axis = 0)
PopV_stdev = PopV_list_internal.std(axis = 0)
VBR_average = VBR_list_internal.mean(axis = 0)
VBR_stdev = VBR_list_internal.std(axis = 0)
Shannon_average = Shannon_list.mean(axis = 0)
Shannon_stdev = Shannon_list.std(axis = 0)
P_average = P_list_internal.mean(axis = 0)
P_stdev = P_list_internal.std(axis = 0)
    


print("The Temperature is: {} and time is: {} \n alpha is: {}, q is: {} \n mu x {}, phi x {}, gamma is: mu/mu+0.1, s is: {}, and P is {} ".format(temperature, time, alpha, Q, mu_x, phi_x,s, P))

print ("Simulation is finished")


# ------- Plotting -------- #

# Population plots; VBR plot
f = plt.figure(figsize=(15, 3), dpi= 600)


ax3 = f.add_subplot(131)
ax4 = f.add_subplot(132)
ax5 = f.add_subplot(133)

ax3.set(yscale="log")
ax4.set(yscale="log")


ax3.title.set_text('Total Bacterial Population')
ax4.title.set_text('Total Virus Population')
ax5.title.set_text('Virus:Bacteria Ratio')   


sns.lineplot(x = timepoints, y = PopB_average, ax = ax3 )
ax3.fill_between(timepoints, PopB_average - PopB_stdev, PopB_average + PopB_stdev, alpha = 0.25)
ax3.set(ylabel='Population mL$^-$$^1$')
#ax3.set_ylim(1e1,1e7)



sns.lineplot(x = timepoints, y = PopV_average, ax = ax4)
ax4.fill_between(timepoints, PopV_average - PopV_stdev, PopV_average + PopV_stdev, alpha = 0.25)
ax4.set(ylabel='Population mL$^-$$^1$')
#ax4.set_ylim(1e2,1e10)

sns.lineplot(x = timepoints, y = PopV_average/PopB_average, ax = ax5 )
#ax5.fill_between(timepoints, VBR_average - VBR_stdev, VBR_average + VBR_stdev, alpha = 0.25)
ax5.set(ylabel='Ratio')
#ax5.set_ylim(0.001,1e7)
ax5.set_yscale('log')



f.tight_layout()



plt.show()


# sea born heat maps of bacteria, virus
#matplotlib.rcParams.update({'font.size': 15})
f2 = plt.figure(figsize=(13, 20), dpi= 600)
ax1 = f2.add_subplot(121)
ax2 = f2.add_subplot(122)

mi, ma = dfB.values.min(), dfB.values.max()
mi_V, ma_V = dfV.values.min(), dfV.values.max()
yticks = np.arange(0,1+tolerance, 0.01).tolist()


ax1.title.set_text('Bacterial Population by Genotype')
ax2.title.set_text('Virus Population by Genotype')


sns.heatmap(dfB, norm=LogNorm(vmin=mi, vmax=ma), yticklabels=['{:,.2f}'.format(y) for y in yticks],cbar_kws={'label': 'Population mL $^-$$^1$'},ax = ax1)
ax1.set(xlabel='time (hrs)')
sns.heatmap(dfV, norm=LogNorm(vmin=mi_V, vmax=ma_V),yticklabels=['{:,.2f}'.format(y) for y in yticks], cbar_kws={'label': 'Population mL $^-$$^1$'}, ax = ax2)
ax2.set(xlabel='time (hrs)')

dfB.to_csv('Bacterial Pop.csv')
dfV.to_csv('Viral Pop.csv')

ax1.invert_yaxis()
ax2.invert_yaxis()
f2.tight_layout()

plt.show()



# I have been saving select values to an excel sheet here:
import xlwt

wb = xlwt.Workbook()
ws = wb.add_sheet('Sheet 1')

first_row = 0

# write each item in the list to consecutive columns on the first row
for i, item in enumerate(P_average):
        ws.write(i, 0, item)
        
for i, item in enumerate(P_stdev):
        ws.write(i, 1, item)

for i, item in enumerate(Shannon_average):
        ws.write(i, 2, item) 
        
for i, item in enumerate(Shannon_stdev):
        ws.write(i, 3, item) 

for i, item in enumerate(VBR_average):
        ws.write(i, 4, item) 
        
for i, item in enumerate(VBR_stdev):
        ws.write(i, 5, item) 
         
for i, item in enumerate(PopB_average):
        ws.write(i, 6, item) 
        
for i, item in enumerate(PopB_stdev):
        ws.write(i, 7, item) 

for i, item in enumerate(PopV_average):
        ws.write(i, 8, item) 
        
for i, item in enumerate(PopV_stdev):
        ws.write(i, 9, item) 
         
# saving here:
wb.save('29July_MSruns10_neg15_s_250.xls')
f.savefig('29July_MultispeciesFigure_1_runs10_neg15_s_250.png', dpi=600)
f2.savefig('29July_MultispeciesFigure_2_runs10_neg15_s_250.png', dpi=600)
