# -*- coding: utf-8 -*-
"""
Created on Sun Jun 22 16:23:21 2025
Function to test Monte Carlo / metropolis method for determining SFD

@author: pufall 
"""


import numpy as np
import matplotlib.pyplot as plt
import math
import EBdwclass
import random
import os
import matplotlib.pyplot as plt
from tqdm import tqdm

#### define constants #####
kb = 1.380649E-23  #J/K
Temp = 300
kbT = kb*Temp
fatt = 1E9
mu0 = 4*np.pi*1E-7

# time dependence
samprate = int(2E3)
deltat = 1/samprate  # assume time from sampling rate of 1/yyy kHz 
print(f'deltat = {deltat:.5f}')
# H ranges, in T
Hmin = 0.0001 
Hmax = 0.5
Herr = 0.0005/mu0   # assume 5 Oe max uncertainty on field 
fields = np.linspace(Hmin, Hmax, num=samprate)/mu0 # convert to A/m
 # number of steps set by sample rate for now

#number of switching trials
numTrials = 1000
 
# parameters for energy barrier. Keep them constant for a while
## domain wall stuff
dwmod = 0.0
scale = 1.0
edw = scale*(1-dwmod)*5.5E-3 # in J/m2. From Fig 4
wdw = scale*(12E-9)/(1-dwmod)  # from fig 4
## material stuff
Ms = 1.0*1200*(4*np.pi*1E-4)/mu0  # emu/cc to T to A/m
print(f'Ms = {Ms}A/m => {mu0*Ms:.3f} T')
Aexeff = edw*wdw/(8*np.log(2))

# now calc inputs edw and wdw from film level quantities
# Aex = 10E-12  # pJ/m
# Meff = -0.3/mu0 # tesla/mu0, negative is out of plane, from FMR
# Hk = Ms - Meff
# Keff = mu0*Ms*Hk/2
# edwc = 4*np.sqrt(Keff*Aex)
# wdwc = 2*np.log(2)*np.sqrt(Aex/Keff)
# wdwc2 = 8*np.log(2)*Aex/edwc 
# print(f'energy of dom = {edwc:.5f} J/m2; width = {wdwc:.5g} m, Alt. width = {wdwc2:.5g} m')

#device params
thickness = 1.5E-9 # thickness 
diam = 65E-9 

# define Eb dw class. Inputs are in SI
EbC = EBdwclass.EbclassH(edw,wdw,thickness,diam,Ms)
#EbC.setH(0.001)
#Ebtest = EbC.Eb(200000)
#print(Ebtest)
random.seed(666) # set seed for initial diags
theDist = np.zeros(samprate)
switches = 0
for i in tqdm(range(numTrials)):
    index = 0        
    Switch = False
    while Switch is False:
        theField = fields[index]+ Herr*2*(random.random()-0.5)  #include random field noise here. could put read error in elsewhere
        Eb = EbC.Eb(theField)
        print(Eb/kbT)
            # define switching probability function (eq 10 in Metropolis-using paper:  Breth, L., Suess, D., Vogler, C., Bergmair, B., Fuger, M., Heer, R., & Brueckl, H. (2012). Thermal switching field distribution of a single domain particle for field-dependent attempt frequency. Journal of Applied Physics, 112(2). https://doi.org/10.1063/1.4737413)
        ProbE = 1-np.exp(-deltat*fatt*np.exp(-Eb/kbT))
        ProbR = random.uniform(0,1)
        if ProbE> ProbR:
            #print(f'Eb = {Eb}; prob E= {ProbE}; prob R = {ProbR}; H = {theField*mu0}; index = {index}')
            # increment count on switching at H = theField
#            theDist[index]+=1
            theDist[index:]+=1
            Switch = True
            switches = switches +1
        else:
            index = index +1
# need to set P = 1 for H > Hk or nonphysical?

normDist = theDist/numTrials
SFD = np.gradient(normDist)

fig, ax = plt.subplots()
ax.plot(fields*mu0, theDist, 'ok')
fig.savefig('fieldsxmu0_theDist.svg')
fig, ax = plt.subplots()
ax.plot(fields*mu0, normDist, 'ok')
ax.set_title(f'SFD: dt={deltat:.4f}, Ms={Ms*mu0:.2f}, dia={diam}, #Trials={numTrials}')
#ax.text(0.1,0.1, f'dwmod={dwmod:.2f}', transform=ax.transAxes)
#ax.text(0.7,0.8, f'edw(1-dwmod)', transform=ax.transAxes)
ax.text(0.7,0.75, f'wdw, edw scale = {scale}', transform=ax.transAxes)
fig.savefig('fieldsxmu0_normDist.svg')
fig, ax = plt.subplots()
ax.plot(fields*mu0, SFD, 'ok')
fig.savefig('fieldsxmu0_sfd.svg')
arr = np.stack((fields*mu0,normDist)).T
print(np.shape(arr))
fname = f'sfd2-Ms{Ms*mu0:.2f}-dia{diam*1E9:.0f}-dt{deltat:.2e}-edw{edw*scale*(1-dwmod):.3E}-wdw{wdw*scale/(1-dwmod):.3E}-Herr{Herr*mu0*10000:.0f}Oe-{numTrials}trials.txt'
# print(fname)
# pathstub = r'\\687spintronics1\X\Shared Data\Data\CHIPS exchange\Sims\Metropolis-Ebdw'
# outfile = os.path.join(pathstub, fname)
np.savetxt(fname, arr)
