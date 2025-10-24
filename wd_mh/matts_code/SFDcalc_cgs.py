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
import scienceplots

plt.style.use(['science','nature'])

#### define constants #####
kb = 1.380649e-16 # erg/K
Temp = 300 # K
kbT = kb*Temp # erg
fatt = 1E9 # Hz
mu0 = 4*np.pi*1E-7

#device params
thickness = 1.5E-7 # cm
diam = 65E-7 # cm

# time dependence
samprate = int(2E3)
deltat = 1/samprate  # assume time from sampling rate of 1/yyy kHz 
print(f'deltat = {deltat:.4f}')
# H ranges, in T
Hmin = 1 
Hmax = 3500
Herr = 5
fields = np.linspace(Hmin, Hmax, num=samprate)
# number of steps set by sample rate for now

#number of switching trials
numTrials = 1000
 
# parameters for energy barrier. Keep them constant for a while
## domain wall stuff
edw = 6.6 # erg/cm^2
wdw = 12.7E-7  # cm
## material stuff
Ms = 1495 # emu/cm^3
Hk = 1.81e3
print(f'Ms = {Ms} emu/cm^3')
Aexeff = edw*edw/(8*Ms*Hk)
print(f'Aex = {Aexeff}')
print(f'wdw = {wdw*1e7} nm')

# define Eb dw class. Inputs are in SI
EbC = EBdwclass.EbclassHCGS(edw,wdw,thickness,diam,Ms)
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
            # define switching probability function (eq 10 in Metropolis-using paper:  Breth, L., Suess, D., Vogler, C., Bergmair, B., Fuger, M., Heer, R., & Brueckl, H. (2012). Thermal switching field distribution of a single domain particle for field-dependent attempt frequency. Journal of Applied Physics, 112(2). https://doi.org/10.1063/1.4737413)
        ProbE = 1-np.exp(-deltat*fatt*np.exp(-Eb/kbT))
        ProbR = random.uniform(0,1)
        if ProbE>ProbR:
            # print(f'Eb = {Eb}; prob E= {ProbE}; prob R = {ProbR}; H = {theField*mu0}; index = {index}')
            # increment count on switching at H = theField
            # theDist[index]+=1
            theDist[index:]+=1
            Switch = True
            switches = switches+1
        else:
            index = index+1
# need to set P = 1 for H > Hk or nonphysical?

normDist = theDist/numTrials
SFD = np.gradient(normDist)

fig, ax = plt.subplots(1,1,figsize=(3,2))
ax.plot(fields, theDist,'^',color='tab:purple',mfc='none')
fig.savefig('cgs_fieldsxmu0_theDist.svg')
fig, ax = plt.subplots(1,1,figsize=(3,2))
ax.plot(fields, normDist,'^',color='tab:red',mfc='none')
# ax.set_title(f'SFD: dt={deltat:.4f}, Ms={Ms:.2f}, dia={diam}, #Trials={numTrials}')
#ax.text(0.1,0.1, f'dwmod={dwmod:.2f}', transform=ax.transAxes)
#ax.text(0.7,0.8, f'edw(1-dwmod)', transform=ax.transAxes)
fig.savefig('cgs_fieldsxmu0_normDist.svg')
fig, ax = plt.subplots(1,1,figsize=(3,2))
ax.plot(fields, SFD,'o',color='tab:green',mfc='none')
fig.savefig('cgs_fieldsxmu0_sfd.svg')
arr = np.stack((fields,normDist)).T
print(np.shape(arr))
# fname = f'sfd2-Ms{Ms:.2f}-dia{diam*1E9:.0f}-dt{deltat:.2e}-edw{edw:.3E}-wdw{wdw:.3E}-Herr{Herr:.0f}Oe-{numTrials}trials.txt'
fname = f'ph.npy'
# print(fname)
# pathstub = r'\\687spintronics1\X\Shared Data\Data\CHIPS exchange\Sims\Metropolis-Ebdw'
# outfile = os.path.join(pathstub, fname)
np.save(fname, arr)
