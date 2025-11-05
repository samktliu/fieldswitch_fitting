import numpy as np
from modules import gen_mh,calc_ph
from ph_fit_mmfatt import ph_fit
import matplotlib.pyplot as plt

data = np.load('pmm_A1.776.npy')
Hz = data[:,0]
# hfit = Hz
hfit = np.linspace(0,np.max(Hz),1000)
ph_p = data[:,1]

# cgs units
D = 65e-7           # cm
Ms = 1495           # emu/cm^3
Hk = 1.81e3         # Oe
t_fl = 1.5e-7       # cm

Aex,wdw,Delta,fatt,popt,pcov = ph_fit(Hz,hfit,ph_p,True)