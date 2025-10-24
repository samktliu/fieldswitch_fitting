import numpy as np
from modules import gen_mh,calc_ph
from ph_fit import ph_fit
import matplotlib.pyplot as plt

data = np.load('ph_2.01.npy')
Hz = data[:,0]
ph_p = data[:,1]

# cgs units
D = 65e-7           # cm
Ms = 1495           # emu/cm^3
Hk = 1.81e3         # Oe
t_fl = 1.5e-7       # cm

Aex,wdw,Delta,popt,pcov = ph_fit(Hz,ph_p,True)