import numpy as np
from modules import gen_mh,calc_ph
from ph_fit import ph_fit
import matplotlib.pyplot as plt

# set average Hc
Hc = 2000           # Oe
sig_Hc = 50         # Oe
Hrange = 3500       # Oe
samples = 50

# generate hypothetical ph curve
M_p_arr,M_n_arr,Hz = gen_mh(Hc,sig_Hc,Hrange,samples,True)
# Hz = Hz*79.577    # A/m
ph_p,ph_n = calc_ph(M_p_arr,M_n_arr,Hz,True)

# ph fitting constants

# SI units
# D = 65e-9         # m
# Ms = 1.495e6      # A/m
# Hk = 144035       # A/m
# t_fl = 1.5e-9 

# cgs units
D = 65e-7           # cm
Ms = 1495           # emu/cm^3
Hk = 1.81e3         # Oe
t_fl = 1.5e-7       # cm

ph = np.where(Hz<0,1-ph_n,ph_p)
Aex,wdw,Delta,popt,pcov = ph_fit(Hz,ph_p,True)
Aex,wdw,Delta,popt,pcov = ph_fit(-Hz,1-ph_n,True)