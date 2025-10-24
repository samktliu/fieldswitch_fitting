import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import scienceplots

plt.style.use(['science','nature'])

#### define constants #####
kb = 1.380649e-16 # erg/K
Temp = 300
kbT = kb*Temp
fatt = 1E9
mu0 = 4*np.pi*1E-7
samprate = int(2E3)
deltat = 1/samprate 
D = 65e-7 # cm
Ms = 1495 # emu/cm^3
Hk = 1.81e3 # Oe
t_fl = 1.5e-7 # cm

def ph_func(h,Aex,delta):
    # edw calculation
    edw = np.sqrt(8*Ms*Hk*Aex)
    # q1 calculation
    eps = edw/(Ms*np.abs(h)*D)
    q1 = 1 + eps - np.sqrt(1+eps*eps)
    q1_p = q1 + delta
    q1_n = q1 - delta
    # thetaq calculation
    thetaq1 = np.atan(q1*(1-q1/2)/(1-q1))
    thetaq1_p = np.atan(q1_p*(1-q1_p/2)/(1-q1_p))
    thetaq1_n = np.atan(q1_n*(1-q1_n/2)/(1-q1_n))

    thetaq1 = np.where(q1>1,np.pi-np.atan(q1*(1-q1/2)/(q1-1)),np.atan(q1*(1-q1/2)/(1-q1)))
    thetaq1_p = np.where(q1_p>1,np.pi-np.atan(q1_p*(1-q1_p/2)/(q1_p-1)),np.atan(q1_p*(1-q1_p/2)/(1-q1_p)))
    thetaq1_n = np.where(q1_n>1,np.pi-np.atan(q1_n*(1-q1_n/2)/(q1_n-1)),np.atan(q1_n*(1-q1_n/2)/(1-q1_n)))
    # Ad calculation
    # Adq1 = 0.25*D*2/(thetaq1-np.tan(thetaq1)+(np.pi/2-thetaq1)*np.tan(thetaq1)*np.tan(thetaq1))
    Adq1_p = 0.25*D*D*(thetaq1_p-np.tan(thetaq1_p)+(np.pi/2-thetaq1_p)*np.tan(thetaq1_p)*np.tan(thetaq1_p))
    Adq1_n = 0.25*D*D*(thetaq1_n-np.tan(thetaq1_n)+(np.pi/2-thetaq1_n)*np.tan(thetaq1_n)*np.tan(thetaq1_n))
    # Ldw calculation
    Ldq1 = D*(np.pi/2-thetaq1)*np.tan(thetaq1)
    # Eq calculation
    Eq1 = Ldq1*edw*t_fl+np.abs(h)*Ms*t_fl*(((D*D*np.pi)/4)-Adq1_p-Adq1_n)
    Eb = Eq1-((np.pi/4)*D*D*Ms*h*t_fl)
    # Calculate P(H)
    # print(-(Eb/kbT))
    ph = 1-np.exp(-fatt*deltat*np.exp(-(Eb/kbT)))
    # fig,ax = plt.subplots(1,1,figsize=(3,2))
    # ax.plot(h,Adq1_p,color='black')
    # fig.savefig('temp.svg')
    return ph

def ph_fit(h,prob,save_fig):
    popt,pcov = curve_fit(ph_func,h,prob,p0=[1e-6,0.3])
    # popt,pcov = curve_fit(ph_func,h,prob,p0=[1e-6,0.3],bounds=(0, [100e-6, 1]))
    Aex = popt[0]
    delta = popt[1]
    wdw = np.abs(D*delta)
    edw = np.sqrt(8*Ms*Hk*Aex)
    Delta = D*t_fl*edw/kbT
    print(popt)
    print(pcov)
    print(f'Aex = {Aex} | w_dw = {wdw*1e7} nm | Delta = {Delta}')
    if save_fig:
        fig,ax = plt.subplots(1,1,figsize=(3,2))
        ax.plot(h,prob,'^',color='black',mfc='none')
        ax.plot(h,ph_func(h,popt[0],popt[1]),color='tab:red')
        fig.savefig('fitting.svg')
    return Aex,wdw,Delta,popt,pcov
