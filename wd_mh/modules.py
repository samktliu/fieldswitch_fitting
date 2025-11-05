import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import scienceplots

plt.style.use(['science','nature'])

def gen_mh(Hc=2000,sig_Hc=50,Hrange=3500,samples=100,save_fig=True):

    Hc_arr_p = np.random.normal(Hc,sig_Hc,samples)
    Hc_arr_n = np.random.normal(-Hc,sig_Hc,samples)
    Hz = np.linspace(-Hrange,Hrange,10000)

    M_p_arr = []
    M_n_arr = []
    for i in tqdm(range(samples)):
        M_p = np.zeros_like(Hz)
        M_p[Hz>Hc_arr_p[i]] = 1
        M_p_arr.append(M_p)
        M_n = np.zeros_like(Hz)
        M_n[Hz>Hc_arr_n[i]] = 1
        M_n_arr.append(M_n)

    M_p_arr = np.array(M_p_arr).T
    M_n_arr = np.array(M_n_arr).T
    if save_fig:
        fig,ax = plt.subplots(1,1,figsize=(3,2))
        ax.plot(Hz,M_p_arr,color='tab:blue',alpha=0.4)
        ax.plot(Hz,M_n_arr,color='tab:red',alpha=0.4)
        fig.savefig('mh.svg')
    return M_p_arr,M_n_arr,Hz

def calc_ph(M_p_arr,M_n_arr,Hz,save_fig):
    M_p_norm = norm_mh(M_p_arr)
    M_n_norm = norm_mh(M_n_arr)
    ph_p = np.mean(M_p_norm,axis=1)
    ph_n = np.mean(M_n_norm,axis=1)
    if save_fig:
        fig,ax = plt.subplots(1,1,figsize=(3,2))
        ax.plot(Hz,ph_p,color='tab:blue')
        ax.plot(Hz,ph_n,color='tab:red')
        fig.savefig('ph.svg')
    return ph_p,ph_n

def norm_mh(m):
    temp = m-np.min(m)
    return temp/np.max(temp)