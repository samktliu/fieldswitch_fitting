import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl    
import scienceplots

plt.style.use(['science','nature'])

fig_m,ax_m = plt.subplots(1,1,figsize=(3,2))
fig_dw,ax_dw = plt.subplots(1,1,figsize=(3,2))
fig_nrg,ax_nrg = plt.subplots(1,1,figsize=(3,2))

folder = 'dwlen_4100'
field = int(folder.replace('dwlen_',''))

tstep = 1e-10
means = []
stds = []
times = []

cmap = mpl.colormaps['rainbow']
colors = cmap(np.linspace(0,1,10))

dwl_arrs_t = []
maxlen = 0
for seed in range(10):
    dwl_arr = np.load(f'./{folder}/s{seed}/dwl_arr.npy')
    mz = np.load(f'./{folder}/s{seed}/mzavg_arr.npy')
    t_arr = np.load(f'./{folder}/s{seed}/t_arr.npy')
    table = np.loadtxt(f'./{folder}/s{seed}/mprof_npy/mprof_table.txt')
    t_table = table[:,0]
    nrg_table = table[:,7]
    ax_m.plot(t_arr*1e9,mz,color=colors[seed])
    ax_dw.plot(t_arr*1e9,dwl_arr*1e9,color=colors[seed])
    ax_nrg.plot(t_table*1e9,nrg_table,color=colors[seed])

ax_m.set_xlabel('Time (ns)')
ax_m.set_ylabel('Unit magnetization $\it{m_z}$')
fig_m.savefig('mzavg_raw.svg')
ax_dw.set_xlabel('Time (ns)')
ax_dw.set_ylabel('Domain wall length (nm)')
fig_dw.savefig('dwl_raw.svg')
ax_nrg.set_xlabel('Time (ns)')
ax_nrg.set_ylabel('Total energy (J)')
fig_nrg.savefig('nrg_raw.svg')