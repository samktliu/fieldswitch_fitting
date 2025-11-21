import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl    
import scienceplots

plt.style.use(['science','nature'])

fig_m,ax_m = plt.subplots(1,1,figsize=(3,2))
fig_dw,ax_dw = plt.subplots(1,1,figsize=(3,2))

legtype = 'even'

folders = ['dwlen_3600','dwlen_3650','dwlen_3675','dwlen_3685','dwlen_3700','dwlen_3800']
fields = []
legend = []
for f in folders:
    legend.append(f.replace('dwlen_','')+' Oe')
    fields.append(int(f.replace('dwlen_','')))
fields = np.array(fields)
cmap = mpl.colormaps['viridis']

if legtype == 'prop':
    frange = np.max(fields)-np.min(fields)
    fmid = 0.5*(np.max(fields)+np.min(fields))
    fcent = fields-fmid
    fscale = 0.9
    f_cmap = (fscale*fcent)+fmid
    colors = cmap((f_cmap-(np.min(fields)))/frange)

if legtype == 'even':
    colors = cmap(np.linspace(0,1,len(fields)))

for idx,f in enumerate(folders):
    t_arr = np.load(f'./{f}/t_arr.npy')
    dwl_arr = np.load(f'./{f}/dwl_arr.npy')
    mzavg_arr = np.load(f'./{f}/mzavg_arr.npy')
    ax_m.plot(t_arr*1e9,mzavg_arr,color=colors[idx])
    ax_dw.plot(t_arr*1e9,dwl_arr*1e9,color=colors[idx])
ax_m.legend(legend)
ax_dw.legend(legend)
ax_m.set_xlabel('Time (ns)')
ax_dw.set_xlabel('Time (ns)')
ax_m.set_ylabel('Unit magnetization \it{m}')
ax_dw.set_ylabel('Domain wall length (nm)')
fig_m.savefig('mzavg_t.svg')
fig_dw.savefig('dwl_t.svg')