import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl    
import scienceplots

plt.style.use(['science','nature'])

# fig_m,ax_m = plt.subplots(1,1,figsize=(3,2))
fig_dw,ax_dw = plt.subplots(1,1,figsize=(3,2))

legtype = 'even'

folders = ['dwlen_3800','dwlen_3900','dwlen_4000','dwlen_4100']
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

tstep = 1e-10
means = []
stds = []
times = []
for idx,f in enumerate(folders):
    dwl_arrs_t = []
    maxlen = 0
    for seed in range(10):
        dwl_arr = np.load(f'./{f}/s{seed}/dwl_arr.npy')
        if np.sum(np.abs(dwl_arr)) == 0:
            print(f'series {f}, seed {seed} no switch')
        elif dwl_arr[-1] != 0:
            print(f'series {f}, seed {seed} partial switch')
        else:
            dwl_trimmed = np.trim_zeros(dwl_arr)
            dwl_arrs_t.append(dwl_trimmed)
            if len(dwl_trimmed) > maxlen:
                maxlen = len(dwl_trimmed)
    dwl_arrs = []
    for dwl in dwl_arrs_t:
        temp_dwl = np.pad(dwl,(0,maxlen-len(dwl)),'constant',constant_values=(0,0))
        dwl_arrs.append(temp_dwl)
    dwl_arrs = np.array(dwl_arrs).T
    print(dwl_arrs.shape)
    time = np.linspace(0,tstep*(dwl_arrs.shape[0]-1),dwl_arrs.shape[0])
    times.append(time)
    means.append(np.mean(dwl_arrs,axis=1))
    stds.append(np.std(dwl_arrs,axis=1))

for idx,f in enumerate(folders):
#     ax_m.plot(t_arr*1e9,mzavg_arr,color=colors[idx])
    ax_dw.plot(times[idx]*1e9,means[idx]*1e9,color=colors[idx])
for idx,f in enumerate(folders):
    ax_dw.fill_between(times[idx]*1e9,(means[idx]-stds[idx])*1e9,(means[idx]+stds[idx])*1e9,color=colors[idx],alpha=0.3,linewidth=0)
    print(f'avg std for {f}: {np.mean(stds[idx])*1e9}')
# ax_m.legend(legend)
ax_dw.legend(legend)
# ax_m.set_xlabel('Time (ns)')
ax_dw.set_xlabel('Time (ns)')
# ax_m.set_ylabel('Unit magnetization $\it{m_z}$')
ax_dw.set_ylabel('Domain wall length (nm)')
# fig_m.savefig('mzavg_t.svg')
fig_dw.savefig('dwl_t.svg')