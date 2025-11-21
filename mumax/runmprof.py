import subprocess
import shutil 
import os
from tqdm import tqdm
from glob import glob
import numpy as np
import matplotlib.pyplot as plt
from img_proc import dw_length
import scienceplots

plt.style.use(['science','nature'])

run_f = False
copyt_f = False
convert_f = False
process_f = True
plot_f = True

Hfield = 3800 # Oe
filtype = 'canny'
tstep = 1e-10

basedir = f'./dwlen_{Hfield}/'
mmoutputdir = f'./fieldtest_mprof.out/'
npoutputdir = f'{basedir}mprof_npy/'
pngoutputdir = f'{basedir}mprof_png/'
plotdir = f'{basedir}mprof_plots/'
dwoutputdir = f'{basedir}mprof_dw/'

try:
    os.mkdir(basedir)
    print(f"Directory '{basedir}' created successfully.")
except FileExistsError:
    print(f"Directory '{basedir}' already exists.")

try:
    os.mkdir(npoutputdir)
    print(f"Directory '{npoutputdir}' created successfully.")
except FileExistsError:
    print(f"Directory '{npoutputdir}' already exists.")

try:
    os.mkdir(pngoutputdir)
    print(f"Directory '{pngoutputdir}' created successfully.")
except FileExistsError:
    print(f"Directory '{pngoutputdir}' already exists.")

try:
    os.mkdir(dwoutputdir)
    print(f"Directory '{dwoutputdir}' created successfully.")
except FileExistsError:
    print(f"Directory '{dwoutputdir}' already exists.")

try:
    os.mkdir(plotdir)
    print(f"Directory '{plotdir}' created successfully.")
except FileExistsError:
    print(f"Directory '{plotdir}' already exists.")

try:
    os.mkdir(f'{dwoutputdir}{filtype}/')
    print(f"Directory '{f'{dwoutputdir}{filtype}/'}' created successfully.")
except FileExistsError:
    print(f"Directory '{f'{dwoutputdir}{filtype}/'}' already exists.")

if run_f:
    file = open('fieldtest_mprof.mx3','r')
    content = file.readlines()
    content[20] = f'Hmax        := {Hfield} // Oe\n'
    file = open('fieldtest_mprof.mx3','w')
    file.writelines(content)
    file.close()
    subprocess.run(['mumax3','fieldtest_mprof.mx3'], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

if copyt_f:
    shutil.copyfile(f'{mmoutputdir}table.txt',f'{npoutputdir}mprof_table.txt')

if convert_f:
    subprocess.run(['mumax3-convert','-numpy',mmoutputdir+'*.ovf'])
    subprocess.run(['mumax3-convert','-png',mmoutputdir+'*.ovf'])
    for npyfile in tqdm(sorted(glob(mmoutputdir+'*.npy'))):
        fname = os.path.split(npyfile)[1]
        shutil.move(npyfile,npoutputdir+fname)
    for npyfile in tqdm(sorted(glob(mmoutputdir+'*.png'))):
        fname = os.path.split(npyfile)[1]
        shutil.move(npyfile,pngoutputdir+fname)

if process_f:
    npy_paths = sorted(glob(npoutputdir+'*.npy'))
    t_arr = np.arange(0,len(npy_paths)*tstep,tstep)
    mzavg_arr = []
    dwl_arr = []
    for idx,npyfile in enumerate(tqdm(npy_paths)):
        data_np = np.load(npyfile)
        mz = data_np[2,0,:,:]
        mz[mz==0] = np.nan
        mzavg_arr.append(np.nanmean(mz))
        dwl,_ = dw_length(mz,dwoutputdir,idx,filtype)
        dwl_arr.append(dwl)
    mzavg_arr = np.array(mzavg_arr)
    dwl_arr = np.array(dwl_arr)
    print(f'Domain wall length array (nm):\n{dwl_arr*1e9}')
    np.save(f'{basedir}t_arr',t_arr)
    np.save(f'{basedir}mzavg_arr',mzavg_arr)
    np.save(f'{basedir}dwl_arr',dwl_arr)

if plot_f:
    fig,ax = plt.subplots(1,1,figsize=(3,2))
    fig2,ax2 = plt.subplots(1,1,figsize=(3,2))
    ax.plot(t_arr*1e9,mzavg_arr)
    ax2.plot(t_arr*1e9,dwl_arr*1e9)
    fig.savefig(f'{plotdir}mzavg.svg')
    fig2.savefig(f'{plotdir}dwl_{filtype}.svg')