import matplotlib.pyplot as plt
import numpy as np
import glob
import scienceplots

plt.style.use(['science','nature'])

folders = ['./A1.600/','./A1.776/','./A1.900/']
# folders = ['./A1.776/']

for f in folders:
    mzs = []
    for fname in glob.iglob(f+'*.txt'):
        data = np.genfromtxt(fname)
        mz = data[:,3]
        mzs.append(mz)
        Bz = data[:,-1]*1e4
    mzs_arr = np.array(mzs).T
    p = []
    for i in range(len(Bz)):
        p.append(np.sum(mzs_arr[i,:]>0)/len(mzs))
    p = np.array(p)
    p = np.concatenate([Bz,p]).reshape(2,len(Bz)).T
    np.save(f'pmm_{f[2:-1]}',p)
    plt.plot(p[:,0],p[:,1],'^-')

plt.savefig('test.svg')