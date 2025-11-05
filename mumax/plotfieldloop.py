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
        Bz = data[:,-1]
    mzs_arr = np.array(mzs).T
    p = []
    for i in range(len(Bz)):
        p.append(np.sum(mzs_arr[i,:]>0)/len(mzs))
    p = np.array(p)
    np.save(f'p_{f[2:-1]}',p)
    plt.plot(Bz,p,'^-')
np.save(f'Bz.npy',Bz)

plt.savefig('test.svg')