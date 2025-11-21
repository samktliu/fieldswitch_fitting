import numpy as np
from scipy import ndimage
import matplotlib.pyplot as plt
from skimage import feature
import os

cellsize = 2.03125e-9

def dw_length(mz_arr,folder,ind=0,filter='canny'):

    mz_dwprof = np.sign(mz_arr)
    # mz_dwprof = mz_arr
    if filter == 'sobel':
        sobel_h = ndimage.sobel(mz_dwprof, 0)
        sobel_v = ndimage.sobel(mz_dwprof, 1)
        magnitude = np.sqrt(sobel_h**2 + sobel_v**2)
    if filter == 'canny':
        magnitude = feature.canny(mz_dwprof,sigma=0.05)
    magnitude = np.float64(magnitude)
    magnitude[np.isnan(mz_arr)] = np.nan
    # fig, axs = plt.subplots(2, 2, figsize=(8, 8))
    # plt.gray()
    # axs[0, 0].imshow(mz_arr)
    # axs[0, 1].imshow(sobel_h)
    # axs[1, 0].imshow(sobel_v)
    # axs[1, 1].imshow(magnitude)
    # fig.savefig('test_sobel.svg')
    plt.imshow(magnitude,origin='lower')
    plt.savefig(f'{folder}{filter}/{ind}.png')
    if filter == 'sobel':
        magnitude = magnitude/np.nanmax(magnitude)
    dwl = np.nansum(magnitude)*cellsize
    return dwl,mz_dwprof