#!/usr/bin/env python3
import os
import sys
import netCDF4
import numpy as np
import matplotlib.pyplot as plt

def quasilinear_estimate(gs2File, fractionToConsider = 0.3, show=False, savefig=False, out_folder=''):
    # Read GS2 output file
    file2read = netCDF4.Dataset(gs2File,'r')
    # Index for the start of the fitting domain (fraction of time to consider)
    time = file2read.variables['t'][()]
    startIndex  = int(len(time)*(1-fractionToConsider))
    # find the growth rate gamma for each ky by fitting |phi|^2 to an exponential(2 gamma)
    phi2_by_ky  = np.array(file2read.variables['phi2_by_ky'][()])
    ky = file2read.variables['ky'][()]
    assert phi2_by_ky.shape[0]==len(time)
    assert phi2_by_ky.shape[1]==len(ky)
    growthRates  = []
    for i in range(len(ky)):
        mask  = np.isfinite(phi2_by_ky[:,i])
        data_x = time[mask]
        data_y = phi2_by_ky[mask,i]
        fitX  = np.polyfit(data_x[startIndex:], np.log(data_y[startIndex:]), 1)
        thisGrowthRate  = fitX[0]/2
        growthRates.append(thisGrowthRate)
    # jacobian = nabla x \times nabla y . nabla z
    jacob = np.array(file2read.variables['jacob'][()])
    # gds2 = |\nabla y|^2
    gds2 = np.array(file2read.variables['gds2'][()])
    # Electrostatic potential phi eigenfunction
    phi = np.array(file2read.variables['phi'][()]) # indices ky, kx, z, real/imaginary
    # Assuming that there is only one kx, calculate |phi|^2
    phi2_by_ky_of_z = phi[:,0,:,0]**2+phi[:,0,:,1]**2
    z = file2read.variables['theta'][()]
    assert phi2_by_ky_of_z.shape[0]==len(ky)
    assert phi2_by_ky_of_z.shape[1]==len(z)
    # kperp^2 = ky^2 * |\nabla y|^2
    # kperp2 = np.array([ky_each**2 * gds2 for ky_each in ky])
    # weighted_kperp2 = integral( kperp2 * |phi|^2 * jacob )/integral(|phi|^2 * jacob )
    weighted_kperp2 = np.array([np.sum(ky_each*ky_each*gds2*phi2_by_ky_of_z[i]*jacob)/np.sum(phi2_by_ky_of_z[i]*jacob) for i, ky_each in enumerate(ky)])
    # weighted_gamma = gamma/weighted_kperp2
    weighted_gamma = growthRates/weighted_kperp2
    if show:
        fig, ax = plt.subplots(nrows=1, ncols=2)
        fig.set_size_inches(8.5, 5)
        # Plot growth rates
        ax[0].plot(ky, growthRates, label='$\gamma$')
        ax[0].plot(ky, growthRates/ky, label='$\gamma/k_y$')
        ax[0].plot(ky, growthRates/ky/ky, label='$\gamma/k_y^2$')
        ax[0].plot(ky, weighted_gamma, label=r'$\gamma/\langle k_{\perp}^2 \rangle$')
        # Plot eigenfunctions
        max_index = np.nanargmax(growthRates)
        indices = [0,-1,max_index]
        labels = ['lowest $k_y$', 'highest $k_y$', 'max($\gamma$)']
        for index, label in zip(indices, labels):
            phi2_by_ky_of_z0 = phi2_by_ky_of_z[index, int(len(z)/2)]
            ax[1].plot(z, phi2_by_ky_of_z[index,:]/phi2_by_ky_of_z0, label=rf'$|\phi|^2/|\phi_0|^2$ for {label}')
        ## Plot parameters
        ax[0].set_xlabel('$k_y$', fontsize=14)
        ax[0].set_xscale('log')
        ax[0].legend(fontsize=12)
        ax[1].set_xlabel('$z$', fontsize=14)
        ax[1].legend(fontsize=10)
        plt.tight_layout()
        if savefig: 
            plt.savefig(os.path.join(out_folder,gs2File.split("/")[-1]+'_quasilinear.png'))
        else:
            plt.show()
        plt.close()
    return weighted_gamma

if __name__ == "__main__":
    weighted_gamma = quasilinear_estimate(sys.argv[1], show=True)
    print('Resulting gamma/<k_perp^2>:')
    print(weighted_gamma)