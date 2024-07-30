import glob
import netCDF4
import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
from tempfile import mkstemp
from os import fdopen, remove
from shutil import move, copymode

from simsopt.mhd import vmec_fieldlines
def to_gs2(filename, vs, s=0.5, alpha=0, nlambda=30, theta1d=None, phi1d=None, phi_center=0):
    r"""
    Compute field lines and geometric quantities along the
    field lines in a vmec configuration needed to run the
    gyrokinetic GS2 code.
    Args:
        vs: Either an instance of :obj:`simsopt.mhd.vmec.Vmec`
          or the structure returned by :func:`vmec_splines`.
        s: Values of normalized toroidal flux on which to construct the field lines.
          You can give a single number, or a list or numpy array.
        alpha: Values of the field line label :math:`\alpha` on which to construct the field lines.
          You can give a single number, or a list or numpy array.
        theta1d: 1D array of :math:`\theta_{pest}` values, setting the grid points
          along the field line and the parallel extent of the field line.
        phi1d: 1D array of :math:`\phi` values, setting the grid points along the
          field line and the parallel extent of the field line.
        phi_center: :math:`\phi_{center}`, an optional shift to the toroidal angle
          in the definition of :math:`\alpha`.
    """
    try: assert not isinstance(alpha, (list, tuple, np.ndarray))
    except Exception: raise ValueError("Only working for a single field line, alpha should be a scalar quantity")
    arrays = vmec_fieldlines(vs, s, alpha, theta1d=theta1d, phi1d=phi1d, phi_center=phi_center, plot=False, show=True)
    nperiod = 1
    drhodpsi = 1.0
    rmaj     = 1.0
    kxfac    = 1.0
    shat = arrays.shat[0]
    q = 1/arrays.iota[0]
    phi = arrays.phi[0,0]
    ntheta=len(phi)-1
    ntgrid=int(np.floor(len(phi)/2))
    bMax=np.max(arrays.bmag[0])
    bMin=np.min(arrays.bmag[0])
    with open(filename,'w') as f:
        f.write(f"nlambda\n{nlambda}\nlambda")
        for i in range(nlambda):
            f.write(f"\n{(bMax - bMax*i + bMin*i - bMin*nlambda)/(bMax*bMin - bMax*bMin*nlambda)}")
        f.write("\nntgrid nperiod ntheta drhodpsi rmaj shat kxfac q")
        f.write(f"\n{ntgrid} {nperiod} {ntheta} {drhodpsi} {rmaj} {shat} {kxfac} {q}")
        f.write("\ngbdrift gradpar grho tgrid")
        for gbdrift, gradpar, tgrid in zip(arrays.gbdrift[0,0], arrays.gradpar_phi[0,0], phi):
            f.write(f"\n{gbdrift} {gradpar} 1.0 {tgrid}")
        f.write("\ncvdrift gds2 bmag tgrid")
        for cvdrift, gds2, bmag, tgrid in zip(arrays.cvdrift[0,0], arrays.gds2[0,0], arrays.bmag[0,0], phi):
            f.write(f"\n{cvdrift} {gds2} {bmag} {tgrid}")
        f.write("\ngds21 gds22 tgrid")
        for gds21, gds22, tgrid in zip(arrays.gds21[0,0], arrays.gds22[0,0], phi):
            f.write(f"\n{gds21} {gds22} {tgrid}")
        f.write("\ncvdrift0 gbdrift0 tgrid")
        for cvdrift0, gbdrift0, tgrid in zip(arrays.cvdrift0[0,0], arrays.gbdrift0[0,0], phi):
            f.write(f"\n{cvdrift0} {gbdrift0} {tgrid}")
        f.write("\nRplot Rprime tgrid")
        for tgrid in phi:
            f.write(f"\n0.0 0.0 {tgrid}")
        f.write("\nZplot Rprime tgrid")
        for tgrid in phi:
            f.write(f"\n0.0 0.0 {tgrid}")
        f.write("\naplot Rprime tgrid")
        for tgrid in phi:
            f.write(f"\n0.0 0.0 {tgrid}")

# Get growth rates
def getgamma(stellFile, fractionToConsider=0.3, savefig=False):
    f = netCDF4.Dataset(stellFile,'r',mmap=False)
    phi2 = np.log(f.variables['phi2'][()])
    t = f.variables['t'][()]
    startIndex = int(len(t)*(1-fractionToConsider))
    mask = np.isfinite(phi2)
    data_x = t[mask]
    data_y = phi2[mask]
    fit = np.polyfit(data_x[startIndex:], data_y[startIndex:], 1)
    poly = np.poly1d(fit)
    GrowthRate = fit[0]/2
    omega_average_array = np.array(f.variables['omega_average'][()])
    omega_average_array_omega = np.mean(omega_average_array[startIndex:,:,0,0],axis=0)
    omega_average_array_gamma = np.mean(omega_average_array[startIndex:,:,0,1],axis=0)
    max_index = np.nanargmax(omega_average_array_gamma)
    gamma = omega_average_array_gamma[max_index]
    omega = omega_average_array_omega[max_index]
    kyX  = f.variables['ky'][()]
    ky_max = kyX[max_index]
    # gamma  = np.mean(f.variables['omega'][()][startIndex:,0,0,1])
    # omega  = np.mean(f.variables['omega'][()][startIndex:,0,0,0])
    #fitRes = np.poly1d(coeffs)
    # if not os.path.exists(stellFile+'_phi2.pdf'):
    if savefig:
        fig = plt.figure(figsize=(7.5,4.0))
        ##############
        plt.plot(t, phi2,'.', label=r'data - $\gamma_{GS2} = $'+str(gamma))
        plt.plot(t, poly(t),'-', label=r'fit - $\gamma = $'+str(GrowthRate))
        ##############
        plt.legend(loc=0,fontsize=14)
        plt.xlabel(r'$t$');plt.ylabel(r'$\ln |\hat \phi|^2$')
        plt.subplots_adjust(left=0.16, bottom=0.19, right=0.98, top=0.97)
        plt.savefig(stellFile+'_phi2.png')
        plt.close()
    return GrowthRate, omega, ky_max

# Save final eigenfunction
def eigenPlot(stellFile):
    f = netCDF4.Dataset(stellFile,'r',mmap=False)
    y = f.variables['phi'][()]
    x = f.variables['theta'][()]
    plt.figure()
    omega_average_array = np.array(f.variables['omega_average'][()])
    fractionToConsider=0.3
    tX   = f.variables['t'][()]
    startIndexX  = int(len(tX)*(1-fractionToConsider))
    omega_average_array_gamma = np.mean(omega_average_array[startIndexX:,:,0,1],axis=0)
    max_index = np.nanargmax(omega_average_array_gamma)
    phiR0= y[max_index,0,int((len(x))/2),0]
    phiI0= y[max_index,0,int((len(x))/2),1]
    phi02= phiR0**2+phiI0**2
    phiR = (y[max_index,0,:,0]*phiR0+y[max_index,0,:,1]*phiI0)/phi02
    phiI = (y[max_index,0,:,1]*phiR0-y[max_index,0,:,0]*phiI0)/phi02
    ##############
    plt.plot(x, phiR, label=r'Re($\hat \phi/\hat \phi_0$)')
    plt.plot(x, phiI, label=r'Im($\hat \phi/\hat \phi_0$)')
    ##############
    plt.xlabel(r'$\theta$');plt.ylabel(r'$\hat \phi$')
    plt.legend(loc="upper right")
    plt.subplots_adjust(left=0.16, bottom=0.19, right=0.98, top=0.93)
    plt.savefig(stellFile+'_eigenphi.png')
    plt.close()
    return 0
##### Function to obtain gamma and omega for each ky
def gammabyky(stellFile,fractionToConsider=0.6, savefig=False):
    # Compute growth rate:
    fX   = netCDF4.Dataset(stellFile,'r',mmap=False)
    tX   = fX.variables['t'][()]
    kyX  = fX.variables['ky'][()]
    phi2_by_kyX  = fX.variables['phi2_by_ky'][()]
    omegaX  = fX.variables['omega'][()]
    startIndexX  = int(len(tX)*(1-fractionToConsider))
    growthRateX  = []
    ## assume that kyX=kyNA
    for i in range(len(kyX)):
        maskX  = np.isfinite(phi2_by_kyX[:,i])
        data_xX = tX[maskX]
        data_yX = phi2_by_kyX[maskX,i]
        fitX  = np.polyfit(data_xX[startIndexX:], np.log(data_yX[startIndexX:]), 1)
        thisGrowthRateX  = fitX[0]/2
        growthRateX.append(thisGrowthRateX)
    # Compute real frequency:
    realFreqVsTimeX  = []
    realFrequencyX   = []
    for i in range(len(kyX)):
        realFreqVsTimeX.append(omegaX[:,i,0,0])
        realFrequencyX.append(np.mean(realFreqVsTimeX[i][startIndexX:]))
    if savefig:
        fig = plt.figure()
        
        numRows = 1
        numCols = 2

        plt.subplot(numRows, numCols, 1)
        plt.plot(kyX,growthRateX,'.-')
        plt.xlabel(r'$k_y$')
        plt.ylabel(r'$\gamma$')
        plt.xscale('log')
        plt.rc('font', size=8)
        plt.rc('axes', labelsize=8)
        plt.rc('xtick', labelsize=8)
        # plt.legend(frameon=False,prop=dict(size='xx-small'),loc=0)

        plt.subplot(numRows, numCols, 2)
        plt.plot(kyX,realFrequencyX,'.-')
        plt.xlabel(r'$k_y$')
        plt.ylabel(r'$\omega$')
        plt.xscale('log')
        plt.rc('font', size=8)
        plt.rc('axes', labelsize=8)
        plt.rc('xtick', labelsize=8)
        # plt.legend(frameon=False,prop=dict(size=12),loc=0)

        plt.tight_layout()
        #plt.subplots_adjust(left=0.14, bottom=0.15, right=0.98, top=0.96)
        plt.savefig(stellFile+"_GammaOmegaKy.png")
    plt.close()
    return np.array(kyX), np.array(growthRateX), np.array(realFrequencyX)

# Function to replace text in a file
def replace(file_path, pattern, subst):
    fh, abs_path = mkstemp()
    with fdopen(fh,'w') as new_file:
        with open(file_path) as old_file:
            for line in old_file:
                new_file.write(line.replace(pattern, subst))
    copymode(file_path, abs_path)
    remove(file_path)
    move(abs_path, file_path)

def output_to_csv(output_csv, growth_rate, omega, ky, weighted_growth_rate, ln, lt):
    keys=np.concatenate([['ln'],['lt'],['growth_rate'],['omega'],['ky'], ['weighted_growth_rate']])
    values=np.concatenate([[ln],[lt],[growth_rate],[omega],[ky],[weighted_growth_rate]])
    dictionary = dict(zip(keys, values))
    df = pd.DataFrame(data=[dictionary])
    if not os.path.exists(output_csv): pd.DataFrame(columns=df.columns).to_csv(output_csv, index=False)
    df.to_csv(output_csv, mode='a', header=False, index=False)


def remove_files():
    for f in glob.glob('*.amoments'): remove(f)
    for f in glob.glob('*.eigenfunc'): remove(f)
    for f in glob.glob('*.error'): remove(f)
    for f in glob.glob('*.fields'): remove(f)
    for f in glob.glob('*.g'): remove(f)
    for f in glob.glob('*.lpc'): remove(f)
    for f in glob.glob('*.mom2'): remove(f)
    for f in glob.glob('*.moments'): remove(f)
    for f in glob.glob('*.vres'): remove(f)
    for f in glob.glob('*.vres2'): remove(f)
    for f in glob.glob('*.exit_reason'): remove(f)
    for f in glob.glob('*.optim'): remove(f)
    for f in glob.glob('*.out'): remove(f)
    for f in glob.glob('*.scratch'): remove(f)
    for f in glob.glob('*.used_inputs.in'): remove(f)
    for f in glob.glob('*.vspace_integration_error'): remove(f)
    ## THIS SHOULD ONLY REMOVE FILES STARTING WTH .gs2
    for f in glob.glob('.gs2*'): remove(f)
    ## REMOVE ALSO INPUT FILES
    for f in glob.glob('*.in'): remove(f)
    ## REMOVE ALSO OUTPUT FILES
    for f in glob.glob('*.out.nc'):
        if f not in 'gs2Input-LN1.0-LT3.0.out.nc': remove(f)