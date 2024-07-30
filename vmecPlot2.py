#!/usr/bin/env python3
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
from scipy.io import netcdf_file
import os
import math
import sys
from pathlib import Path

def main(file,name='',figures_folder='.', coils_curves=None, s_plot_ignore=0.2,savefig=True):

    filename = file
    if name=='': name=os.path.basename(filename)[5:-3]
    f = netcdf_file(filename,'r',mmap=False)
    phi = f.variables['phi'][()]
    iotaf = f.variables['iotaf'][()]
    presf = f.variables['presf'][()]
    iotas = f.variables['iotas'][()]
    pres = f.variables['pres'][()]
    ns = f.variables['ns'][()]
    nfp = f.variables['nfp'][()]
    xn = f.variables['xn'][()]
    xm = f.variables['xm'][()]
    xn_nyq = f.variables['xn_nyq'][()]
    xm_nyq = f.variables['xm_nyq'][()]
    rmnc = f.variables['rmnc'][()]
    zmns = f.variables['zmns'][()]
    bmnc = f.variables['bmnc'][()]
    raxis_cc = f.variables['raxis_cc'][()]
    zaxis_cs = f.variables['zaxis_cs'][()]
    buco = f.variables['buco'][()]
    bvco = f.variables['bvco'][()]
    jcuru = f.variables['jcuru'][()]
    jcurv = f.variables['jcurv'][()]
    lasym = f.variables['lasym__logical__'][()]
    if lasym==1:
        rmns = f.variables['rmns'][()]
        zmnc = f.variables['zmnc'][()]
        bmns = f.variables['bmns'][()]
        raxis_cs = f.variables['raxis_cs'][()]
        zaxis_cc = f.variables['zaxis_cc'][()]
    else:
        rmns = 0*rmnc
        zmnc = 0*rmnc
        bmns = 0*bmnc
        raxis_cs = 0*raxis_cc
        zaxis_cc = 0*raxis_cc

    try:
        ac = f.variables['ac'][()]
    except:
        ac = []

    try:
        pcurr_type = f.variables['pcurr_type'][()]
    except:
        pcurr_type = ""

    ac_aux_s = f.variables['ac_aux_s'][()]
    ac_aux_f = f.variables['ac_aux_f'][()]

    #print type(pcurr_type)
    #print pcurr_type
    #print str(pcurr_type)
    #exit(0)

    print("nfp: ",nfp)
    print("ns: ",ns)

    mpol = f.variables['mpol'][()]
    print("mpol: ",mpol)

    ntor = f.variables['ntor'][()]
    print("ntor: ",ntor)

    Aminor_p = f.variables['Aminor_p'][()]
    print("Aminor_p: ",Aminor_p)

    Rmajor_p = f.variables['Rmajor_p'][()]
    print("Rmajor_p: ",Rmajor_p)

    data = f.variables['aspect'][()]
    print("aspect:            ",data)
    print("Rmajor_p/Aminor_p: ",Rmajor_p/Aminor_p)

    data = f.variables['betatotal'][()]
    print("betatotal: ",data)

    data = f.variables['betapol'][()]
    print("betapol:   ",data)

    data = f.variables['betator'][()]
    print("betator:   ",data)

    data = f.variables['betaxis'][()]
    print("betaxis:   ",data)

    ctor = f.variables['ctor'][()]
    print("ctor:   ",ctor)

    DMerc = f.variables['DMerc'][()]
    #print("DMerc:   ",DMerc)

    f.close()
    nmodes = len(xn)

    s = np.linspace(0,1,ns)
    s_half = [(i-0.5)/(ns-1) for i in range(1,ns)]

    phiedge = phi[-1]
    phi_half = [(i-0.5)*phiedge/(ns-1) for i in range(1,ns)]

    ntheta = 200
    nzeta = 8
    theta = np.linspace(0,2*np.pi,num=ntheta)
    zeta = np.linspace(0,2*np.pi/nfp,num=nzeta,endpoint=False)
    iradius = ns-1
    R = np.zeros((ntheta,nzeta))
    Z = np.zeros((ntheta,nzeta))
    for itheta in range(ntheta):
        for izeta in range(nzeta):
            for imode in range(nmodes):
                angle = xm[imode]*theta[itheta] - xn[imode]*zeta[izeta]
                R[itheta,izeta] = R[itheta,izeta] + rmnc[iradius,imode]*math.cos(angle) + rmns[iradius,imode]*math.sin(angle)
                Z[itheta,izeta] = Z[itheta,izeta] + zmns[iradius,imode]*math.sin(angle) + zmnc[iradius,imode]*math.cos(angle)

    Raxis = np.zeros(nzeta)
    Zaxis = np.zeros(nzeta)
    for izeta in range(nzeta):
        for n in range(ntor+1):
            angle = -n*nfp*zeta[izeta]
            Raxis[izeta] += raxis_cc[n]*math.cos(angle) + raxis_cs[n]*math.sin(angle)
            Zaxis[izeta] += zaxis_cs[n]*math.sin(angle) + zaxis_cc[n]*math.cos(angle)

    xLabel = r'$s = \psi/\psi_b$'


    fig = plt.figure()
    fig.set_size_inches(14,7)
    fig.patch.set_facecolor('white')

    numCols = 3
    numRows = 3
    plotNum = 1

    plt.subplot(numRows,numCols,plotNum)
    plotNum += 1
    #plt.plot(phi, iotaf, '.-',label='iotaf')
    #plt.plot(phi_half, iotas[1:],'.-',label='iotas')
    plt.plot(s, iotaf, '.-',label='iota')
    # plt.plot(s_half, iotas[1:],'.-',label='iotas')
    # plt.legend(fontsize='x-small')
    plt.xlabel(xLabel)
    plt.ylabel(f'$\iota$')

    plt.subplot(numRows,numCols,plotNum)
    plotNum += 1
    #plt.plot(phi, presf, '.-',label='presf')
    #plt.plot(phi_half, pres[1:], '.-',label='pres')
    plt.plot(s, presf, '.-',label='presf')
    plt.plot(s_half, pres[1:], '.-',label='pres')
    plt.legend(fontsize='x-small')
    plt.xlabel(xLabel)

    plt.subplot(numRows,numCols,plotNum)
    plotNum += 1
    #plt.plot(phi_half, buco[1:], '.-',label='buco')
    plt.plot(s_half, buco[1:], '.-',label='buco')
    plt.title('buco')
    plt.xlabel(xLabel)

    plt.subplot(numRows,numCols,plotNum)
    plotNum += 1
    #plt.plot(phi_half, bvco[1:], '.-',label='bvco')
    plt.plot(s_half, bvco[1:], '.-',label='bvco')
    plt.title('bvco')
    plt.xlabel(xLabel)

    plt.subplot(numRows,numCols,plotNum)
    plotNum += 1
    #plt.plot(phi, jcuru, '.-',label='jcuru')
    plt.plot(s, jcuru, '.-',label='jcuru')
    plt.title('jcuru')
    plt.xlabel(xLabel)

    plt.subplot(numRows,numCols,plotNum)
    plotNum += 1
    #plt.plot(phi, jcurv, '.-',label='jcurv')
    plt.plot(s, jcurv, '.-',label='jcurv')
    plt.title('jcurv')
    plt.xlabel(xLabel)

    plt.subplot(numRows,numCols,plotNum)
    plotNum += 1
    # if 'power_series' in pcurr_type:
    # 	ac_profile = phi*0.0
    # 	for i in range(len(ac)):
    # 		ac_profile += ac[i]*(s**i)
    # 	plt.plot(s, ac_profile, '.-')
    # else:
    # 	mask = (ac_aux_s >= 0)
    # 	plt.plot(ac_aux_s[mask], ac_aux_f[mask],'.-')
    # plt.title('ac profile')
    # plt.xlabel(xLabel)
    plt.plot(s[int(s_plot_ignore*len(s)):-2],DMerc[int(s_plot_ignore*len(s)):-2],'.-')
    plt.title('DMerc')
    plt.xlabel(xLabel)

    titles = ['|B| at half radius','|B| at LCFS']
    iradii = [int((ns*0.25).round()), ns-1]
    print("bmnc.shape:",bmnc.shape)
    print("bmns.shape:",bmns.shape)
    for i in range(2):
        iradius = iradii[i]
        Ntheta = 30
        Nzeta = 65
        theta = np.linspace(0,2*np.pi,num=Ntheta)
        zeta = np.linspace(0,2*np.pi,num=Nzeta)
        b = np.zeros([Ntheta,Nzeta])
        zeta2D,theta2D = np.meshgrid(zeta,theta)
        iota = iotaf[iradius]
        for imode in range(len(xn_nyq)):
            angle = xm_nyq[imode]*theta2D - xn_nyq[imode]*zeta2D
            b += bmnc[iradius,imode]*np.cos(angle) + bmns[iradius,imode]*np.sin(angle)

        plt.subplot(numRows,numCols,plotNum)
        plotNum += 1
        plt.contourf(zeta2D,theta2D,b,20)
        plt.title(titles[i]+'\n1-based index='+str(iradius+1))
        plt.xlabel('zeta')
        plt.ylabel('theta')
        plt.colorbar()
        # Plot a field line:
        if iota>0:
            plt.plot([0,zeta.max()],[0,zeta.max()*iota],'k')
        else:
            plt.plot([0,zeta.max()],[-zeta.max()*iota,0],'k')
        plt.xlim([0,2*np.pi])
        plt.ylim([0,2*np.pi])

    plt.tight_layout()
    plt.figtext(0.5,0.99,os.path.abspath(filename),ha='center',va='top',fontsize=6)

    if savefig: plt.savefig(os.path.join(figures_folder, name+'_VMECparams.pdf'), bbox_inches = 'tight', pad_inches = 0)

    ########################################################
    # Now make plot of flux surface shapes
    ########################################################

    fig = plt.figure()
    fig.set_size_inches(6,6)
    ax=fig.add_subplot(111, label="1")
    fig.patch.set_facecolor('white')
    plt.plot(R[:,0], Z[:,0], '-',label=r'$\phi$=0')
    # plt.plot(R[:,1], Z[:,1], '-',label='_nolegend_')
    plt.plot(R[:,2], Z[:,2], '-',label=r'$\phi=\pi/2$')
    # plt.plot(R[:,3], Z[:,3], '-',label='_nolegend_')
    plt.plot(R[:,4], Z[:,4], '-',label=r'$\phi=\pi$')
    # plt.plot(R[:,5], Z[:,5], '-',label='_nolegend_')
    plt.plot(R[:,6], Z[:,6], '-',label=r'$\phi=3\pi/2$')
    # plt.plot(R[:,7], Z[:,7], '-',label='_nolegend_')
    plt.gca().set_aspect('equal',adjustable='box')
    plt.legend(fontsize=18)
    plt.xlabel('R', fontsize=22)
    plt.ylabel('Z', fontsize=22)
    ax.tick_params(axis='x', labelsize=16)
    ax.tick_params(axis='y', labelsize=16)
    plt.tight_layout()
    if savefig: plt.savefig(os.path.join(figures_folder, name+'_poloidal_plot.png'))
    R_boundary = R
    Z_boundary = Z

    fig = plt.figure()
    fig.set_size_inches(14,7)
    fig.patch.set_facecolor('white')
    # plt.subplot(numRows,numCols,plotNum)
    numCols = 4
    numRows = 2
    plotNum = 1
    ntheta = 200
    nzeta = 8
    nradius = 8
    theta = np.linspace(0,2*np.pi,num=ntheta)
    zeta = np.linspace(0,2*np.pi/nfp,num=nzeta,endpoint=False)
    iradii = np.linspace(0,ns-1,num=nradius).round()
    iradii = [int(i) for i in iradii]
    R = np.zeros((ntheta,nzeta,nradius))
    Z = np.zeros((ntheta,nzeta,nradius))
    for itheta in range(ntheta):
        for izeta in range(nzeta):
            for iradius in range(nradius):
                for imode in range(nmodes):
                    angle = xm[imode]*theta[itheta] - xn[imode]*zeta[izeta]
                    R[itheta,izeta,iradius] = R[itheta,izeta,iradius] + rmnc[iradii[iradius],imode]*math.cos(angle) \
                                                                    + rmns[iradii[iradius],imode]*math.sin(angle)
                    Z[itheta,izeta,iradius] = Z[itheta,izeta,iradius] + zmns[iradii[iradius],imode]*math.sin(angle) \
                                                                    + zmnc[iradii[iradius],imode]*math.cos(angle)

    for izeta in range(nzeta):
        plt.subplot(numRows,numCols,plotNum)
        plotNum += 1
        for iradius in range(nradius):
            plt.plot(R[:,izeta,iradius], Z[:,izeta,iradius], '-')
        plt.plot(Raxis[izeta],Zaxis[izeta],'xr')
        plt.gca().set_aspect('equal',adjustable='box')
        plt.xlabel('R', fontsize=10)
        plt.ylabel('Z', fontsize=10)
        plt.title(r'$\phi$ = '+str(round(zeta[izeta],2)))

    plt.tight_layout()
    # plt.subplots_adjust(wspace=0, hspace=0)
    # plt.figtext(0.5,0.99,os.path.abspath(filename),ha='center',va='top',fontsize=6)
    if savefig: plt.savefig(os.path.join(figures_folder, name+'_VMECsurfaces.pdf'), bbox_inches = 'tight', pad_inches = 0)

    ########################################################
    # Now make 3D surface plot
    ########################################################

    fig = plt.figure(figsize=(5, 4), frameon=False)

    ntheta = 80
    nzeta = np.max((500,int(150*nfp)))
    theta1D = np.linspace(0,2*np.pi,num=ntheta)
    zeta1D = np.linspace(0,2*np.pi,num=nzeta)
    zeta2D, theta2D = np.meshgrid(zeta1D,theta1D)
    iradius = ns-1
    R = np.zeros((ntheta,nzeta))
    Z = np.zeros((ntheta,nzeta))
    B = np.zeros((ntheta,nzeta))
    for imode in range(nmodes):
        angle = xm[imode]*theta2D - xn[imode]*zeta2D
        R = R + rmnc[iradius,imode]*np.cos(angle) + rmns[iradius,imode]*np.sin(angle)
        Z = Z + zmns[iradius,imode]*np.sin(angle) + zmnc[iradius,imode]*np.cos(angle)

    for imode in range(len(xn_nyq)):
        angle = xm_nyq[imode]*theta2D - xn_nyq[imode]*zeta2D
        B = B + bmnc[iradius,imode]*np.cos(angle) + bmns[iradius,imode]*np.sin(angle)

    X = R * np.cos(zeta2D)
    Y = R * np.sin(zeta2D)
    # Rescale to lie in [0,1]:
    B_rescaled = (B - B.min()) / (B.max() - B.min())

    ax = fig.add_subplot(111, projection='3d')  # Create 3D axes
    surf = ax.plot_surface(X, Y, Z, facecolors = cm.jet(B_rescaled), rstride=1, cstride=1, antialiased=False)
    ax.auto_scale_xyz([0.7*X.min(), 0.7*X.max()], [0.7*X.min(), 0.7*X.max()], [0.7*X.min(), 0.7*X.max()])
    # plt.tight_layout()
    ax.set_box_aspect([1, 1, 1])
    
    # Remove axis
    ax.axis('off')
    
    from matplotlib.colors import Normalize
    cax = fig.add_axes([0.21, 0.8, 0.6, 0.03])  # Adjust position and size as needed
    norm = Normalize(vmin=B.min(), vmax=B.max())
    sm = cm.ScalarMappable(cmap=cm.jet, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, orientation='horizontal', cax=cax)
    cbar.set_label('|B| [T]')

    # plt.figtext(0.5,0.99,os.path.abspath(filename),ha='center',va='top',fontsize=6)
    if savefig: plt.savefig(os.path.join(figures_folder, name+'_VMEC_3Dplot.png'), bbox_inches = 'tight', pad_inches = 0, dpi=400)
    else: plt.show()

    plt.close()
    # exit()
if __name__ == "__main__":
    # Create results folders if not present
    try:
        Path(sys.argv[2]).mkdir(parents=True, exist_ok=True)
        figures_results_path = str(Path(sys.argv[2]).resolve())
        main(sys.argv[1], sys.argv[2])
    except:
        main(sys.argv[1])