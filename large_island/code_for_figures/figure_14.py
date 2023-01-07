
import sys
sys.path.append('../src')

import os
import matplotlib.pyplot as plt
import numpy as np

from ms_toolbox import gfd
from ms_toolbox import rgmodel as rgm
from pycurrents.num import spectra, interp1
from netCDF4 import Dataset

from pycurrents.file.matfile import loadmatbunch
from scipy.interpolate import interp2d as int2d

import matplotlib
csfont = {'fontname':'Times New Roman'}

matplotlib.rcParams['mathtext.fontset'] = 'stix'

figpath = "../large_island/figures/" #path to git repository


dt = 120
tstart=0   # seconds
tend=52 * 86400 # Number of days to run model>> seconds
theta0 = 8

times=np.arange(tstart, tend+dt, dt)

nstart = 2

tdays = times / 86400
dtdays = tdays[1] - tdays[0]
theta0=8
fcor = gfd.coriolis(theta0)

amp1 = .05
amp2 = .05 / 3

nevent1 = 1  # number of inertial cycles
nevent2 = 3  # number of inertial cycles


omega = fcor[0]
taux1,tauy1 = rgm.make_wind(nevent1,omega, amp1, times,nstart,taper=True,wind_type="CW")
taux2,tauy2 = rgm.make_wind(nevent2,omega, amp2, times,nstart,taper=True,wind_type="CW")

spec1 = spectra.spectrum(taux1+1j*tauy1,dt=dtdays, nfft=None,window='quadratic')
spec2 = spectra.spectrum(taux2+1j*tauy2,dt=dtdays, nfft=None,window='quadratic')

fig,ax = plt.subplots(figsize=(10,8))
ax.loglog(spec1.cwfreqs,spec1.cwps,'m',label="1IPF (CW)")
ax.loglog(spec2.cwfreqs,spec2.cwps,'b',label="3IPF (CW)")
ax.loglog(spec1.ccwfreqs,spec1.ccwps,'m--',label="1IPF (CCW)")
ax.loglog(spec2.ccwfreqs,spec2.ccwps,'b--',label="3IPF (CCW)")
ax.plot([fcor[1],fcor[1]],[1e-33,1e-2],'--k')
# ax.set_xlim([np.min(spec1.freqs),np.max(spec1.freqs)])
ax.set_ylim([1e-15,1e-2])
ax.text(fcor[1],8e-15,"$f$",fontsize=20)
# ax.plot([0.1118,0.1118],[1e-8,1e-1],'--k',label="$f_{tw}$")
ax.set_ylabel("$\\mathrm{N^{2}m^{-4}cpd^{-1}}$",fontsize=16)
ax.legend(fontsize=16)
ax.tick_params("both",labelsize=16)
ax.set_xlabel("cpd",fontsize=16)
ax.set_title("Power Spectral Density: $\\vec{\\tau}$ ",fontsize=16)
ax.set_xlim([2e-2,1e0])

fig.savefig(figpath + 'FIGURE_14' + ".png", dpi=300, bbox_inches='tight',orientation="horizontal")
