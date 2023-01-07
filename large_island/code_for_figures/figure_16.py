
import sys
sys.path.append('../src')

import matplotlib.pyplot as plt
import numpy as np

from ms_toolbox import gfd
from ms_toolbox import rgmodel as rgm
from pycurrents.num import spectra
from pycurrents.file.matfile import loadmatbunch

from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy.interpolate import interp2d as int2d

from netCDF4 import Dataset
from pycurrents.system import Bunch
from matplotlib.ticker import FormatStrFormatter

import cmocean
import matplotlib.cm as cm
import matplotlib
csfont = {'fontname':'Times New Roman'}
matplotlib.rcParams['mathtext.fontset'] = 'stix'

"""

"""

ff = "../large_island/data_sub/"
p1 = ff + 'r300/'
p2 = ff + 'r300_longforce/'

figpath = '../large_island/figures/'


def load_var(varpath,n_filename):
    aN = np.load(varpath + n_filename)

    b = Bunch()
    b.uE = aN["uE"][:]
    b.vE = aN["vE"][:]
    b.hE = aN["hE"][:]
    b.etaE = aN["etaE"][:]

    b.uff1= aN["uff1"][:]
    b.vff1= aN["vff1"][:]
    b.hff1 = aN["hff1"][:]
    b.etaff1 = aN["etaff1"][:]

    b.uff2= aN["uff2"][:]
    b.vff2= aN["vff2"][:]
    b.hff2 = aN["hff2"][:]
    b.etaff2 = aN["etaff2"][:]

    b.t = aN["t"][:]

    fcor = gfd.coriolis(8)
    b.fcor = fcor

    return b

################# Read Grid ###################

file_name1 = "complex_demod_ke_x.npz"

BA = load_var(p1,file_name1)
r75 = load_var(p2,file_name1)


def spec_r(TS,dt,nsmooth):
    spec = spectra.spectrum(TS,dt=dt, nfft=None, window='quadratic')
    nmaxi = np.argmax(spec.cwpsd)
    fmax = spec.freqs[nmaxi]  ### tw frequency

    return spec,fmax


dtdays = BA.t[1] -  BA.t[0]
nsmooth=3


svname = "FIGURE_16.png"
fig,axs = plt.subplots(figsize=(10,10),ncols=2,nrows=2)

TS1 = BA.uE+1j*BA.vE
spec,fmax = spec_r(TS1,dtdays,nsmooth)
TS2 = r75.uE+1j*r75.vE
spec2,fmax = spec_r(TS2,dtdays,nsmooth)


ax = axs[0,1]
ax.loglog(spec.cwfreqs,spec.cwpsd,'m',label="1 Inertial Period Forcing (CW)")
ax.loglog(spec2.cwfreqs,spec2.cwpsd,'b',label="3 Inertial Period Forcing (CW)")
ax.loglog(spec.ccwfreqs,spec.ccwpsd,'m--',label="1 Inertial Period Forcing (CCW)")
ax.loglog(spec2.ccwfreqs,spec2.ccwpsd,'b--',label="3 Inertial Period Forcing (CCW)")

ax.plot([BA.fcor[1],BA.fcor[1]],[1e-8,1e0],'--k',label="$f_{in}$")
ax.plot([0.1118,0.1118],[1e-8,1e0],'--k',label="$f_{tw}$")
ax.set_xlabel("cpd",fontsize=16)
ax.set_title("Rotary Spectrum (island boundary)",fontsize=16)
ax.set_ylabel("$\\mathrm{m^2s^{-2}} \\mathrm{cpd^{-1}}$",fontsize=16)
ax.set_ylim([1e-8,1e0])
ax.set_xlim([2e-2,1e0])
ax.text(BA.fcor[1],1e-7,"$f$",fontsize=16)
ax.text(.11,1e-7,"$f_{tw}$",fontsize=16)
ax.tick_params("both",labelsize=14)
ax.text(9e-3,2e0,"b)",fontsize=16)

TS1 = BA.uff2+1j*BA.vff2
spec,fmax = spec_r(TS1,dtdays,nsmooth)
TS2 = r75.uff2+1j*r75.vff2
spec2,fmax = spec_r(TS2,dtdays,nsmooth)


ax = axs[0,0]
ax.loglog(spec.cwfreqs,spec.cwpsd,'m',label="1IPF (CW)")
ax.loglog(spec2.cwfreqs,spec2.cwpsd,'b',label="3IPF (CW)")
ax.loglog(spec.ccwfreqs,spec.ccwpsd,'m--',label="1IPF (CCW)")
ax.loglog(spec2.ccwfreqs,spec2.ccwpsd,'b--',label="3IPF (CCW)")
ax.plot([BA.fcor[1],BA.fcor[1]],[1e-8,1e0],'--k')
ax.plot([0.1118,0.1118],[1e-8,1e0],'--k')
ax.set_xlabel("cpd",fontsize=16)
ax.set_title("Rotary Spectrum (far-field)",fontsize=16)
ax.set_ylim([1e-8,1e0])
ax.set_xlim([2e-2,1e0])
ax.text(BA.fcor[1],1e-7,"$f$",fontsize=16)
ax.text(.11,1e-7,"$f_{tw}$",fontsize=16)
ax.set_ylabel("$\\mathrm{m^2s^{-2}} \\mathrm{cpd^{-1}}$",fontsize=16)
ax.tick_params("both",labelsize=14)
ax.text(9e-3,2e0,"a)",fontsize=16)
ax.legend(loc="upper left",fontsize=12)



def spec_r(TS,dt,nsmooth):
    spec = spectra.spectrum(TS,dt=dt, nfft=None, window='quadratic')
    nmaxi = np.argmax(spec.psd)
    fmax = spec.freqs[nmaxi]  ### tw frequency

    return spec,fmax


TS1 = BA.etaE
spec,fmax = spec_r(TS1,dtdays,nsmooth)
TS2 = r75.etaE
spec2,fmax = spec_r(TS2,dtdays,nsmooth)


ax = axs[1,1]

ax.loglog(spec.freqs,spec.psd,'m',label="1 Inertial Period Forcing")
ax.loglog(spec2.freqs,spec2.psd,'b',label="3 Inertial Period Forcing")

ax.plot([BA.fcor[1],BA.fcor[1]],[1e-12,1e-1],'--k')
ax.plot([0.1118,0.1118],[1e-12,1e-1],'--k')
ax.set_xlabel("cpd",fontsize=16)
ax.set_title("$\\eta$ (island boundary)",fontsize=16)
ax.set_ylabel("$\\mathrm{m^2} \\mathrm{cpd^{-1}}$",fontsize=16)
ax.set_ylim([1e-12,1e-1])
ax.set_xlim([2e-2,1e0])
ax.tick_params("both",labelsize=14)
ax.text(9e-3,2e-1,"d)",fontsize=16)
ax.text(BA.fcor[1],5e-12,"$f$",fontsize=16)
ax.text(.11,5e-12,"$f_{tw}$",fontsize=16)

ax = axs[1,0]

TS1 = BA.etaff2
spec,fmax = spec_r(TS1,dtdays,nsmooth)
TS2 = r75.etaff2
spec2,fmax = spec_r(TS2,dtdays,nsmooth)


ax.loglog(spec.freqs,spec.psd,'m',label="1IPF")
ax.loglog(spec2.freqs,spec2.psd,'b',label="3IPF")
ax.plot([BA.fcor[1],BA.fcor[1]],[1e-12,1e-1],'--k')
ax.plot([0.1118,0.1118],[1e-12,1e-1],'--k')
ax.set_xlabel("cpd",fontsize=16)
ax.set_title("$\\eta$ (far-field)",fontsize=16)
ax.set_ylabel("$\\mathrm{m^2} \\mathrm{cpd^{-1}}$",fontsize=16)
ax.set_ylim([1e-12,1e-1])
ax.set_xlim([2e-2,1e0])
ax.tick_params("both",labelsize=14)
ax.text(9e-3,2e-1,"c)",fontsize=16)
ax.text(BA.fcor[1],5e-12,"$f$",fontsize=16)
ax.text(.11,5e-12,"$f_{tw}$",fontsize=16)
ax.legend(loc="upper left",fontsize=12)

fig.tight_layout(pad=1)


fig.savefig(figpath + svname, dpi=300, bbox_inches='tight',orientation="horizontal",pad_inches=.3)
