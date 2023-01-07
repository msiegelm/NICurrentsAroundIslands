
import sys
sys.path.append('../src')

import rg_plotting as rgp
import animate_plotting as ap
import a_tools as at

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
matplotlib.rcParams['mathtext.fontset'] = 'stix'
"""

"""

ff = '../small_islands/compare/data_subsection/'
p1 = ff + 'r50/'
p2 = ff + 'r75/'
p3 = ff + 'r100/'
p4 = ff + 'r150/'
p5 = ff + 'r200/'

figpath = '../small_islands/compare/figures/'


def load_var(varpath,n_filename):
    aN = np.load(varpath + n_filename)

    b = Bunch()
    b.uN = aN["uN"][:431]
    b.vN = aN["vN"][:431]
    b.hN = aN["hN"][:431]
    b.etaN = aN["etaN"][:431]

    b.uff1= aN["uff1"][:431]
    b.vff1= aN["vff1"][:431]
    b.hff1 = aN["hff1"][:431]
    b.etaff1 = aN["etaff1"][:431]

    b.uff2= aN["uff2"][:431]
    b.vff2= aN["vff2"][:431]
    b.hff2 = aN["hff2"][:431]
    b.etaff2 = aN["etaff2"][:431]

    b.t = aN["t"][:431]

    fcor = gfd.coriolis(8)
    b.fcor = fcor

    return b

################# Read Grid ###################

file_name1 = "North_ts_new.npz"

varpath = p1
BA = load_var(varpath,file_name1)

varpath = p2
r75 = load_var(varpath,file_name1)

varpath = p3
r100 = load_var(varpath,file_name1)

varpath = p4
r150 = load_var(varpath,file_name1)

varpath = p5
r200 = load_var(varpath,file_name1)




H = 100
gp=.09
fcor = gfd.coriolis(8)
Ld=np.sqrt(gp*H) / fcor[0]

eps50 = (50e3**2) / (Ld**2)
eps75 = (75e3**2) / (Ld**2)
eps100 = (100e3**2) / (Ld**2)
eps150 = (150e3**2) / (Ld**2)
eps200 = (200e3**2) / (Ld**2)




def spec_r(TS,dt,nsmooth):
    spec = spectra.spectrum(TS,dt=dt, nfft=None,window='quadratic')
    nmaxi = np.argmax(spec.psd)
    fmax = spec.freqs[nmaxi]  ### peak frequency

    return spec,fmax






dtdays = BA.t[1] -  BA.t[0]
nsmooth=1


#%%

svname = "FIGURE_12.png"
fig,ax = plt.subplots(figsize=(7,7),ncols=2,nrows=2)

TS1 = BA.uN
spec,fmax = spec_r(TS1,dtdays,nsmooth)
TS2 = r75.uN
spec2,fmax2 = spec_r(TS2,dtdays,nsmooth)
TS3 = r100.uN
spec3,fmax3 = spec_r(TS3,dtdays,nsmooth)
TS4 = r150.uN
spec4,fmax4 = spec_r(TS4,dtdays,nsmooth)
TS5 = r200.uN
spec5,fmax5 = spec_r(TS5,dtdays,nsmooth)


ax[0,0].loglog(spec.freqs,spec.psd,'r',label=("$\\epsilon_L =$ %.02f" % eps50 ),color="skyblue")
ax[0,0].loglog(spec2.freqs,spec2.psd,'g',label=("$\\epsilon_L =$ %.02f" % eps75 ),color="dodgerblue")
ax[0,0].loglog(spec3.freqs,spec3.psd,'b',label=("$\\epsilon_L =$ %.02f" % eps100 ),color="royalblue")
ax[0,0].loglog(spec4.freqs,spec4.psd,'c',label=("$\\epsilon_L =$ %.02f" % eps150 ),color="blue")
ax[0,0].loglog(spec5.freqs,spec5.psd,'m',label=("$\\epsilon_L =$ %.02f" % eps200 ),color="darkblue")
ax[0,0].plot([BA.fcor[1],BA.fcor[1]],[1e-12,1],'--k')
ax[0,0].legend(frameon=False,bbox_to_anchor=(0.54, 0.48))
ax[0,0].set_xlabel("$\\mathrm{cpd}$",fontsize=14)
ax[0,0].set_title("Northern Island Boundary: $u$",fontsize=16)
ax[0,0].set_ylabel("$\\mathrm{m^2 \ s^{-2} cpd^{-1}}$",fontsize=14)
ax[0,0].set_xlim([spec.freqs[0],1])
ax[0,0].set_ylim([1e-6,1])
ax[0,0].text(BA.fcor[1]+.2e-1,3e-6,"$f$",fontsize=14)
ax[0,0].tick_params(axis="both",labelsize=14)
ax[0,0].text(2e-2,3,"a)",fontsize=16)
TS1 = BA.uff2
spec,fmax = spec_r(TS1,dtdays,nsmooth)
TS2 = r75.uff2
spec2,fmax2 = spec_r(TS2,dtdays,nsmooth)
TS3 = r100.uff2
spec3,fmax3 = spec_r(TS3,dtdays,nsmooth)

TS4 = r150.uff2
spec4,fmax4 = spec_r(TS4,dtdays,nsmooth)

TS5 = r200.uff2
spec5,fmax5 = spec_r(TS5,dtdays,nsmooth)

ax[0,1].loglog(spec.freqs,spec.psd,'r',label=("$\\epsilon_L =$ %.02f" % eps50 ),color="skyblue")
ax[0,1].loglog(spec2.freqs,spec2.psd,'g',label=("$\\epsilon_L =$ %.02f" % eps75 ),color="dodgerblue")
ax[0,1].loglog(spec3.freqs,spec3.psd,'b',label=("$\\epsilon_L =$ %.02f" % eps100 ),color="royalblue")
ax[0,1].loglog(spec4.freqs,spec4.psd,'c',label=("$\\epsilon_L =$ %.02f" % eps150 ),color="blue")
ax[0,1].loglog(spec5.freqs,spec5.psd,'m',label=("$\\epsilon_L =$ %.02f" % eps200 ),color="darkblue")
ax[0,1].plot([BA.fcor[1],BA.fcor[1]],[1e-12,1],'--k')
ax[0,1].legend(frameon=False,bbox_to_anchor=(0.525, 0.536))
ax[0,1].set_xlabel("$\\mathrm{cpd}$",fontsize=14)
ax[0,1].set_title("Far-field: $u$",fontsize=16)
ax[0,1].set_ylabel("$\\mathrm{m^2 \ s^{-2} cpd^{-1}}$",fontsize=14)
ax[0,1].set_xlim([spec.freqs[0],1])
ax[0,1].set_ylim([1e-6,1])
ax[0,1].text(BA.fcor[1]+.2e-1,3e-6,"$f$",fontsize=14)
ax[0,1].tick_params(axis="both",labelsize=14)
ax[0,1].text(2e-2,3,"b)",fontsize=16)

TS1 = BA.etaN
spec,fmax = spec_r(TS1,dtdays,nsmooth)
TS2 = r75.etaN
spec2,fmax2 = spec_r(TS2,dtdays,nsmooth)
TS3 = r100.etaN
spec3,fmax3 = spec_r(TS3,dtdays,nsmooth)
TS4 = r150.etaN
spec4,fmax4 = spec_r(TS4,dtdays,nsmooth)
TS5 = r200.etaN
spec5,fmax5 = spec_r(TS5,dtdays,nsmooth)


ax[1,0].loglog(spec.freqs,spec.psd,'r',label=("$\\epsilon_L =$ %.02f" % eps50 ),color="skyblue")
ax[1,0].loglog(spec2.freqs,spec2.psd,'g',label=("$\\epsilon_L =$ %.02f" % eps75 ),color="dodgerblue")
ax[1,0].loglog(spec3.freqs,spec3.psd,'b',label=("$\\epsilon_L =$ %.02f" % eps100 ),color="royalblue")
ax[1,0].loglog(spec4.freqs,spec4.psd,'c',label=("$\\epsilon_L =$ %.02f" % eps150 ),color="blue")
ax[1,0].loglog(spec5.freqs,spec5.psd,'m',label=("$\\epsilon_L =$ %.02f" % eps200 ),color="darkblue")
ax[1,0].plot([BA.fcor[1],BA.fcor[1]],[1e-12,500],'--k')
ax[1,0].legend(frameon=False)
ax[1,0].set_xlabel("$\\mathrm{cpd}$",fontsize=14)
ax[1,0].set_title("Northern Island Boundary: $\\eta$",fontsize=16)
ax[1,0].set_ylabel("$\\mathrm{m^2 \ cpd^{-1}}$",fontsize=14)
ax[1,0].set_xlim([spec.freqs[0],1])
ax[1,0].set_ylim([1e-12,1e-1])
ax[1,0].text(BA.fcor[1]+.2e-1,5e-12,"$f$",fontsize=14)
ax[1,0].tick_params(axis="both",labelsize=14)
ax[1,0].text(2e-2,3e-1,"c)",fontsize=16)

TS1 = BA.etaff2
spec,fmax = spec_r(TS1,dtdays,nsmooth)
TS2 = r75.etaff2
spec2,fmax2 = spec_r(TS2,dtdays,nsmooth)
TS3 = r100.etaff2
spec3,fmax3 = spec_r(TS3,dtdays,nsmooth)
TS4= r150.etaff2
spec4,fmax4 = spec_r(TS4,dtdays,nsmooth)
TS5 = r200.etaff2
spec5,fmax5 = spec_r(TS5,dtdays,nsmooth)

ax[1,1].loglog(spec.freqs,spec.psd,'r',label=("$\\epsilon_L =$ %.02f" % eps50 ),color="skyblue")
ax[1,1].loglog(spec2.freqs,spec2.psd,'g',label=("$\\epsilon_L =$ %.02f" % eps75 ),color="dodgerblue")
ax[1,1].loglog(spec3.freqs,spec3.psd,'b',label=("$\\epsilon_L =$ %.02f" % eps100 ),color="royalblue")
ax[1,1].loglog(spec4.freqs,spec4.psd,'c',label=("$\\epsilon_L =$ %.02f" % eps150 ),color="blue")
ax[1,1].loglog(spec5.freqs,spec5.psd,'m',label=("$\\epsilon_L =$ %.02f" % eps200 ),color="darkblue")
ax[1,1].plot([BA.fcor[1],BA.fcor[1]],[1e-12,1e-1],'--k')
ax[1,1].legend(frameon=False)
ax[1,1].set_xlabel("$\\mathrm{cpd}$",fontsize=14)
ax[1,1].set_title("Far-field: $\\eta$",fontsize=16)
ax[1,1].set_ylabel("$\\mathrm{m^2 \ cpd^{-1}}$",fontsize=14)
ax[1,1].set_xlim([spec.freqs[0],1])
ax[1,1].set_ylim([1e-12,1e-1])
ax[1,1].text(BA.fcor[1]+.2e-1,5e-12,"$f$",fontsize=14)
ax[1,1].tick_params(axis="both",labelsize=14)
ax[1,1].text(2e-2,3e-1,"d)",fontsize=16)

fig.tight_layout(pad=1)

fig.savefig(figpath + svname, dpi=300, bbox_inches='tight',orientation="horizontal",pad_inches=.3)
