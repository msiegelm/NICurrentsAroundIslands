
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
    b.vresidE = aN["vresidE"]
    b.vtwE = aN["vtwE"]

    b.t = aN["t"][:]

    fcor = gfd.coriolis(8)
    b.fcor = fcor

    return b

################# Read Grid ###################

file_name1 = "complex_demod_ke_x.npz"

BA = load_var(p1,file_name1)

r75 = load_var(p2,file_name1)
ftw =  0.1118 # LH69 solution

ind = np.argmin(np.abs(BA.t[:]*ftw-1))

##################################
svname="FIGURE_17.png"
fig,ax = plt.subplots(ncols=1,nrows=3,figsize=(10,9),sharex=True)
ax[0].plot(BA.t*BA.fcor[1],.5*(BA.uff2**2+BA.vff2**2),'m',label="1 Inertial Period Forcing")
ax[0].plot(BA.t*BA.fcor[1],.5*(r75.uff2**2+r75.vff2**2),'b',label="3 Inertial Period Forcing")
ax[0].arrow(3.57,.0013,0,.0005,color="m",width = 0.05,head_width=.2,head_length=.0001)
ax[0].arrow(4.8,.0013,0,.0005,color="b",width = 0.05,head_width=.2,head_length=.0001)

ax[0].legend(loc="lower right",fontsize=13)
ax[0].set_title("Kinetic Energy: 3000 km from eastern island boundary",fontsize=16)
ax[0].set_xlim([0,11])
ax[0].set_ylim([0,.0035])
ax[0].set_ylabel("$\\mathrm{m^2s^{-2}}$",fontsize=16)
ax[0].tick_params("both",labelsize=16)
ax[0].text(-1.5,.0035,"a)",fontsize=16)


ax[1].plot(BA.t*BA.fcor[1],BA.vresidE,'m',label="1 Inertial Period Forcing")
ax[1].plot(BA.t*BA.fcor[1],r75.vresidE,'b',label="3 Inertial Period Forcing")
ax[1].fill_between((BA.t[0]*BA.fcor[1],BA.t[ind]*BA.fcor[1]),(-.05,-.05),(.05,.05),facecolor="k",alpha=.5)
ax[1].fill_between((BA.t[-1]*BA.fcor[1],BA.t[-ind]*BA.fcor[1]),(-.05,-.05),(.05,.05),facecolor="k",alpha=.5)
ax[1].legend(loc="upper right",fontsize=13)
ax[1].set_title("Eastern Boundary: Residual $v$",fontsize=16)
ax[1].set_ylabel("$\\mathrm{ms^{-1}}$",fontsize=16)
ax[1].tick_params("both",labelsize=16)
ax[1].text(-1.5,.052,"b)",fontsize=16)
ax[1].set_xlim([0,BA.t[1035]*BA.fcor[1]])
ax[1].set_ylim([-.05,.05])
ftw =  0.1118
ax[1].set_xlim([0,11])



def in2tw(t):
    return (t/BA.fcor[1])*ftw


def tw2in(t):
    return (t/ftw)*BA.fcor[1]


ax[2].plot(BA.t*BA.fcor[1],BA.vtwE,'m',label="1 Inertial Period Forcing")
ax[2].plot(BA.t*BA.fcor[1],r75.vtwE,'b',label="3 Inertial Period Forcing")
ax[2].fill_between((BA.t[0]*BA.fcor[1],BA.t[ind]*BA.fcor[1]),(-.1,-.1),(.1,.1),facecolor="k",alpha=.5)
ax[2].fill_between((BA.t[-1]*BA.fcor[1],BA.t[-ind]*BA.fcor[1]),(-.1,-.1),(.1,.1),facecolor="k",alpha=.5)

ax[2].set_title("Eastern Boundary: Trapped Wave $v$",fontsize=16)
ax[2].set_ylabel("$\\mathrm{ms^{-1}}$",fontsize=16)
ax[2].set_xlabel("Inertial Periods",fontsize=16)
ax[2].tick_params("both",labelsize=16)
ax[2].text(-1.5,.12,"c)",fontsize=16)
ax[2].set_ylim([-.1,.1])
secax = ax[2].secondary_xaxis(-.4, functions=(in2tw,tw2in))
ax[2].set_xlim([0,BA.t[1035]*BA.fcor[1]])
secax.set_xlabel("Trapped Wave Periods",fontsize=16)
secax.tick_params(labelsize=16)
ax[2].set_xlim([0,11])

fig.tight_layout(pad=1)

fig.savefig(figpath + svname, dpi=300, bbox_inches='tight',orientation="horizontal",pad_inches=.3)
