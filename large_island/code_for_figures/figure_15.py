
import sys
sys.path.append('../src')

from ms_toolbox import tools as tls

import matplotlib.pyplot as plt
import numpy as np

from ms_toolbox import gfd
from ms_toolbox import rgmodel as rgm
from pycurrents.num import spectra

from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy.interpolate import interp2d as int2d

from netCDF4 import Dataset
from pycurrents.system import Bunch
from matplotlib.ticker import FormatStrFormatter
import scipy.special as mb

import cmocean
import matplotlib.cm as cm
from pycurrents.num.harmfit import complex_demodulation
import matplotlib

csfont = {'fontname':'Times New Roman'}
matplotlib.rcParams['mathtext.fontset'] = 'stix'

"""
"""

################# Define Paths ###################

anpath = "../large_island/data_sub/r300/"
figpath = '../large_island/figures/'   #path to save figures to

# ########################################


ftw =  0.1118



data = np.load(anpath + "complex_demod_full_x.npz")
tinertial = data["tinertial"]
etaresid = data["etaresid"]
uresid = data["uresid"]
vresid = data["vresid"]
etacdetw = data["etacdetw"]
urectw = data["urectw"]
vcdetw = data["vcdetw"]
xmath = data["xmath"]
tdays = data["tdays"]


ind = np.argmin(np.abs(tdays[:]*ftw-1))


fig,ax=plt.subplots(ncols=2,nrows=3,figsize=(14,10))
fig.tight_layout(pad=5.7)
pc0=ax[0,0].pcolormesh(xmath/1000,tinertial[:],etaresid,cmap=cmocean.cm.diff,vmin=-.005,vmax=.005)
ax[0,0].set_xlabel("X [KM]",fontsize=14)
ax[0,0].set_ylabel("Inertial Periods",fontsize=14)
ax[0,0].fill_between((-300,300), (0,0),(tinertial.max(),tinertial.max()),facecolor='black')
ax[0,0].set_title("Residual: $\\eta$",fontsize=16)
ax[0,0].plot([xmath[0]/1000,xmath[-1]/1000],[1,1],'--k')
ax[0,0].set_xlim([-5000,5000])
ax[0,0].tick_params(axis="both",labelsize=12)
ax[0,0].text(-6000,12,"a)",fontsize=16)
ax[0,0].set_ylim([tinertial[ind],tinertial[-ind]])

cb1=fig.colorbar(pc0,ax=ax[0,0])
cb1.ax.tick_params(axis="both",labelsize=13)
cb1.set_label('$\\eta \ [\\mathrm{m}]$',fontsize=13)

pc1=ax[1,0].pcolormesh(xmath/1000,tinertial[:],uresid,cmap=cmocean.cm.diff,vmin=-.05,vmax=.05)
ax[1,0].set_xlabel("X [KM]",fontsize=14)
ax[1,0].set_ylabel("Inertial Periods",fontsize=14)
ax[1,0].fill_between((-300,300), (0,0),(tinertial.max(),tinertial.max()),facecolor='black')
ax[1,0].set_title("Residual: $u$",fontsize=16)
ax[1,0].plot([xmath[0]/1000,xmath[-1]/1000],[1,1],'--k')
ax[1,0].tick_params(axis="both",labelsize=12)
ax[1,0].set_xlim([-1000,1000])
cb2=fig.colorbar(pc1,ax=ax[1,0])
cb2.ax.tick_params(axis="both",labelsize=13)
cb2.set_label('$u \ [\\mathrm{ms^{-1}}]$',fontsize=13)
ax[1,0].text(-1200,12,"c)",fontsize=16)
ax[1,0].set_ylim([tinertial[ind],tinertial[-ind]])

pc2=ax[2,0].pcolormesh(xmath/1000,tinertial[:],vresid,cmap=cmocean.cm.diff,vmin=-.05,vmax=.05)
ax[2,0].set_xlabel("X [KM]",fontsize=14)
ax[2,0].set_ylabel("Inertial Periods",fontsize=14)
ax[2,0].fill_between((-300,300), (0,0),(tinertial.max(),tinertial.max()),facecolor='black')
ax[2,0].set_title("Residual: $v$",fontsize=16)
ax[2,0].plot([xmath[0]/1000,xmath[-1]/1000],[1,1],'--k')
ax[2,0].tick_params(axis="both",labelsize=12)
ax[2,0].set_xlim([-1000,1000])
ax[2,0].text(-1200,12,"e)",fontsize=16)
cb3=fig.colorbar(pc2,ax=ax[2,0])
cb3.ax.tick_params(axis="both",labelsize=13)
cb3.set_label('$v \ [\\mathrm{ms^{-1}}]$',fontsize=13)
ax[2,0].set_ylim([tinertial[ind],tinertial[-ind]])

pc0=ax[0,1].pcolormesh(xmath/1000,tdays[:]*ftw,etacdetw,cmap=cmocean.cm.diff,vmin=-.005,vmax=.005)
ax[0,1].set_xlabel("X [KM]",fontsize=14)
ax[0,1].set_ylabel("Trapped Wave Periods",fontsize=14)
ax[0,1].fill_between((-300,300), (0,0),(tinertial.max(),tinertial.max()),facecolor='black')
ax[0,1].set_title("Trapped Wave: $\\eta$",fontsize=16)
ax[0,1].plot([xmath[0]/1000,xmath[-1]/1000],[.4,.4],'--k')
ax[0,1].set_xlim([-1000,1000])
ax[0,1].tick_params(axis="both",labelsize=12)
ax[0,1].text(-1200,5,"b)",fontsize=16)
ax[0,1].set_ylim([1,np.max(tdays[:]*ftw)-1])
cb1=fig.colorbar(pc0,ax=ax[0,1])
cb1.ax.tick_params(axis="both",labelsize=13)
cb1.set_label('$\\eta \ [\\mathrm{m}]$',fontsize=13)

pc1=ax[1,1].pcolormesh(xmath/1000,tdays[:]*ftw,urectw,cmap=cmocean.cm.diff,vmin=-.05,vmax=.05)
ax[1,1].set_xlabel("X [KM]",fontsize=14)
ax[1,1].set_ylabel("Trapped Wave Periods",fontsize=14)
ax[1,1].fill_between((-300,300), (0,0),(tinertial.max(),tinertial.max()),facecolor='black')
ax[1,1].set_title("Trapped Wave: $u$",fontsize=16)
ax[1,1].plot([xmath[0]/1000,xmath[-1]/1000],[.4,.4],'--k')
ax[1,1].tick_params(axis="both",labelsize=12)
ax[1,1].set_xlim([-1000,1000])
ax[1,1].text(-1200,5,"d)",fontsize=16)
ax[1,1].set_ylim([1,np.max(tdays[:]*ftw)-1])

cb2=fig.colorbar(pc1,ax=ax[1,1])
cb2.ax.tick_params(axis="both",labelsize=13)
cb2.set_label('$u \ [\\mathrm{ms^{-1}}]$',fontsize=13)

pc2=ax[2,1].pcolormesh(xmath/1000,tdays[:]*ftw,vcdetw,cmap=cmocean.cm.diff,vmin=-.05,vmax=.05)
ax[2,1].set_xlabel("X [KM]",fontsize=14)
ax[2,1].set_ylabel("Trapped Wave Periods",fontsize=14)
ax[2,1].fill_between((-300,300), (0,0),(tinertial.max(),tinertial.max()),facecolor='black')
ax[2,1].set_title("Trapped Wave: $v$",fontsize=16)
ax[2,1].plot([xmath[0]/1000,xmath[-1]/1000],[.4,.4],'--k')
ax[2,1].tick_params(axis="both",labelsize=12)
ax[2,1].set_xlim([-1000,1000])
ax[2,1].text(-1200,5,"f)",fontsize=16)
ax[2,1].set_ylim([1,np.max(tdays[:]*ftw)-1])

cb3=fig.colorbar(pc2,ax=ax[2,1])
cb3.ax.tick_params(axis="both",labelsize=13)
cb3.set_label('$v \ [\\mathrm{ms^{-1}}]$',fontsize=13)

fig.savefig(figpath + 'FIGURE_15' + ".png", dpi=300, bbox_inches='tight',orientation="horizontal")
