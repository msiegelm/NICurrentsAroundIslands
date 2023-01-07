
import sys
apath = "/Users/miks/Desktop/Palau/reduced_gravity_model/rg_tools/"
sys.path.append(apath)
import rg_plotting as rgp
import animate_plotting as ap
import a_tools as at
from ms_toolbox import tools as tls

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

"""
"""
################# Define Paths ###################

fpath = '/Users/miks/Desktop/Palau/reduced_gravity_model/Runs/circle_PW/base/'
svpath = fpath + 'output/'

ncname = 'C_LR_FS_r50_Ah00_r0_gp09'

ncFname = svpath + ncname + '.nc'
ncFnameG = svpath + ncname + '_grid.nc'

figpath = fpath + 'files_for_publication/figures/'   #path to save figures to


####################################################



################# Read Grid ###################
xrng = ()
yrng = ()

g = rgm.read_grid(xrng,yrng,ncFnameG)
g.xx,g.yy = np.meshgrid(g.x,g.y)
g.fcor = gfd.coriolis(8)
g.gp=.09
g.Ld=np.sqrt(g.gp*g.H) / g.fcor[0]
g.epsilon_L = (g.Xrad**2) / (g.Ld**2)
#################################################

######## Check to confirm output slice ##########
fig,ax = plt.subplots()
pc = ax.pcolormesh(g.el)
ax.set_aspect("equal")
#################################################

################# Read Output ###################

mo = rgm.read_output(xrng,yrng,ncFname)
mo.tinertial = mo.tdays * g.fcor[1]
#################################################




Ny = len(g.y)
Nx = len(g.x)
indy = int(np.floor(Ny/2))
indx = int(np.floor(Nx/2))

xmath = g.x - g.x[indx]
ymath = g.y - g.y[indy]

xmath_g,ymath_g = np.meshgrid(xmath,ymath)

inmid_y = int(np.floor(Ny / 2))
inmid_x = int(np.floor(Nx / 2))

hhtran = mo.eta[:431,inmid_y,:]
uutran = mo.uu[:431,inmid_y,:]
vvtran = mo.vv[:431,inmid_y,:]


(nt,nx) = np.shape(hhtran)

dtdays=mo.tdays[1]-mo.tdays[0]

ht = hhtran[:,10]
nsmooth=3
spec = spectra.spectrum(ht,dt=dtdays, nfft=None, window='quadratic')

hhspec_x = np.ones((len(spec.freqs),nx))
uuspec_x = np.ones((len(spec.freqs),nx))
vvspec_x = np.ones((len(spec.freqs),nx))
# nsmooth=3

for jx in range(nx):
    ht = hhtran[:,jx]
    spec = spectra.spectrum(ht,dt=dtdays, nfft=None, window='quadratic')
    hhspec_x[:,jx] = spec.psd

    ht = uutran[:,jx]
    spec = spectra.spectrum(ht,dt=dtdays, nfft=None,window='quadratic')
    uuspec_x[:,jx] = spec.psd

    ht = vvtran[:,jx]
    spec = spectra.spectrum(ht,dt=dtdays, nfft=None,window='quadratic')
    vvspec_x[:,jx] = spec.psd

    if jx == 1:
        freqx = spec.freqs
    print(jx)

import matplotlib
matplotlib.rcParams['mathtext.fontset'] = 'stix'

cx,cy=np.meshgrid(xmath/1000,mo.tinertial[:431])
fig,ax=plt.subplots(ncols=2,nrows=3,figsize=(13,10))
fig.tight_layout(pad=7.3,w_pad=8)
pc0=ax[0,0].pcolormesh(xmath/1000,mo.tinertial[:432],hhtran,cmap=cmocean.cm.diff,vmin=-.005,vmax=.005)
# ax[0,0].contour(cx,cy,hhtran,levels=[-.004,.004],colors="c",linestyles="dashed")
ax[0,0].set_xlabel("X [$\\mathrm{km}$]",fontsize=14)
ax[0,0].set_ylabel("Inertial Periods",fontsize=14)
ax[0,0].fill_between((-50,50), (0,0),(mo.tinertial.max(),mo.tinertial.max()),facecolor='black')
ax[0,0].set_title("$\\eta$",fontsize=16)
ax[0,0].plot([xmath[0]/1000,xmath[-1]/1000],[1,1],'--m')
ax[0,0].text(-2000,.3,"Wind off",fontsize=14,color="m")
ax[0,0].text(-2500,5.5,"a)",fontsize=14)
ax[0,0].set_ylim([mo.tinertial[0],mo.tinertial[431]])
ax[0,0].tick_params(axis="both",labelsize=12)
cb1=fig.colorbar(pc0,ax=ax[0,0],extend="both")
cb1.ax.tick_params(axis="both",labelsize=13)
# cb3.ax.xaxis.offsetText.set_fontsize(16)
cb1.set_label('$\\eta [\\mathrm{m}]$',fontsize=13)
ax[0,0].set_xlim([-2000,2000])

pc1=ax[1,0].pcolormesh(xmath/1000,mo.tinertial[:432],uutran,cmap=cmocean.cm.diff,vmin=-.07,vmax=.07)
# ax[1,0].contour(cx,cy,uutran,levels=[-.06,.06],colors="c",linestyles="dashed")
ax[1,0].text(-2500,5.5,"c)",fontsize=14)
ax[1,0].set_xlabel("X [$\\mathrm{km}$]",fontsize=14)
ax[1,0].set_ylabel("Inertial Periods",fontsize=14)
ax[1,0].fill_between((-50,50), (0,0),(mo.tinertial.max(),mo.tinertial.max()),facecolor='black')
ax[1,0].set_title("$u$",fontsize=16)
ax[1,0].plot([xmath[0]/1000,xmath[-1]/1000],[1,1],'--m')
ax[1,0].text(-2000,.3,"Wind off",fontsize=14,color="m")
ax[1,0].set_ylim([mo.tinertial[0],mo.tinertial[431]])
ax[1,0].tick_params(axis="both",labelsize=12)
ax[1,0].set_xlim([-2000,2000])
cb2=fig.colorbar(pc1,ax=ax[1,0],extend="both")
cb2.ax.tick_params(axis="both",labelsize=13)
# cb3.ax.xaxis.offsetText.set_fontsize(16)
cb2.set_label('$u \ [\\mathrm{ms^{-1}]}$',fontsize=13)

pc2=ax[2,0].pcolormesh(xmath/1000,mo.tinertial[:432],vvtran,cmap=cmocean.cm.diff,vmin=-.07,vmax=.07)
# ax[2,0].contour(cx,cy,vvtran,levels=[-.06,.06],colors="c",linestyles="dashed")
ax[2,0].text(-2500,5.5,"e)",fontsize=14)
ax[2,0].set_xlabel("X [$\\mathrm{km}$]",fontsize=14)
ax[2,0].set_ylabel("Inertial Periods",fontsize=14)
ax[2,0].fill_between((-50,50), (0,0),(mo.tinertial.max(),mo.tinertial.max()),facecolor='black')
ax[2,0].set_title("$v$",fontsize=16)
ax[2,0].plot([xmath[0]/1000,xmath[-1]/1000],[1,1],'--m')
ax[2,0].text(-2000,.3,"Wind off",fontsize=14,color="m")
ax[2,0].set_ylim([mo.tinertial[0],mo.tinertial[431]])
ax[2,0].tick_params(axis="both",labelsize=12)
cb3=fig.colorbar(pc2,ax=ax[2,0],extend="both")
cb3.ax.tick_params(axis="both",labelsize=13)
# cb3.ax.xaxis.offsetText.set_fontsize(16)
cb3.set_label('$v \ [\\mathrm{ms^{-1}}]$',fontsize=13)
ax[2,0].set_xlim([-2000,2000])


xpc,ypc = tls.pcm_xy(xmath/1000,freqx)

pc0=ax[0,1].pcolormesh(xpc,np.log10(ypc),np.log10(hhspec_x),vmin=-10,vmax=-5,cmap=cmocean.cm.matter)
ax[0,1].plot((-3000,3000),(np.log10(g.fcor[1]),np.log10(g.fcor[1])),'--',color="deepskyblue")
ax[0,1].text(-2500,.15,"b)",fontsize=14)
cb0=fig.colorbar(pc0,ax=ax[0,1],extend="both")
cb0.set_label("$log_{10}(\\mathrm{m^2 cpd^{-1}})$",fontsize=13)
cb0.ax.tick_params(axis="both",labelsize=13)
ax[0,1].set_ylabel("$log_{10}(\\mathrm{cpd})$",fontsize=14)
ax[0,1].set_xlabel("X [$\\mathrm{km}$]",fontsize=14)
ax[0,1].set_title("Power Spectral Density: $\\eta$",fontsize=16)
ax[0,1].set_ylim([-1,0])
ax[0,1].set_xlim([-2000,2000])
ax[0,1].text(-2000 + 150,np.log10(g.fcor[1])  + .06,"$f$",fontsize=16,color="deepskyblue")
ax[0,1].fill_between((-50,50), (-2.4,-2.4),(0,0),facecolor='black')
ax[0,1].tick_params(axis="both",labelsize=12)

pc1=ax[1,1].pcolormesh(xpc,np.log10(ypc),np.log10(uuspec_x),vmin=-4,vmax=-1,cmap=cmocean.cm.matter)
ax[1,1].text(-2500,.15,"d)",fontsize=14)
ax[1,1].plot((-3000,3000),(np.log10(g.fcor[1]),np.log10(g.fcor[1])),'--',color="deepskyblue")
ax[1,1].text(-2000 + 150,np.log10(g.fcor[1]) + .06,"$f$",fontsize=16,color="deepskyblue")
cb1=fig.colorbar(pc1,ax=ax[1,1],extend="both")
cb1.set_label("$log_{10}(\\mathrm{m^2 s^{-2} cpd^{-1}})$",fontsize=13)
cb1.ax.tick_params(axis="both",labelsize=13)
ax[1,1].set_ylabel("$log_{10}(\\mathrm{cpd})$",fontsize=14)
ax[1,1].set_xlabel("X [$\\mathrm{km}$]",fontsize=14)
ax[1,1].set_title("Power Spectral Density: $u$",fontsize=16)
ax[1,1].set_ylim([-1,0])
ax[1,1].set_xlim([-2000,2000])
ax[1,1].fill_between((-50,50), (-2.4,-2.4),(0,0),facecolor='black')
ax[1,1].tick_params(axis="both",labelsize=12)

pc2=ax[2,1].pcolormesh(xpc,np.log10(ypc),np.log10(vvspec_x),vmin=-4,vmax=-1,cmap=cmocean.cm.matter)
ax[2,1].text(-2500,.15,"f)",fontsize=14)
ax[2,1].plot((-3000,3000),(np.log10(g.fcor[1]),np.log10(g.fcor[1])),'--',color="deepskyblue")
ax[2,1].text(-2000 + 150,np.log10(g.fcor[1])  + .06,"$f$",fontsize=16,color="deepskyblue")
cb2=fig.colorbar(pc2,ax=ax[2,1],extend="both")
cb2.set_label("$log_{10}(\\mathrm{m^2 s^{-2} cpd^{-1}})$",fontsize=13)
cb2.ax.tick_params(axis="both",labelsize=13)
ax[2,1].set_ylabel("$log_{10}(\\mathrm{cpd})$",fontsize=14)
ax[2,1].set_xlabel("X [$\\mathrm{km}$]",fontsize=14)
ax[2,1].set_title("Power Spectral Density: $v$",fontsize=16)
ax[2,1].set_ylim([-1,0])
ax[2,1].set_xlim([-2000,2000])
ax[2,1].fill_between((-50,50), (-2.4,-2.4),(0,0),facecolor='black')
ax[2,1].tick_params(axis="both",labelsize=12)


fig.savefig(figpath + 'FIGURE_09' + ".png", dpi=300, bbox_inches='tight',orientation="horizontal")
