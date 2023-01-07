


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

from ms_toolbox import tools as tls

ff = '../small_islands/compare/data_subsection/'
p1 = ff + 'r50/'
p2 = ff + 'r75/'
p3 = ff + 'r100/'
p4 = ff + 'r150/'
p5 = ff + 'r200/'

figpath = '../small_islands/compare/figures/'


def load_var(bpath,n_filename):


    aN = np.load(bpath + n_filename)



    b = Bunch()
    b.eps = aN["eps"]
    b.Xrad = aN["Xrad"]
    b.r = aN["r"]
    b.tin = aN["tin"]
    b.urpA = aN["urpA"]
    b.uthpA = aN["uthpA"]
    b.fcor = aN["fcor"]
    b.ke = aN["ke"]
    b.pe = aN["pe"]
    return b
fname="EF_ke_pe.npz"

r50 = load_var(p1,fname)
r75 = load_var(p2,fname)
r100 = load_var(p3,fname)
r150 = load_var(p4,fname)
r200 = load_var(p5,fname)


def in2beat1(tin):
    beat = (tin/r50.fcor[1])/4769
    return beat

def beat2in1(beat):
    tin = (beat * 4769) * r50.fcor[1]
    return tin

def in2beat2(tin):
    beat = (tin/r50.fcor[1])/80.8
    return beat

def beat2in2(beat):
    tin = (beat * 80.8) * r50.fcor[1]
    return tin


def in2beat3(tin):
    beat = (tin/r50.fcor[1])/25.12
    return beat

def beat2in3(beat):
    tin = beat * 25.12 * r50.fcor[1]
    return tin

def in2beat4(tin):
    beat = (tin/r50.fcor[1])/11.11
    return beat

def beat2in4(beat):
    tin = beat * 11.11 * r50.fcor[1]
    return tin

def in2beat5(tin):
    beat = (tin/r50.fcor[1])/8.02
    return beat

def beat2in5(beat):
    tin = beat * 8.02 * r50.fcor[1]
    return tin

svname="FIGURE_13.png"
fig,ax=plt.subplots(nrows=5,ncols=3,figsize=(10,10),sharex=True, sharey=False, constrained_layout=True)

axi = ax[0,0]
pc0=axi.pcolormesh(r50.r/1000 - (r50.Xrad/1000),r50.tin,r50.ke/1e2,cmap=cmocean.cm.amp,vmin = 0,vmax=5)
axi.text(1500,.3,"$\\epsilon_L$ = %.02f" % r50.eps,fontsize=18)
axi.set_ylabel("Inertial Periods",fontsize=14)
secax_y = axi.secondary_yaxis(-.2, functions=(in2beat1, beat2in1))
secax_y.set_ylabel('Beat Periods',fontsize=14)
axi.tick_params(axis="both",labelsize=14)
axi.set_xlim([0,np.max(r200.r/1000)-(r200.Xrad/1000)])
axi.set_title("Kinetic Energy",fontsize=16)
axi = ax[1,0]
pc0=axi.pcolormesh(r75.r/1000 - (r75.Xrad/1000),r50.tin,r75.ke/1e2,cmap=cmocean.cm.amp,vmin = 0,vmax=5)
axi.text(1500,.3,"$\\epsilon_L$ = %.02f" % r75.eps,fontsize=18)
axi.set_ylabel("Inertial Periods",fontsize=14)
axi.tick_params(axis="both",labelsize=14)
axi.set_xlim([0,np.max(r200.r/1000)-(r200.Xrad/1000)])
secax_y = axi.secondary_yaxis(-.2, functions=(in2beat2, beat2in2))
secax_y.set_ylabel('Beat Periods',fontsize=14)

axi = ax[2,0]
pc0=axi.pcolormesh(r100.r/1000 - (r100.Xrad/1000),r50.tin,r100.ke/1e2,cmap=cmocean.cm.amp,vmin = 0,vmax=5)
axi.text(1500,.3,"$\\epsilon_L$ = %.02f" % r100.eps,fontsize=18)
axi.set_ylabel("Inertial Periods",fontsize=14)
axi.tick_params(axis="both",labelsize=14)
axi.set_xlim([0,np.max(r200.r/1000)-(r200.Xrad/1000)])
secax_y = axi.secondary_yaxis(-.2, functions=(in2beat3, beat2in3))
secax_y.set_ylabel('Beat Periods',fontsize=14)


axi = ax[3,0]
pc0=axi.pcolormesh(r150.r/1000 - (r150.Xrad/1000),r50.tin,r150.ke/1e2,cmap=cmocean.cm.amp,vmin = 0,vmax=5)
axi.text(1500,.3,"$\\epsilon_L$ = %.02f" % r150.eps,fontsize=18)
axi.set_ylabel("Inertial Periods",fontsize=14)
axi.tick_params(axis="both",labelsize=14)
axi.set_xlim([0,np.max(r200.r/1000)-(r200.Xrad/1000)])
secax_y = axi.secondary_yaxis(-.2, functions=(in2beat4, beat2in4))
secax_y.set_ylabel('Beat Periods',fontsize=14)

axi = ax[4,0]
pc0=axi.pcolormesh(r200.r/1000 - (r200.Xrad/1000),r50.tin,r200.ke/1e2,cmap=cmocean.cm.amp,vmin = 0,vmax=5)
axi.text(1500,.3,"$\\epsilon_L$ = %.02f" % r200.eps,fontsize=18)
axi.set_ylabel("Inertial Periods",fontsize=14)
axi.tick_params(axis="both",labelsize=14)
axi.set_xlim([0,np.max(r200.r/1000)-(r200.Xrad/1000)])
axi.set_xlabel("Distance from Island Boundary [km]",fontsize=10)
secax_y = axi.secondary_yaxis(-.2, functions=(in2beat5, beat2in5))
secax_y.set_ylabel('Beat Periods',fontsize=14)
cax= inset_axes(axi,
                   width="100%",  # width = 5% of parent_bbox width
                   height="5%",  # height : 50%
                   loc='lower left',
                   bbox_to_anchor=(0, -.36, 1, 1),
                   bbox_transform=axi.transAxes,
                   borderpad=0,
                   )

cb1=plt.colorbar(pc0,ax=axi,cax=cax,orientation="horizontal")
cb1.set_label("KE [$10^{2} J m^{-2}$]",fontsize=20)
cb1.ax.tick_params(labelsize=14)
cb1.formatter.set_powerlimits((0, 0))


axi = ax[0,1]
pc0=axi.pcolormesh(r50.r/1000 - (r50.Xrad/1000),r50.tin,r50.pe/1e-1,cmap=cmocean.cm.amp,vmin = 0,vmax=5)
axi.text(1500,.3,"$\\epsilon_L$ = %.02f" % r50.eps,fontsize=18)
axi.tick_params(axis="both",labelsize=14)
axi.set_xlim([0,np.max(r200.r/1000)-(r200.Xrad/1000)])
axi.axes.yaxis.set_ticklabels([])
axi.set_title("Potential Energy",fontsize=16)

axi = ax[1,1]
pc0=axi.pcolormesh(r75.r/1000 - (r75.Xrad/1000),r50.tin,r75.pe/1e-1,cmap=cmocean.cm.amp,vmin = 0,vmax=5)
axi.text(1500,.3,"$\\epsilon_L$ = %.02f" % r75.eps,fontsize=18)
axi.tick_params(axis="both",labelsize=14)
axi.set_xlim([0,np.max(r200.r/1000)-(r200.Xrad/1000)])
axi.axes.yaxis.set_ticklabels([])

axi = ax[2,1]
pc0=axi.pcolormesh(r100.r/1000 - (r100.Xrad/1000),r50.tin,r100.pe/1e-1,cmap=cmocean.cm.amp,vmin = 0,vmax=5)
axi.text(1500,.3,"$\\epsilon_L$ = %.02f" % r100.eps,fontsize=18)
axi.tick_params(axis="both",labelsize=14)
axi.set_xlim([0,np.max(r200.r/1000)-(r200.Xrad/1000)])
axi.axes.yaxis.set_ticklabels([])


axi = ax[3,1]
pc0=axi.pcolormesh(r150.r/1000 - (r150.Xrad/1000),r50.tin,r150.pe/1e-1,cmap=cmocean.cm.amp,vmin = 0,vmax=5)
axi.text(1500,.3,"$\\epsilon_L$ = %.02f" % r150.eps,fontsize=18)
axi.tick_params(axis="both",labelsize=14)
axi.set_xlim([0,np.max(r200.r/1000)-(r200.Xrad/1000)])
axi.axes.yaxis.set_ticklabels([])

axi = ax[4,1]
pc0=axi.pcolormesh(r200.r/1000 - (r200.Xrad/1000),r50.tin,r200.pe/1e-1,cmap=cmocean.cm.amp,vmin = 0,vmax=5)
axi.text(1500,.3,"$\\epsilon_L$ = %.02f" % r200.eps,fontsize=18)
axi.tick_params(axis="both",labelsize=14)
axi.set_xlim([0,np.max(r200.r/1000)-(r200.Xrad/1000)])
axi.axes.yaxis.set_ticklabels([])
axi.set_xlabel("Distance from Island Boundary [km]",fontsize=10)

cax= inset_axes(axi,
                   width="100%",  # width = 5% of parent_bbox width
                   height="5%",  # height : 50%
                   loc='lower left',
                   bbox_to_anchor=(0, -.36, 1, 1),
                   bbox_transform=axi.transAxes,
                   borderpad=0,
                   )

cb1=plt.colorbar(pc0,ax=axi,cax=cax,orientation="horizontal")
cb1.set_label("PE [$10^{-1} J m^{-2}$]",fontsize=20)
cb1.ax.tick_params(labelsize=14)
cb1.formatter.set_powerlimits((0, 0))


axi = ax[0,2]
pc0=axi.pcolormesh(r50.r/1000 - (r50.Xrad/1000),r50.tin,r50.urpA/1e6,cmap=cmocean.cm.diff,vmin = -5,vmax=5)
axi.text(1500,.3,"$\\epsilon_L$ = %.02f" % r50.eps,fontsize=18)
axi.tick_params(axis="both",labelsize=14)
axi.set_xlim([0,np.max(r200.r/1000)-(r200.Xrad/1000)])
axi.axes.yaxis.set_ticklabels([])
axi.set_title("Radial Energy Flux",fontsize=16)

axi = ax[1,2]

pc0=axi.pcolormesh(r75.r/1000 - (r75.Xrad/1000),r50.tin,r75.urpA/1e6,cmap=cmocean.cm.diff,vmin = -5,vmax=5)
axi.text(1500,.3,"$\\epsilon_L$ = %.02f" % r75.eps,fontsize=18)
axi.tick_params(axis="both",labelsize=14)
axi.set_xlim([0,np.max(r200.r/1000)-(r200.Xrad/1000)])
axi.axes.yaxis.set_ticklabels([])


axi = ax[2,2]

pc0=axi.pcolormesh(r100.r/1000 - (r75.Xrad/1000),r50.tin,r100.urpA/1e6,cmap=cmocean.cm.diff,vmin = -5,vmax=5)
axi.text(1500,.3,"$\\epsilon_L$ = %.02f" % r100.eps,fontsize=18)
axi.tick_params(axis="both",labelsize=14)
axi.set_xlim([0,np.max(r200.r/1000)-(r200.Xrad/1000)])
axi.axes.yaxis.set_ticklabels([])

axi = ax[3,2]

pc0=axi.pcolormesh(r150.r/1000 - (r150.Xrad/1000),r50.tin,r150.urpA/1e6,cmap=cmocean.cm.diff,vmin = -5,vmax=5)
axi.text(1500,.3,"$\\epsilon_L$ = %.02f" % r150.eps,fontsize=18)
axi.tick_params(axis="both",labelsize=14)
axi.set_xlim([0,np.max(r200.r/1000)-(r200.Xrad/1000)])
axi.axes.yaxis.set_ticklabels([])

axi = ax[4,2]

pc0=axi.pcolormesh(r200.r/1000 - (r200.Xrad/1000),r50.tin,r200.urpA/1e6,cmap=cmocean.cm.diff,vmin = -5,vmax=5)
axi.text(1500,.3,"$\\epsilon_L$ = %.02f" % r200.eps,fontsize=18)
axi.set_xlabel("Distance from Island Boundary [km]",fontsize=10)
axi.tick_params(axis="both",labelsize=14)
axi.set_xlim([0,np.max(r200.r/1000)-(r200.Xrad/1000)])
axi.axes.yaxis.set_ticklabels([])

cax= inset_axes(axi,
                   width="100%",  # width = 5% of parent_bbox width
                   height="5%",  # height : 50%
                   loc='lower left',
                   bbox_to_anchor=(0, -.36, 1, 1),
                   bbox_transform=axi.transAxes,
                   borderpad=0,
                   )

cb1=plt.colorbar(pc0,ax=axi,cax=cax,orientation="horizontal")
cb1.set_label("$10^{6} W m^{-1}$",fontsize=20)
cb1.ax.tick_params(labelsize=14)
cb1.formatter.set_powerlimits((0, 0))

fig.savefig(figpath + svname, dpi=300, bbox_inches='tight',orientation="horizontal",pad_inches=0.2)
