
import sys
sys.path.append('../src')

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
from pycurrents.plot.mpltools import axes_inches as axi

import matplotlib
matplotlib.rcParams['mathtext.fontset'] = 'stix'

"""
"""
################# Define Paths ###################

fpath = '../model_output/remote_wind_forcing/r50/'

ncname = 'C_R5_FS_r50_Ah00_r0_gp09'

ncFname = fpath + ncname + '.nc'
ncFnameG = fpath + ncname + '_grid.nc'


####################################################

figpath = "../remote_wind_forcing/figures/"

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

################# Read Output ###################

mo = rgm.read_output(xrng,yrng,ncFname)
mo.tinertial = mo.tdays * g.fcor[1]
#################################################

######## Find closest to point to make TS ##########
mindist_n = rgm.min_dist(g.xx,g.yy,4000e3,3050e3)
mindist_s = rgm.min_dist(g.xx,g.yy,4000e3,2950e3)
mindist_e = rgm.min_dist(g.xx,g.yy,4050e3,3000e3)
mindist_w = rgm.min_dist(g.xx,g.yy,3950e3,3000e3)
mindist_ff = rgm.min_dist(g.xx,g.yy,7700e3,5700e3)
mindist_ff2 = rgm.min_dist(g.xx,g.yy,4000e3,5500e3)

ele_n = [mindist_n[0],mindist_n[1]] #y,x
ele_e = [mindist_e[0],mindist_e[1]] #y,x
ele_s = [mindist_s[0],mindist_s[1]] #y,x
ele_w = [mindist_w[0],mindist_w[1]] #y,x
ele_ff = [mindist_ff[0],mindist_ff[1]] #y,x
ele_ff2 = [mindist_ff2[0],mindist_ff2[1]] #y,x

###################################################
Ny = len(g.y)
Nx = len(g.x)


inmid_y = int(np.floor(Ny / 2))
inmid_x = int(np.floor(Nx / 2))

hhtran = mo.eta[:,inmid_y,:]
uutran = mo.uu[:,inmid_y,:]
vvtran = mo.vv[:,inmid_y,:]


indy = int(np.floor(Ny/2))
indx = int(np.floor(Nx/2))

xmath = g.x - g.x[indx]
ymath = g.y - g.y[indy]
Xmath,Ymath = np.meshgrid(xmath,ymath)

scale = 1
nskip = 20
nr1=0
nr2 = 481
hcblim = (-.005,.005)
vellim = (-.2,.2)

ss=np.s_[::nskip]
xk = 0.9
yk = 1.03
xpc,ypc = tls.pcm_xy(xmath,ymath)

width = .005
color = "dodgerblue"


fig = plt.figure(figsize=(7,12))


ax0 = axi(fig,[1,8.5,5,3])
ax2 = axi(fig,[1,1,5,3],sharex = ax0)
ax1 = axi(fig,[1,5.5,5,3],sharex = ax0)

axcb0 = axi(fig,[1,.025,5,.3])
axcb1 = axi(fig,[1,5,5,.3])

ax0.set_ylabel("Y [km]",fontsize=16,labelpad=-10)
ax0.set_xlabel("X [km]",fontsize=16)
ax0.tick_params(axis="both",labelsize=14)
ccm0 = cmocean.cm.curl
ccm0.set_bad(color='gray')
p0=ax0.pcolormesh(xpc/1000,ypc/1000,np.ma.array(mo.eta[302,:,:],mask=(g.mask_h==0)),vmin=hcblim[0],vmax=hcblim[-1],cmap=ccm0)
fig.suptitle("$\\epsilon_L = $ %.2f " % (g.epsilon_L),fontsize=20)
q0 = ax0.quiver(xmath[ss]/1000,ymath[ss]/1000,
               np.ma.array(mo.uu[302,ss,ss],mask=(g.mask_u[ss,ss]==0)),
               np.ma.array(mo.vv[302,ss,ss],mask=(g.mask_v[ss,ss]==0)),scale=scale,color=color,width=width)
qk0 = ax0.quiverkey(q0,xk,yk,.1, '$0.1 \ \\mathrm{ms^{-1}}$',coordinates = 'axes',fontproperties={'size':16})
ax0.set_xlim([-2000,3000])
ax0.set_ylim([-1000,1000])
ax0.set_aspect("equal")
ax0.set_title("%.2f Inertial Periods" % (mo.tdays[302]*g.fcor[1]),fontsize=20)
ax0.text(-2200,1250,"a)",fontsize=20)
ax1.set_ylabel("Y [km]",fontsize=16,labelpad=-10)
ax1.set_xlabel("X [km]",fontsize=16)
ax1.tick_params(axis="both",labelsize=14)
ccm0 = cmocean.cm.curl
ccm0.set_bad(color='gray')
pc0=ax1.pcolormesh(xpc/1000,ypc/1000,np.ma.array(mo.eta[323,:,:],mask=(g.mask_h==0)),vmin=hcblim[0],vmax=hcblim[-1],cmap=ccm0)
q0 = ax1.quiver(xmath[ss]/1000,ymath[ss]/1000,
               np.ma.array(mo.uu[323,ss,ss],mask=(g.mask_u[ss,ss]==0)),
               np.ma.array(mo.vv[323,ss,ss],mask=(g.mask_v[ss,ss]==0)),scale=scale,color=color,width=width)
qk0 = ax1.quiverkey(q0,xk,yk,.1, '$0.1 \ \\mathrm{ms^{-1}}$',coordinates = 'axes',fontproperties={'size':16})
ax1.set_xlim([-2000,3000])
ax1.set_ylim([-1000,1000])
ax1.set_aspect("equal")
ax1.set_title("%.2f Inertial Periods" % (mo.tdays[323]*g.fcor[1]),fontsize=20)
ax1.text(-2200,1250,"b)",fontsize=20)

cb1=fig.colorbar(pc0,cax=axcb1,orientation="horizontal")
cb1.ax.tick_params(axis="both",labelsize=13)
cb1.set_label('$\\eta \ [\\mathrm{m}]$',fontsize=13)

pc1=ax2.pcolormesh(xmath/1000,mo.tinertial,hhtran,cmap=cmocean.cm.curl,vmin=-.005,vmax=.005)
ax2.set_xlabel("X [KM]",fontsize=14)
ax2.set_xlim([-2000,3000])
ax2.set_ylabel("Inertial Periods",fontsize=14)
ax2.fill_between((-50,50), (0,0),(mo.tinertial.max(),mo.tinertial.max()),facecolor='black')
ax2.set_title("$\\eta$",fontsize=16)
ax2.plot([xmath[0]/1000,xmath[-1]/1000],[1,1],'--k')
ax2.tick_params(axis="both",labelsize=12)
ax2.get_shared_x_axes().join(ax1, ax2)
ax2.fill_between((-50,50), (0,0),(mo.tinertial.max(),mo.tinertial.max()),facecolor='black')
ax2.fill_between((1000,2000), (0,0),(mo.tinertial.max(),mo.tinertial.max()),facecolor='cyan', alpha=0.3)
ax2.set_title("$\\eta$",fontsize=24)
ax2.plot([xmath[0]/1000,xmath[-1]/1000],[1,1],'--k')
ax2.text(-1900,1.06,"Wind off",fontsize=14)
ax2.text(1000,.05,"Wind forcing  \nregion",fontsize=8)

cb1=fig.colorbar(pc1,cax=axcb0,orientation="horizontal")
cb1.ax.tick_params(axis="both",labelsize=13)
cb1.set_label('$\\eta \ [\\mathrm{m}]$',fontsize=13)
ax2.text(-2200,6,"c)",fontsize=20)

fig.savefig(figpath + 'FIGURE_18' + ".png", dpi=300, bbox_inches='tight',orientation="horizontal")
