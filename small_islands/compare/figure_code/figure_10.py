import sys
sys.path.append('../src')

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

from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib
matplotlib.rcParams['mathtext.fontset'] = 'stix'

################# Define Paths ###################

fpath1 = '../model_output/small_islands/r50/' #path to files from UCSD Library Collection

ncname1 = 'C_LR_FS_r50_Ah00_r0_gp09'

ncFname1 = fpath1 + ncname1 + '.nc'
ncFnameG1 = fpath1 + ncname1 + '_grid.nc'

####################################################


################# Define Paths ###################

fpath2 = '../model_output/small_islands/r150/' #path to files from UCSD Library Collection

ncname2 = 'C_LR_FS_r150_Ah00_r0_gp09'

ncFname2 = fpath2 + ncname2 + '.nc'
ncFnameG2 = svpath2 + ncname2 + '_grid.nc'

####################################################

################# Define Paths ###################

fpath3 = '../model_output/small_islands/r150/' #path to files from UCSD Library Collection

ncname3 = 'C_LR_FS_r200_Ah00_r0_gp09'

ncFname3 = fpath3 + ncname3 + '.nc'
ncFnameG3 = fpath3 + ncname3 + '_grid.nc'

####################################################

svpath="../small_islands/compare/figures/" #path to git repository

################# Read Grid ###################
xrng = ()
yrng = ()

g1 = rgm.read_grid(xrng,yrng,ncFnameG1)
g1.xx,g1.yy = np.meshgrid(g1.x,g1.y)
g1.fcor = gfd.coriolis(8)
g1.gp=.09
g1.Ld=np.sqrt(g1.gp*g1.H) / g1.fcor[0]
g1.epsilon_L = (g1.Xrad**2) / (g1.Ld**2)
#################################################

################# Read Grid ###################
xrng = ()
yrng = ()

g2 = rgm.read_grid(xrng,yrng,ncFnameG2)
g2.xx,g2.yy = np.meshgrid(g2.x,g2.y)
g2.fcor = gfd.coriolis(8)
g2.gp=.09
g2.Ld=np.sqrt(g2.gp*g2.H) / g2.fcor[0]
g2.epsilon_L = (g2.Xrad**2) / (g2.Ld**2)
#################################################

################# Read Grid ###################
xrng = ()
yrng = ()

g3 = rgm.read_grid(xrng,yrng,ncFnameG3)
g3.xx,g3.yy = np.meshgrid(g3.x,g3.y)
g3.fcor = gfd.coriolis(8)
g3.gp=.09
g3.Ld=np.sqrt(g3.gp*g3.H) / g3.fcor[0]
g3.epsilon_L = (g3.Xrad**2) / (g3.Ld**2)
#################################################

jt = 86
trng=(jt,jt+1)
mo50 = rgm.read_output_ts(trng,xrng,yrng,ncFname1)
mo150 = rgm.read_output_ts(trng,xrng,yrng,ncFname2)
mo200 = rgm.read_output_ts(trng,xrng,yrng,ncFname3)



Ny = len(g1.y)
Nx = len(g1.x)
ind = int(np.floor(Ny/2))
scale = 1
nskip = 30
hcblim = (-.005,.005)
vellim = (-.2,.2)
rolim=(-8,8)


Ny = len(g1.y)
Nx = len(g1.x)
indy = int(np.floor(Ny/2))
indx = int(np.floor(Nx/2))


xmath = g1.x - g1.x[indx]
ymath = g1.y - g1.y[indy]
Xmath,Ymath = np.meshgrid(xmath,ymath)
Zmath = Xmath + 1j * Ymath
r = np.abs(Zmath)
th = np.angle(Zmath)

ticks = (-.005,0,.005)

ss=np.s_[::nskip]
xk = 0.9
yk = 1.02
xpc,ypc = tls.pcm_xy(xmath,ymath)

width = .005
color = "dodgerblue"
svname = "FIGURE_10.png"
fig, ax = plt.subplots(nrows=1,ncols=3,figsize=(14,5))


ax[0].set_xlim(-1000,1000)
ax[0].set_ylim(-1000,1000)
ax[0].set_aspect("equal")
ax[0].set_ylabel("Y [km]",fontsize=16,labelpad=-10)
ax[0].set_xlabel("X [km]",fontsize=16)
ax[0].tick_params(axis="both",labelsize=14)
ccm0 = cmocean.cm.curl
ccm0.set_bad(color='gray')
p0=ax[0].pcolormesh(xpc/1000,ypc/1000,np.ma.array(mo50.eta[0,:,:],mask=(g1.mask_h==0)),vmin=hcblim[0],vmax=hcblim[-1],cmap=ccm0)

ax[0].set_title("$\\epsilon_L = $ %.2f " % (g1.epsilon_L),fontsize=20)
q0 = ax[0].quiver(xmath[ss]/1000,ymath[ss]/1000,
               np.ma.array(mo50.uu[0,ss,ss],mask=(g1.mask_u[ss,ss]==0)),
               np.ma.array(mo50.vv[0,ss,ss],mask=(g1.mask_v[ss,ss]==0)),scale=scale,color=color,width=width)
qk0 = ax[0].quiverkey(q0,xk,yk,.2, '$0.2 \ m s^{-1}$',coordinates = 'axes',fontproperties={'size':16})
ax[0].text(-1500,1300,"a)",fontsize=16)


ax[1].set_xlim(-1000,1000)
ax[1].set_ylim(-1000,1000)
ax[1].set_aspect("equal")
ax[1].set_ylabel("Y [km]",fontsize=16,labelpad=-10)
ax[1].set_xlabel("X [km]",fontsize=16)
ax[1].tick_params(axis="both",labelsize=14)
ax[1].text(-1500,1300,"b)",fontsize=16)
ccm0 = cmocean.cm.curl
ccm0.set_bad(color='gray')
p0=ax[1].pcolormesh(xpc/1000,ypc/1000,np.ma.array(mo150.eta[0,:,:],mask=(g2.mask_h==0)),vmin=hcblim[0],vmax=hcblim[-1],cmap=ccm0)

ax[1].set_title("$\\epsilon_L = $ %.2f " % (g2.epsilon_L),fontsize=20)
q0 = ax[1].quiver(xmath[ss]/1000,ymath[ss]/1000,
               np.ma.array(mo150.uu[0,ss,ss],mask=(g2.mask_u[ss,ss]==0)),
               np.ma.array(mo150.vv[0,ss,ss],mask=(g2.mask_v[ss,ss]==0)),scale=scale,color=color,width=width)
qk0 = ax[1].quiverkey(q0,xk,yk,.2, '$0.2 \ m s^{-1}$',coordinates = 'axes',fontproperties={'size':16})

plt.subplots_adjust(wspace=.3)
cax= inset_axes(ax[1],
                   width="100%",  # width = 5% of parent_bbox width
                   height="5%",  # height : 50%
                   loc='lower left',
                   bbox_to_anchor=(0, -.28, 1, 1),
                   bbox_transform=ax[1].transAxes,
                   borderpad=0,
                   )

cb=plt.colorbar(p0,cax=cax,orientation="horizontal",extend="both")
cb.set_label('$\eta$ [m]', fontsize=16)
cb.ax.tick_params(labelsize=14)
cb.set_ticks(ticks)

ax[2].set_xlim(-1000,1000)
ax[2].set_ylim(-1000,1000)
ax[2].set_aspect("equal")
ax[2].set_ylabel("Y [km]",fontsize=16,labelpad=-10)
ax[2].set_xlabel("X [km]",fontsize=16)
ax[2].tick_params(axis="both",labelsize=14)
ax[2].text(-1500,1300,"c)",fontsize=16)

ccm0 = cmocean.cm.curl
ccm0.set_bad(color='gray')
p0=ax[2].pcolormesh(xpc/1000,ypc/1000,np.ma.array(mo200.eta[0,:,:],mask=(g3.mask_h==0)),vmin=hcblim[0],vmax=hcblim[-1],cmap=ccm0)
q0 = ax[2].quiver(xmath[ss]/1000,ymath[ss]/1000,
               np.ma.array(mo200.uu[0,ss,ss],mask=(g3.mask_u[ss,ss]==0)),
               np.ma.array(mo200.vv[0,ss,ss],mask=(g3.mask_v[ss,ss]==0)),scale=scale,color=color,width=width)
qk0 = ax[2].quiverkey(q0,xk,yk,.2, '$0.2 \ m s^{-1}$',coordinates = 'axes',fontproperties={'size':16})
ax[2].set_title("$\\epsilon_L = $ %.2f " % (g3.epsilon_L),fontsize=20)
fig.suptitle("One Inertial Period" % (mo50.tdays*g1.fcor[1]),fontsize=20)

fig.savefig(figpath + svname, dpi=300, bbox_inches='tight',orientation="horizontal",pad_inches=0.3)
