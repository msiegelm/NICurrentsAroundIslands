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

figpath='../remote_wind_forcing/figures/'

################# Define Paths ###################

fpath1 = '../model_output/remote_wind_forcing/r50/' #path to files from UCSD Library Collection

ncname = 'C_R5_FS_r50_Ah00_r0_gp09'

ncFname = fpath1 + ncname + '.nc'
ncFnameG = fpath1 + ncname + '_grid.nc'


####################################################


################# Define Paths ###################

fpath2 = '../model_output/remote_wind_forcing/r150/' #path to files from UCSD Library Collection

ncname2 = 'C_R5_FS_r150_Ah00_r0_gp09'

ncFname2 = fpath2 + ncname2 + '.nc'
ncFnameG2 = fpath2 + ncname2 + '_grid.nc'

####################################################

################# Define Paths ###################

fpath3 = '../model_output/remote_wind_forcing/r150/' #path to files from UCSD Library Collection

ncname3 = 'C_R5_FS_r200_Ah00_r0_gp09'

ncFname3 = fpath3 + ncname3 + '.nc'
ncFnameG3 = fpath3 + ncname3 + '_grid.nc'

####################################################


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

jt = 397
trng=(jt,jt+1)
mo50 = rgm.read_output_ts(trng,xrng,yrng,ncFname1)
mo150 = rgm.read_output_ts(trng,xrng,yrng,ncFname2)
mo200 = rgm.read_output_ts(trng,xrng,yrng,ncFname3)


Ny = len(g1.y)
Nx = len(g1.x)
indy = int(np.floor(Ny/2))
indx = int(np.floor(Nx/2))


scale = .5
nskip = 35
hcblim = (-.005,.005)
vellim = (-.2,.2)
rolim=(-8,8)

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
svname="FIGURE_19.png"

fig, ax = plt.subplots(nrows=2,ncols=3,figsize=(14,10))
plt.subplots_adjust(wspace=.3,hspace=.35)

ax[0,0].set_xlim(-2000,2000)
ax[0,0].set_ylim(-2000,2000)
ax[0,0].text(-2500,2300,"a)",fontsize=16)
ax[0,0].set_aspect("equal")
ax[0,0].set_ylabel("Y [km]",fontsize=16,labelpad=-10)
ax[0,0].set_xlabel("X [km]",fontsize=16)
ax[0,0].tick_params(axis="both",labelsize=14)
ccm0 = cmocean.cm.curl
ccm0.set_bad(color='gray')
p0=ax[0,0].pcolormesh(xpc/1000,ypc/1000,np.ma.array(mo50.eta[0,:,:],mask=(g1.mask_h==0)),vmin=hcblim[0],vmax=hcblim[-1],cmap=ccm0)

ax[0,0].set_title("$\\epsilon_L = $ %.2f " % (g1.epsilon_L),fontsize=20)
q0 = ax[0,0].quiver(xmath[ss]/1000,ymath[ss]/1000,
               np.ma.array(mo50.uu[0,ss,ss],mask=(g1.mask_u[ss,ss]==0)),
               np.ma.array(mo50.vv[0,ss,ss],mask=(g1.mask_v[ss,ss]==0)),scale=scale,color=color,width=width)
qk0 = ax[0,0].quiverkey(q0,xk,yk,.1, '$0.1 \ \\mathrm{m \ s^{-1}}$',coordinates = 'axes',fontproperties={'size':16})


ax[0,1].text(-2500,2300,"b)",fontsize=16)
ax[0,1].set_xlim(-2000,2000)
ax[0,1].set_ylim(-2000,2000)
ax[0,1].set_aspect("equal")
ax[0,1].set_ylabel("Y [km]",fontsize=16,labelpad=-10)
ax[0,1].set_xlabel("X [km]",fontsize=16)
ax[0,1].tick_params(axis="both",labelsize=14)
ccm0 = cmocean.cm.curl
ccm0.set_bad(color='gray')
p0=ax[0,1].pcolormesh(xpc/1000,ypc/1000,np.ma.array(mo150.eta[0,:,:],mask=(g2.mask_h==0)),vmin=hcblim[0],vmax=hcblim[-1],cmap=ccm0)

ax[0,1].set_title("$\\epsilon_L = $ %.2f " % (g2.epsilon_L),fontsize=20)
q0 = ax[0,1].quiver(xmath[ss]/1000,ymath[ss]/1000,
               np.ma.array(mo150.uu[0,ss,ss],mask=(g2.mask_u[ss,ss]==0)),
               np.ma.array(mo150.vv[0,ss,ss],mask=(g2.mask_v[ss,ss]==0)),scale=scale,color=color,width=width)
qk0 = ax[0,1].quiverkey(q0,xk,yk,.1, '$0.1 \ \\mathrm{m \ s^{-1}}$',coordinates = 'axes',fontproperties={'size':16})

ax[0,2].text(-2500,2300,"c)",fontsize=16)
ax[0,2].set_xlim(-2000,2000)
ax[0,2].set_ylim(-2000,2000)
ax[0,2].set_aspect("equal")
ax[0,2].set_ylabel("Y [km]",fontsize=16,labelpad=-10)
ax[0,2].set_xlabel("X [km]",fontsize=16)
ax[0,2].tick_params(axis="both",labelsize=14)
ccm0 = cmocean.cm.curl
ccm0.set_bad(color='gray')
p0=ax[0,2].pcolormesh(xpc/1000,ypc/1000,np.ma.array(mo200.eta[0,:,:],mask=(g3.mask_h==0)),vmin=hcblim[0],vmax=hcblim[-1],cmap=ccm0)
q0 = ax[0,2].quiver(xmath[ss]/1000,ymath[ss]/1000,
               np.ma.array(mo200.uu[0,ss,ss],mask=(g3.mask_u[ss,ss]==0)),
               np.ma.array(mo200.vv[0,ss,ss],mask=(g3.mask_v[ss,ss]==0)),scale=scale,color=color,width=width)
qk0 = ax[0,2].quiverkey(q0,xk,yk,.1, '$0.1 \ \\mathrm{m \ s^{-1}}$',coordinates = 'axes',fontproperties={'size':16})
ax[0,2].set_title("$\\epsilon_L = $ %.2f " % (g3.epsilon_L),fontsize=20)
cax= inset_axes(ax[0,2],
                   width="5%",  # width = 5% of parent_bbox width
                   height="100%",  # height : 50%
                   loc='lower right',
                   bbox_to_anchor=(.13, 0, 1, 1),
                   bbox_transform=ax[0,2].transAxes,
                   borderpad=0,
                   )
cb=plt.colorbar(p0,cax=cax,orientation="vertical",extend="both")
cb.set_label('$\eta$ [m]', fontsize=16)
cb.ax.tick_params(labelsize=14)
cb.set_ticks(ticks)
fig.suptitle("%.2f Inertial Periods" % (mo50.tdays*g1.fcor[1]),fontsize=20,y=.94)


################################ Spectra #######################################
ff = "../remote_wind_forcing/subset_data/"
p1 = ff + 'r50/'
p4 = ff + 'r150/'
p5 = ff + 'r200/'



def spec_r(TS,dt,nsmooth):
    spec = spectra.spectrum(TS,dt=dt, nfft=None,window='quadratic')
    nmaxi = np.argmax(spec.psd)
    fmax = spec.freqs[nmaxi]  ### peak frequency

    return spec,fmax



def load_var(bpath,n_filename):

    aN = np.load(bpath + n_filename)



    b = Bunch()
    b.u = aN["u"]
    b.v = aN["v"]
    b.h = aN["h"]
    b.eta = aN["eta"]
    b.tt = aN["t"]

    return b
fname="E_ts.npz"

r50 = load_var(p1,fname)
r150 = load_var(p4,fname)
r200 = load_var(p5,fname)


dtdays = r50.tt[1] -  r50.tt[0]
nsmooth=1


TS1 = r50.eta[:431]
spec,fmax = spec_r(TS1,dtdays,nsmooth)
TS2 = r150.eta[:431]
spec2,fmax2 = spec_r(TS2,dtdays,nsmooth)
TS3 = r200.eta[:431]
spec3,fmax3 = spec_r(TS3,dtdays,nsmooth)

def spec_r(TS,dt,nsmooth):
    spec = spectra.spectrum(TS,dt=dt, nfft=None,window='quadratic')
    nmaxi = np.argmax(spec.cwpsd)
    fmax = spec.freqs[nmaxi]  ### peak frequency

    return spec,fmax


TS1 = r50.u[:431] + 1j*r50.v[:431]
specv,fmax = spec_r(TS1,dtdays,nsmooth)
TS2 = r150.u[:431] + 1j*r150.v[:431]
spec2v,fmax2 = spec_r(TS2,dtdays,nsmooth)
TS3 = r200.u[:431]+ 1j*r200.v[:431]
spec3v,fmax3 = spec_r(TS3,dtdays,nsmooth)

ax[1,0].text(3e-2,1e-3,"d)",fontsize=16)
ax[1,0].set_xlabel("cpd",fontsize=14)
ax[1,0].set_ylabel("$\\mathrm{m^2 \ cpd^{-1}}$ / $\\mathrm{m^2 \ s^{-2} cpd^{-1}}$",fontsize=14)
ax[1,0].loglog(spec.freqs,spec.psd,'r',label="$\\eta$")
ax[1,0].loglog(specv.cwfreqs,specv.cwpsd,'b',label="CW")
ax[1,0].loglog(specv.ccwfreqs,specv.ccwpsd,'b--',label="CCW")

ax[1,0].tick_params(axis="both",labelsize=12)
ax[1,0].set_title("$\\epsilon_L = $ %.2f " % (g1.epsilon_L),fontsize=20)
ax[1,0].plot([g1.fcor[1],g1.fcor[1]],[1e-12,1e-3],'--k')
ax[1,0].set_xlim([spec.freqs[0],spec.freqs[-1]])
ax[1,0].set_ylim([1e-12,1e-3])
ax[1,0].text(g1.fcor[1]+.004,4e-12,"$f$",fontsize=18)
ax[1,0].text(1.5e1,8e-3,"Power Spectral Density: $\\eta$ and $\\vec{u}$",fontsize=20)
ax[1,0].legend(fontsize=16)

ax[1,1].text(3e-2,1e-3,"e)",fontsize=16)
ax[1,1].set_xlabel("cpd",fontsize=14)
ax[1,1].set_ylabel("$\\mathrm{m^2 \ cpd^{-1}}$ / $\\mathrm{m^2 \ s^{-2} cpd^{-1}}$",fontsize=14)
ax[1,1].tick_params(axis="both",labelsize=12)
ax[1,1].set_title("$\\epsilon_L = $ %.2f " % (g2.epsilon_L),fontsize=20)
ax[1,1].loglog(spec2.freqs,spec2.psd,'r',label="$\\eta$")
ax[1,1].loglog(spec2v.cwfreqs,spec2v.cwpsd,'b',label="CW")
ax[1,1].loglog(spec2v.ccwfreqs,spec2v.ccwpsd,'b--',label="CCW")
ax[1,1].legend(fontsize=16)

ax[1,1].plot([g1.fcor[1],g1.fcor[1]],[1e-12,1e-3],'--k')
ax[1,1].set_xlim([spec.freqs[0],spec.freqs[-1]])
ax[1,1].set_ylim([1e-12,1e-3])
ax[1,1].text(g1.fcor[1]+.004,4e-12,"$f$",fontsize=18)
ax[1,2].text(3e-2,1e-3,"f)",fontsize=16)

ax[1,2].loglog(spec3.freqs,spec3.psd,'r',label="$\\eta$")
ax[1,2].loglog(spec3v.cwfreqs,spec3v.cwpsd,'b',label="CW")
ax[1,2].loglog(spec3v.ccwfreqs,spec3v.ccwpsd,'b--',label="CCW")
ax[1,2].set_xlabel("cpd",fontsize=14)
ax[1,2].set_ylabel("$\\mathrm{m^2 \ cpd^{-1}}$ / $\\mathrm{m^2 \ s^{-2} cpd^{-1}}$",fontsize=14)
ax[1,2].tick_params(axis="both",labelsize=12)
ax[1,2].set_title("$\\epsilon_L = $ %.2f " % (g3.epsilon_L),fontsize=20)
ax[1,2].plot([g1.fcor[1],g1.fcor[1]],[1e-12,1e-3],'--k')
ax[1,2].set_xlim([spec.freqs[0],spec.freqs[-1]])
ax[1,2].set_ylim([1e-12,1e-3])
ax[1,2].text(g1.fcor[1]+.004,4e-12,"$f$",fontsize=18)
ax[1,2].legend(fontsize=16)

fig.savefig(figpath+ svname, dpi=300, bbox_inches='tight',orientation="horizontal",pad_inches=0.3)
