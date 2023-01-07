
import sys
sys.path.append('../src')
from ms_toolbox import tools as tls

import matplotlib.pyplot as plt
import numpy as np

from ms_toolbox import gfd
from ms_toolbox import rgmodel as rgm
from pycurrents.num import spectra, rangeslice

from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy.interpolate import interp2d as int2d

from netCDF4 import Dataset
from pycurrents.system import Bunch
from matplotlib.ticker import FormatStrFormatter
import matplotlib
import cmocean
import matplotlib.cm as cm
import scipy.stats as ss

"""
Circle with 100 km diameter

"""
################# Define Paths ###################

fpath = '../model_output/small_islands/r50/' #path to files from UCSD Library Collection

ncname = 'C_LR_FS_r50_Ah00_r0_gp09'

ncFname = fpath + ncname + '.nc'
ncFnameG = fpath + ncname + '_grid.nc'

figpath = '../small_islands/r50/figures/' #path from git directory

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


######## Find closest to point to make TS ##########
mindist_n = rgm.min_dist(g.xx,g.yy,3000e3,3050e3)
mindist_s = rgm.min_dist(g.xx,g.yy,3000e3,2950e3)
mindist_e = rgm.min_dist(g.xx,g.yy,3050e3,3000e3)
mindist_w = rgm.min_dist(g.xx,g.yy,2950e3,3000e3)
mindist_ff = rgm.min_dist(g.xx,g.yy,7700e3,5700e3)
mindist_ff2 = rgm.min_dist(g.xx,g.yy,4000e3,5500e3)

ele_n = [mindist_n[0],mindist_n[1]] #y,x
ele_e = [mindist_e[0],mindist_e[1]] #y,x
ele_s = [mindist_s[0],mindist_s[1]] #y,x
ele_w = [mindist_w[0],mindist_w[1]] #y,x
ele_ff = [mindist_ff[0],mindist_ff[1]] #y,x
ele_ff2 = [mindist_ff2[0],mindist_ff2[1]] #y,x

##################################################
########################## functions #############################
def calculate_divergence(mo,dx,dy):
    dudx = ((mo.uu[:,:,1:] - mo.uu[:,:,:-1]) / dx)
    dvdy = ((mo.vv[:,1:,:] - mo.vv[:,:-1,:]) / dy)

    dudx_m = .5 * (dudx[:,1:,:] + dudx[:,:-1,:])
    dvdy_m = .5 * (dvdy[:,:,1:] + dvdy[:,:,:-1])

    div = dudx_m + dvdy_m
    return div

def dispr(gp,H,k,f):

    omega = np.sqrt((f**2) + (gp*H*(k**2)))
    lamb = np.sqrt(gp * H) / f
    kl = (lamb**2) * (k**2)
    omegf = np.sqrt(1 + kl)

    return omega,omegf,lamb,np.sqrt(kl)

################# Calculate Divergence ##################
dx = np.mean(np.diff(g.x))
dy = np.mean(np.diff(g.y))

div = calculate_divergence(mo,dx,dy)

divtsN = div[:,ele_n[0],ele_n[1]]
divtsE = div[:,ele_e[0],ele_e[1]+1]
etaE = mo.eta[:,ele_e[0],ele_e[1]]

htsN = mo.hh[:,ele_n[0],ele_n[1]]
htsE = mo.hh[:,ele_e[0],ele_e[1]]

Ny = len(g.y)
Nx = len(g.x)
indy = int(np.floor(Ny/2))
indx = int(np.floor(Nx/2))

xmath = g.x - g.x[indx]
ymath = g.y - g.y[indy]


inmid_y = int(np.floor(Ny / 2))
inmid_x = int(np.floor(Nx / 2))

hhtran = mo.eta[:432,inmid_y,:]

t = mo.tdays


######################  Estimate Wavelength and frequency  ###############################
nn = np.shape(mo.hh)
xh = int(np.floor(nn[2]/2))
hr = np.squeeze(mo["hh"][:,indy,xh:]).T
xr = xmath[xh:]


inminh = hr.argmin(axis=0)
inmaxh = hr.argmax(axis=0)

fig,ax = plt.subplots()
ax.set_title(("%d Y index" % indy))
pc=ax.pcolormesh(mo.tinertial,xr / 1000,hr,vmin=-.001,vmax=.001,cmap=cmocean.cm.diff)
ax.plot(mo.tinertial,xr[inminh]/1000,'.r')
ax.plot(mo.tinertial,xr[inmaxh]/1000,'.r')
ax.set_xlabel("Time (days)")
fig.colorbar(pc)


trng1 = rangeslice(mo.tinertial,.21,.57)
trng2 = rangeslice(mo.tinertial,.66,1.02)
trng3 = rangeslice(mo.tinertial,1.1,1.55)

tshort_sec1 = mo.tt[trng1]
xshortMAX1 = xr[inmaxh][trng1]
tshort_ip1 = mo.tinertial[trng1]

tshort_sec2 = mo.tt[trng2]
xshortMIN2 = xr[inminh][trng2]
tshort_ip2 = mo.tinertial[trng2]


tshort_sec3 = mo.tt[trng3]
xshortMAX3 = xr[inmaxh][trng3]
tshort_ip3 = mo.tinertial[trng3]

slope_sec1, intercept_sec1, r, p, stderr = ss.linregress(tshort_sec1, xshortMAX1)
slope_sec2, intercept_sec2, r, p, stderr = ss.linregress(tshort_sec2, xshortMIN2)
slope_sec3, intercept_sec3, r, p, stderr = ss.linregress(tshort_sec3, xshortMAX3)

slope_ip, intercept_ip, r, p, stderr = ss.linregress(tshort_ip2, xshortMIN2)
yfit_ip = slope_ip*np.array([.1,3]) + intercept_ip

yfit1 = slope_sec1*tshort_sec1 + intercept_sec1
yfit2 = slope_sec2*tshort_sec2 + intercept_sec2
yfit3 = slope_sec3*tshort_sec3 + intercept_sec3

yfit_ip1 = slope_sec1*np.array([0,10*86400])+ intercept_sec1
yfit_ip2 = slope_sec2*np.array([0,10*86400])+ intercept_sec2
yfit_ip3 = slope_sec3*np.array([0,10*86400])+ intercept_sec3

#
# svname="Estimate_wavelength_frequency.png"
# fig,ax = plt.subplots(figsize=(10,5))
# ax.set_title(("%d Y index" % indy))
# pc=ax.pcolormesh(mo.tt,xr / 1000,hr,vmin=-.001,vmax=.001,cmap=cmocean.cm.diff)
# ax.plot([0,10*86400],yfit_ip1/1000,'r')
# ax.plot([0,10*86400],yfit_ip2/1000,'r')
# ax.plot([0,10*86400],yfit_ip3/1000,'r')
# ax.plot([300e3,300e3],[0,1500],'r')
# ax.plot([0,400e3],[560,560],'r')
# ax.set_ylim([0,1500])
# ax.set_xlim([0,400000])
# ax.set_xlabel("Time (days)")
# fig.colorbar(pc)

# fig.savefig(figpath+ svname, dpi=300, bbox_inches='tight',orientation="horizontal",pad_inches=0.2)

k_est = (2*np.pi) / (844e3-209e3)
om_est = (2*np.pi) / (390243-211809)

cp = om_est/k_est


k = np.arange(0,.5,.000001)
gp = .09
H = 100
f = g.fcor[0]
omega09,omegaf09,lam09,kl09 = dispr(gp,H,k,f)

# fig,ax=plt.subplots()
# ax.plot(np.log(k),np.log(omega09),label="g' = .09")
# ax.plot([-14,-1],[np.log(g.fcor[0]),np.log(g.fcor[0])],'--k')
# ax.plot(np.log(k_est),np.log(om_est),'*r',markersize=10)
# ax.set_xlim([-14,-1])
# ax.set_ylabel("log($\\omega$)",fontsize=20)
# ax.set_xlabel("log($k$)",fontsize=20)

#################### MAKE FIGURE 6 ##########################
xd = 10
yd = 8

cx = .15
cy = .2
ra = yd / xd
px1 = .8

csfont = {'fontname':'Times New Roman'}

matplotlib.rcParams['mathtext.fontset'] = 'stix'

svname = "FIGURE_06.png"
fig = plt.figure(figsize=(xd,yd))
ax0 = fig.add_axes([cx,.82,px1,.15])
ax1 = fig.add_axes([cx,.6,px1,.15])
ax2 = fig.add_axes([cx,.1,.35,.4])
ax3 = fig.add_axes([cx+.45,.1,.35,.4])

ax0.plot(mo.tinertial[:432],mo.taux[:432],label="$\\tau_x$")
ax0.plot(mo.tinertial[:432],mo.tauy[:432],label="$\\tau_y$")
ax0.legend(loc="upper right",fontsize=14)
ax0.set_xlim([mo.tinertial[0],1])
ax0.set_ylim([-.05,.05])
ax0.tick_params(axis="both",labelsize=12)
ax0.set_xlabel("Inertial Periods",fontsize=12)
ax0.text(.01, .03,"Wind Stress",fontsize=16)
ax0.set_ylabel("$\\mathrm{kg m^{-1}s^{-1}}$",fontsize=12)
ax0.text(-.15,.06,"a)",fontsize=16)
ax0.grid("on")

ax1.plot(mo.tinertial[:432],divtsE[:432],'r')
ax1.set_xlim([mo.tinertial[0],1])
ax1.tick_params(axis="both",labelsize=12)
ax1.set_xlabel("Inertial Periods",fontsize=12)
ax1.text(.01,3e-7,"$\\nabla \\cdot \\vec{u}$",fontsize=16)
ax1.set_ylabel("$\\mathrm{s^{-1}}$",fontsize=12)
ax1.text(-.15,7e-7,"b)",fontsize=16)
ax1.set_ylim([-5e-7,5e-7])
ax1.grid("on")


fig.tight_layout(pad=.8)
pc0=ax2.pcolormesh(mo.tinertial[:432],(xmath/1000-(g.Xrad/1000)),hhtran.T,cmap=cmocean.cm.diff,vmin=-.001,vmax=.001)
ax2.plot([.1,3],(yfit_ip/1000)-(g.Xrad/1000),'--',color="c")
ax2.set_ylabel("Distance from Island Boundary [$\\mathrm{km}$]",fontsize=12)
ax2.set_xlabel("Inertial Periods",fontsize=12)
ax2.set_title("$\\eta$",fontsize=16)
ax2.set_xlim([mo.tinertial[0],3])
ax2.set_ylim([0,1000])
ax2.tick_params(axis="both",labelsize=12)

cax= inset_axes(ax2,
                   width="100%",  # width = 5% of parent_bbox width
                   height="5%",  # height : 50%
                   loc='lower left',
                   bbox_to_anchor=(0, -.25, 1, 1),
                   bbox_transform=ax2.transAxes,
                   borderpad=0,
                   )
cb1=fig.colorbar(pc0,cax=cax,orientation="horizontal",extend="both")
cb1.set_ticks([-.001,0,.001])
cb1.ax.tick_params(axis="both",labelsize=13)
cb1.set_label('$\\eta [m]$',fontsize=13)
ax2.text(-1.01,1100,"c)",fontsize=16)

ax3.loglog(k,omega09,label="g' = .09")
ax3.plot([5e-7,1e-3],[g.fcor[0],g.fcor[0]],'--k')
ax3.plot(k_est,om_est,'*r',markersize=10)
ax3.set_xlim([5e-7,1e-3])
ax3.set_ylim([5e-6,1e-2])
ax3.set_ylabel("$\\omega$ [$\\mathrm{radians \ s^{-1}}$]",fontsize=14)
ax3.set_xlabel("$k$ [$\\mathrm{radians \ m^{-1}}$]",fontsize=14)
ax3.text(1e-7,1.4e-2,"d)",fontsize=16)
ax3.text(6e-7,1.2e-5,"$f$",fontsize=16)
ax3.set_title("Poincaré Wave Dispersion Relationship",fontsize=12)

fig.savefig(figpath+ svname, dpi=300, bbox_inches='tight',orientation="horizontal",pad_inches=0.2)
