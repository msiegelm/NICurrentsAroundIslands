
import sys
sys.path.append('../src')

import matplotlib.pyplot as plt
import numpy as np

from ms_toolbox import gfd
from ms_toolbox import rgmodel as rgm
from ms_toolbox import tools as tls
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
Makes Momentum Balance Plot
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
g.r = 0

#################################################

######## Check to confirm output slice ##########
fig,ax = plt.subplots()
pc = ax.pcolormesh(g.el)
ax.set_aspect("equal")
#################################################

################# Read Output ###################

mo = rgm.read_output(xrng,yrng,ncFname)
#################################################

######## Check to confirm output for TS ##########

fig,ax = plt.subplots()
pc = ax.pcolormesh(g.x/1000,g.y/1000,g.el)
ax.set_aspect("equal")
ax.axes.set_xlabel("km")
ax.axes.set_ylabel("km")
#################################################


######## Find closest to point to make TS ##########
mindist_n = rgm.min_dist(g.xx,g.yy,3000e3,3050e3)
mindist_s = rgm.min_dist(g.xx,g.yy,3000e3,2950e3)
mindist_e = rgm.min_dist(g.xx,g.yy,3050e3,3000e3)
mindist_w = rgm.min_dist(g.xx,g.yy,2950e3,3000e3)
mindist_ff = rgm.min_dist(g.xx,g.yy,5500e3,5500e3)
mindist_ff2 = rgm.min_dist(g.xx,g.yy,4000e3,5500e3)

ele_n = [mindist_n[0],mindist_n[1]] #y,x
ele_e = [mindist_e[0],mindist_e[1]] #y,x
ele_s = [mindist_s[0],mindist_s[1]] #y,x
ele_w = [mindist_w[0],mindist_w[1]] #y,x
ele_ff = [mindist_ff[0],mindist_ff[1]] #y,x
ele_ff2 = [mindist_ff2[0],mindist_ff2[1]] #y,x

###################################################

############# Plot of selected points ############
#
# fig,ax = plt.subplots(figsize=(8,8))
# pc = ax.pcolormesh(g.x/1000,g.y/1000,g.el)
# ax.plot(g.x[ele_n[1]]/1000,g.y[ele_n[0]]/1000,"*",label="N",markersize=15)
# ax.plot(g.x[ele_e[1]]/1000,g.y[ele_e[0]]/1000,"*",label="E",markersize=15)
# ax.plot(g.x[ele_s[1]]/1000,g.y[ele_s[0]]/1000,"*",label="S",markersize=15)
# ax.plot(g.x[ele_w[1]]/1000,g.y[ele_w[0]]/1000,"*",label="W",markersize=15)
# ax.plot(g.x[ele_ff[1]]/1000,g.y[ele_ff[0]]/1000,"*",label="FF",markersize=15)
# ax.set_aspect("equal")
# ax.tick_params(labelsize=16)
# ax.set_ylabel("km",fontsize=16)
# ax.set_xlabel("km",fontsize=16)
# ax.legend(fontsize=16,loc=(.0,1),frameon=False,ncol=5)
# # ax.plot([-10,600],[-10,600])
# fig.savefig(figpath + 'grid' + ".png", dpi=300, bbox_inches='tight',orientation="horizontal")
# ################### PLOTTING FUNCTION ###############################
def plot_momentum_terms(g,xmath,ymath,mtu,mtv,mo,nt,nskip,nskip2,xyrng,xyrng2,ticks,ticks2,vrng,hrng,sca,save=False,svpath=[],svname=[]):
    # ax00.set_aspect("equal")
    sl = np.s_[::nskip]
    sl2 = np.s_[::nskip2]
    xd = 8
    yd = 10 #7
    ra = yd / xd
    px = .24
    cx = .05
    cy = .08

    px2 = .4
    r1h = cy + .28
    xs = .32 # shifts fig to the right

    pp = np.s_[1:]
    mm = np.s_[:-1]

    xlim = (xyrng[0],xyrng[1])
    ylim = (xyrng[0],xyrng[1])

    xlim2 = (xyrng2[0],xyrng2[1])
    ylim2 = (xyrng2[0],xyrng2[1])

    xdx = .5*(g.x[pp] + g.x[mm])
    ydy = .5*(g.y[pp] + g.y[mm])
    xdx2 = .5*(xdx[pp] + xdx[mm])
    ydy2 = .5*(ydy[pp] + ydy[mm])

    fig = plt.figure(figsize=(xd,yd))
    ax10 = fig.add_axes([cx,cy,px,px / ra])
    ax00 = fig.add_axes([cx,r1h,px,px / ra])

    ax11 = fig.add_axes([cx+(xs*1),cy,px,px / ra])
    ax01 = fig.add_axes([cx+(xs*1),r1h,px,px / ra])

    ax12 = fig.add_axes([cx+(2*xs),cy,px,px / ra])
    ax02 = fig.add_axes([cx+(2*xs),r1h,px,px / ra])
    ax11cb = fig.add_axes([cx+xs-.06,cy-.08,px+.12,(px / ra)*.1])

    # ax13 = fig.add_axes([cx+(3*xs),cy,px,px / ra])
    # ax03 = fig.add_axes([cx+(3*xs),r1h,px,px / ra])

    # ax14 = fig.add_axes([cx+(4*xs),cy,px,px / ra])

    axssh = fig.add_axes([cx,r1h + .28,px2,px2 / ra])
    axssh2 = fig.add_axes([cx + .48,r1h + .28,px2,px2 / ra])
    axsshcb = fig.add_axes([cx+xs+ .61,r1h + .28,px2*.08,(px2 / ra)])
    # ax04 = fig.add_axes([cx+(4*xs),r1h,px,px / ra])

    vmin_= vrng[0]
    vmax_ = vrng[1]

    vminh_ = hrng[0]
    vmaxh_ = hrng[1]
    ccm = cm.seismic
    ccm.set_bad(color='dimgray')

    ccm2 = cmocean.cm.curl
    ccm2.set_bad(color='dimgray')

    axssh.pcolormesh(xmath,ymath,np.ma.array(mo.eta[nt,:,:],mask=(g.mask_h==0)),cmap=ccm2,vmin=vminh_,vmax=vmaxh_)
    axssh.set_aspect("equal")
    axssh.set_xlim([xlim[0],xlim[1]])
    axssh.set_ylim([ylim[0],ylim[1]])
    axssh.set_xticks(ticks)
    axssh.set_yticks(ticks)
    axssh.set_ylabel("$ya^{-1}$",fontsize=14)
    axssh.set_xlabel("$xa^{-1}$",fontsize=14)
    axssh.tick_params(axis="both",labelsize=14)
    qq = axssh.quiver(xmath[sl],ymath[sl],np.ma.array(mo.uu[nt,sl,sl],mask=(g.mask_u[sl,sl]==0)),np.ma.array(mo.vv[nt,sl,sl],mask=(g.mask_v[sl,sl]==0)),scale=sca,width=.005,color='dodgerblue')
    axssh.quiverkey(qq,1.05,1.05,.2,"$0.2 \ \\mathrm{ms^{-1}}$",fontproperties={'size': 16})
    # axssh.set_title("%.2f Inertial Periods" % (mo.tdays[nt]*g.fcor[1]),fontsize=20)

    axssh2.pcolormesh(xmath,ymath,np.ma.array(mo.eta[nt,:,:],mask=(g.mask_h==0)),cmap=ccm2,vmin=vminh_,vmax=vmaxh_)
    axssh2.set_aspect("equal")
    axssh2.set_xlim([xlim2[0],xlim2[1]])
    axssh2.set_ylim([ylim2[0],ylim2[1]])
    axssh2.set_xticks(ticks2)
    axssh2.set_yticks(ticks2)
    axssh2.set_ylabel("$ya^{-1}$",fontsize=14)
    axssh2.set_xlabel("$xa^{-1}$",fontsize=14)
    axssh2.yaxis.labelpad = -.01
    axssh2.tick_params(axis="both",labelsize=14)
    qq = axssh2.quiver(xmath[sl2],ymath[sl2],np.ma.array(mo.uu[nt,sl2,sl2],mask=(g.mask_u[sl2,sl2]==0)),np.ma.array(mo.vv[nt,sl2,sl2],mask=(g.mask_v[sl2,sl2]==0)),scale=sca,width=.005,color='dodgerblue')
    axssh2.quiverkey(qq,1.05,1.05,.2,"$0.2 \ \\mathrm{ms^{-1}}$",fontproperties={'size': 16})
    if nt == 86:
        fig.suptitle("One Inertial Period" % (mo.tdays[nt]*g.fcor[1]),fontsize=20,y=1.05)
    else:
        fig.suptitle("%.2f Inertial Periods" % (mo.tdays[nt]*g.fcor[1]),fontsize=20,y=1.05)


    sm7 = plt.cm.ScalarMappable(cmap=cmocean.cm.curl)
    sm7.set_clim(vminh_,vmaxh_)
    sm7.set_array(mtu.viscx)
    cb7 = fig.colorbar(sm7,cax=axsshcb, orientation="vertical")
    # cb7.formatter.set_powerlimits((0,0))
    cb7.ax.tick_params(axis="both",labelsize=18)
    cb7.ax.xaxis.offsetText.set_fontsize(16)
    cb7.set_label('$\\eta [\\mathrm{m}]$',fontsize=16)



    ax00.pcolormesh(xmath,ymath,np.ma.array(mtu.dudt,mask=(mtu.dudt==0)),cmap=ccm,vmin=vmin_,vmax=vmax_)
    ax00.set_aspect("equal")
    ax00.set_xlim([xlim[0],xlim[1]])
    ax00.set_ylim([ylim[0],ylim[1]])
    ax00.text(3.5,3.5,'$\\frac{\\partial u}{\\partial t}$',fontsize=24)
    ax00.set_xticks(ticks)
    ax00.set_yticks(ticks)
    ax00.set_ylabel("$ya^{-1}$",fontsize=14,labelpad=-10)
    ax00.set_xlabel("$xa^{-1}$",fontsize=14,labelpad=-5)
    ax00.tick_params(axis="both",labelsize=14)


    ax01.pcolormesh(xmath,ymath,np.ma.array(mtu.coriolis,mask=(mtu.coriolis==0)),vmin=vmin_,vmax=vmax_,cmap=cm.seismic)
    ax01.set_aspect("equal")
    ax01.set_xlim([xlim[0],xlim[1]])
    ax01.set_ylim([ylim[0],ylim[1]])
    ax01.text(2,3.5,'$-fv$',fontsize=24)
    ax01.set_xticks(ticks)
    ax01.set_yticks(ticks)
    ax01.tick_params(axis="both",labelsize=14)
    ax01.set_title("x-momentum balance",fontsize=20)
    ax01.set_ylabel("$ya^{-1}$",fontsize=14,labelpad=-10)
    ax01.set_xlabel("$xa^{-1}$",fontsize=14,labelpad=-5)
    ax01.tick_params(axis="both",labelsize=14)


    ax02.pcolormesh(xmath,ymath,np.ma.array(mtu.pg,mask=(mtu.pg==0)),cmap=ccm,vmin=vmin_,vmax=vmax_)
    ax02.set_aspect("equal")
    ax02.set_xlim([xlim[0],xlim[1]])
    ax02.set_ylim([ylim[0],ylim[1]])
    ax02.text(1.5,3.5,'$ -g\' \\frac{\\partial h}{\\partial x}$',fontsize=24)
    ax02.set_xticks(ticks)
    ax02.set_yticks(ticks)
    ax02.set_ylabel("$ya^{-1}$",fontsize=14,labelpad=-10)
    ax02.set_xlabel("$xa^{-1}$",fontsize=14,labelpad=-5)
    ax02.tick_params(axis="both",labelsize=14,pad=.08)

    ax10.pcolormesh(xmath,ymath,np.ma.array(mtv.dvdt,mask=(mtv.dvdt==0)),cmap=ccm,vmin=vmin_,vmax=vmax_)
    ax10.set_aspect("equal")
    ax10.set_xlim([xlim[0],xlim[1]])
    ax10.set_ylim([ylim[0],ylim[1]])
    ax10.text(3.5,3.5,'$\\frac{\\partial v}{\\partial t}$',fontsize=24)
    ax10.set_xticks(ticks)
    ax10.set_yticks(ticks)
    ax10.set_ylabel("$ya^{-1}$",fontsize=14,labelpad=-10)
    ax10.set_xlabel("$xa^{-1}$",fontsize=14,labelpad=-5)
    ax10.tick_params(axis="both",labelsize=14,pad=.08)


    ax11.pcolormesh(xmath,ymath,np.ma.array(mtv.coriolis,mask=(mtv.coriolis==0)),vmin=vmin_,vmax=vmax_,cmap=cm.seismic)
    ax11.set_aspect("equal")
    ax11.set_xlim([xlim[0],xlim[1]])
    ax11.set_ylim([ylim[0],ylim[1]])
    ax11.text(3.2,3.5,'$fu$',fontsize=24)
    ax11.set_xticks(ticks)
    ax11.set_yticks(ticks)
    ax11.set_title("y-momentum balance",fontsize=20)
    ax11.set_xlabel("$xa^{-1}$",fontsize=14,labelpad=-5)
    ax11.set_ylabel("$ya^{-1}$",fontsize=14,labelpad=-10)
    ax11.tick_params(axis="both",labelsize=14)

    ax12.pcolormesh(xmath,ymath,np.ma.array(mtv.pg,mask=(mtv.pg==0)),cmap=ccm,vmin=vmin_,vmax=vmax_)
    ax12.set_aspect("equal")
    ax12.set_xlim([xlim[0],xlim[1]])
    ax12.set_ylim([ylim[0],ylim[1]])
    ax12.text(1.5,3.5,'$ -g\' \\frac{\\partial h}{\\partial y}$',fontsize=24)
    ax12.set_xticks(ticks)
    ax12.set_yticks(ticks)
    # ax01.set_xlabel("$xa^{-1}$",fontsize=14)
    ax12.set_ylabel("$ya^{-1}$",fontsize=14,labelpad=-10)
    ax12.set_xlabel("$xa^{-1}$",fontsize=14,labelpad=-5)
    ax12.tick_params(axis="both",labelsize=14)


    sm6 = plt.cm.ScalarMappable(cmap=cm.seismic)
    sm6.set_clim(vmin_,vmax_)
    sm6.set_array(mtu.pg)
    cb6 = fig.colorbar(sm6,cax=ax11cb, orientation="horizontal")
    cb6.formatter.set_powerlimits((0,0))
    cb6.ax.tick_params(axis="both",labelsize=18)
    cb6.ax.xaxis.offsetText.set_fontsize(16)
    cb6.set_label('$[\\mathrm{ms^{-2}}]$',fontsize=16)

    if save:
        fig.savefig(svpath+ svname, dpi=300, bbox_inches='tight',orientation="horizontal",pad_inches=0.2)

##############################################################
nt = 86 ## 305 works
mtu,mtv = at.momentum_terms(nt,g,mo,g.gp,g.r)

nskip=15
nskip2=5
vrng = (-5e-6,5e-6)
hrng = (-.01,.01)
xyrng = (-10,10)
xyrng2 = (-5,5)

ticks = np.arange(-10,xyrng[-1]+5,5)
ticks2 = (-4,-2,0,2,4)

xmid = int(np.floor(len(g.x)/2))
ymid = int(np.floor(len(g.y)/2))

xmath = (g.x-g.x[xmid])/g.Xrad
ymath = (g.y-g.y[ymid])/g.Xrad

sca = 1.5

svnameu = ("FIGURE_07.png")

plot_momentum_terms(g,xmath,ymath,mtu,mtv,mo,nt,nskip,nskip2,xyrng,xyrng2,ticks,ticks2,vrng,hrng,sca,save=True,svpath=figpath,svname=svnameu)
