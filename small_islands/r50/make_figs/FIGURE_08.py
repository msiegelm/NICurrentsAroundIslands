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
import matplotlib as mpl
import matplotlib
matplotlib.rcParams['mathtext.fontset'] = 'stix'


anpath = '../small_islands/r50/data_subsection/' #path from git directory
figpath = '../small_islands/r50/figures/' #path from git directory

####################################################


g09 = np.load(anpath + "normalization_snapshot_FINAL.npz")


rstar = np.arange(1,10+.05,.05)
ur_scale = (1 - (1/(rstar**2)))
utheta_scale = (1 + (1/(rstar**2)))
zeta_scale = (-2/rstar)
#%%

ts = g09["td"][:]
uu = g09["uu"][:]
nt = g09["nt"][:]



cmap = cmocean.cm.amp
norm = mpl.colors.Normalize(vmin=0,vmax=max(ts))
sm = mpl.cm.ScalarMappable(norm=norm,cmap=cmap)
sm.set_array([])

svname = "FIGURE_08.png"
fig,ax=plt.subplots(ncols=1,nrows=4,figsize=(14,10))
ax[0].plot(ts,uu)
sc=ax[0].scatter(ts[nt],uu[nt],c=ts[nt],cmap=cmocean.cm.amp,vmin=0,vmax=max(ts))
cb=fig.colorbar(sc,ax=ax)
cb.set_label("Inertial Periods",fontsize=14)
cb.ax.tick_params(labelsize=16)
ax[0].plot([ts[86],ts[86]],[-.1,.1],'--k')
ax[0].set_title("$u_\\theta \ (r=\\infty)$",fontsize=20)
ax[0].set_ylabel("$u_{\\theta}$ [$\\mathrm{ms^{-1}}$]",fontsize=16)
ax[0].set_xlim([ts[0],ts[-1]])
ax[0].set_ylim([-.1,.1])
ax[0].set_xlabel("Inertial Periods",fontsize=14)
ax[0].tick_params(axis="both",labelsize=16)
ax[0].text(1.02,0,"Wind off",rotation=90)
ax[0].text(-.55,.13,"a)",fontsize=16)
for jt in range(len(nt)):
    ax[1].fill_between((-1+.03,1), (-3,-3),(3,3),facecolor='black',zorder=10)
    d1 = ax[1].plot(g09["rstar_x"],np.ma.array(g09["ur"][:,jt],mask=(g09["mask"][:,jt])),color=cmap(norm(ts[nt[jt]])),zorder=5)
    sc = ax[1].plot(rstar,ur_scale,'--k')
    sc = ax[1].plot(rstar,-1*ur_scale,'--k',label="Theory")
    ax[1].text(.05,1.5,"$u_{r}$",fontsize=24,color="white",zorder=15)
    ax[1].set_ylim((-2,2))
    ax[1].set_xlim((0,10))
    ax[1].axhline(y=0,color="k",linewidth=.75)
    ax[1].tick_params(axis="both",labelsize=16)
    ax[1].set_ylabel("$\\frac{u_{r}}{\\vert u_{\\infty} \\vert}$",fontsize=24)
    ax[1].yaxis.set_major_formatter(FormatStrFormatter('%0.1f'))

    sc = ax[2].plot(g09["rstar_y"]-.03,np.ma.array(-g09["utheta"][:,jt],mask=(g09["mask"][:,jt])),color=cmap(norm(ts[nt[jt]])),zorder=5)
    sc = ax[2].plot(rstar,utheta_scale,'--k',label="Theory")
    sc = ax[2].plot(rstar,-1*utheta_scale,'--k',label="Theory")
    ax[2].text(.05,3.1,"$u_{\\theta}$",fontsize=24,color="white",zorder=15)
    ax[2].fill_between((-1+.03,1), (-4,-4),(4,4),facecolor='black',zorder=10)
    ax[2].set_ylim((-4,4))
    ax[2].set_xlim((0,10))
    ax[2].axhline(y=0,color="k")
    ax[2].axhline(y=0,color="k",linewidth=.5)
    ax[2].tick_params(axis="both",labelsize=16)
    ax[2].set_ylabel("$\\frac{u_{\\theta}}{\\vert u_{\\infty} \\vert}$",fontsize=24)
    ax[2].yaxis.set_major_formatter(FormatStrFormatter('%0.1f'))

    sc = ax[3].plot(g09["rstar_y"]-.03,np.ma.array(-g09["eta"][:,jt],mask=(g09["mask"][:,jt])),color=cmap(norm(ts[nt[jt]])),zorder=5)
    sc = ax[3].plot(rstar,zeta_scale,'--k',label="Theory")
    sc = ax[3].plot(rstar,-1*zeta_scale,'--k',label="Theory")
    ax[3].text(.05,3.1,"$\\eta$",fontsize=24,color="white",zorder=15)
    ax[3].fill_between((-1+.03,1), (-5,-5),(5,5),facecolor='black',zorder=10)
    ax[3].set_ylim((-5,5))
    ax[3].set_xlim((0,10))
    ax[3].set_xlabel("$\\frac{r}{a}$",fontsize=24)
    ax[3].axhline(y=0,color="k")
    ax[3].axhline(y=0,color="k",linewidth=.5)
    ax[3].tick_params(axis="both",labelsize=16)
    ax[3].set_ylabel("$\\frac{g\\eta}{a f\\vert u_{\\infty} \\vert}$",fontsize=24)
    ax[3].yaxis.set_major_formatter(FormatStrFormatter('%0.1f'))

box = ax[0].get_position()
box.x0 = box.x0
box.x1 = box.x1
box.y0 = box.y0 + .05
ax[0].set_position(box)
ax[1].text(-1,2.2,"b)",fontsize=16)
ax[2].text(-1,2.65,"c)",fontsize=16)
ax[3].text(-1,5.25,"d)",fontsize=16)

fig.savefig(figpath + svname, dpi=300, bbox_inches='tight')
