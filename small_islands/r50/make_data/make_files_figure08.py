#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan  5 10:34:16 2020

@author: miks
"""
import sys
apath = "/Users/miks/Desktop/Palau/reduced_gravity_model/rg_tools/"
sys.path.append(apath)
import rg_plotting as rgp
import animate_plotting as ap
import a_tools as at

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

"""

"""
#%%

################# Define Paths ###################

fpath = '/Users/miks/Desktop/Palau/reduced_gravity_model/Runs/circle_PW/base/'
svpath = fpath + 'output/'

ncname = 'C_LR_FS_r50_Ah00_r0_gp09'

ncFname = svpath + ncname + '.nc'
ncFnameG = svpath + ncname + '_grid.nc'

figpath = fpath + 'figures/'   #path to save figures to

anpath = fpath + 'files_for_publication/data_subsection/'

####################################################


################# Read Grid ###################
# xrng = (3,1599)
# yrng = (3,1599)
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

######## Check to confirm output for TS ##########

# fig,ax = plt.subplots()
# pc = ax.pcolormesh(g.x/1000,g.y/1000,g.el)
# ax.set_aspect("equal")
# ax.axes.set_xlabel("km")
# ax.axes.set_ylabel("km")
#################################################

######## Find closest to point to make TS ##########
mindist_n = rgm.min_dist(g.xx,g.yy,3000e3,3050e3)
mindist_s = rgm.min_dist(g.xx,g.yy,3000e3,2950e3)
mindist_e = rgm.min_dist(g.xx,g.yy,3050e3,3000e3)
mindist_w = rgm.min_dist(g.xx,g.yy,2950e3,3000e3)
mindist_ff = rgm.min_dist(g.xx,g.yy,5900e3,5900e3)
mindist_ff2 = rgm.min_dist(g.xx,g.yy,3000e3,5000e3)

ele_n = [mindist_n[0],mindist_n[1]] #y,x
ele_e = [mindist_e[0],mindist_e[1]] #y,x
ele_s = [mindist_s[0],mindist_s[1]] #y,x
ele_w = [mindist_w[0],mindist_w[1]] #y,x
ele_ff = [mindist_ff[0],mindist_ff[1]] #y,x
ele_ff2 = [mindist_ff2[0],mindist_ff2[1]] #y,x

###################################################


def find_extrema(ts):
    du = ts[1:] - ts[:-1]
    kkp =  np.where(du > 0)
    kkm =  np.where(du < 0)


    kk_f = np.ones_like(du)
    kk_f[kkp] = 1
    kk_f[kkm] = 0

    dkk_f = np.diff(kk_f)

    kkp1 = np.where(dkk_f == 1)
    kkm1 = np.where(dkk_f == -1)


    extrema = Bunch()
    extrema.max = kkm1[0] + 1
    extrema.min = kkp1[0] + 1

    eall = np.hstack((extrema.max,extrema.min))
    extrema.all = np.sort(eall)

    return extrema


########################################################
#%% Math coordinates

Ny = len(g.y)
Nx = len(g.x)
indy = int(np.floor(Ny/2))
indx = int(np.floor(Nx/2))

xmath = g.x - g.x[indx]
ymath = g.y - g.y[indy]

xmath_g,ymath_g = np.meshgrid(xmath,ymath)


#%%


mo.ur = np.ones_like(mo.uu)
mo.utheta = np.ones_like(mo.uu)


for nt in range(0,len(mo.tt)):
   ur,utheta = rgm.vel_cart2polar(mo.uu[nt,:,:],mo.vv[nt,:,:],xmath_g,ymath_g)
   mo.ur[nt,:,:] = ur
   mo.utheta[nt,:,:] = utheta

#%%


#############################################%%
Ny = len(g.y)
Nx = len(g.x)
indy = int(np.floor(Ny/2))
indx = int(np.floor(Nx/2))

xmath = g.x - g.x[indx]
ymath = g.y - g.y[indy]

xmath_g,ymath_g = np.meshgrid(xmath,ymath)


ur_tran_y = np.squeeze(mo.ur[:,:,indx])
utheta_tran_y = np.squeeze(mo.utheta[:,:,indx])
eta_tran_y = np.squeeze(mo.eta[:,:,indx])

ur_tran_x = np.squeeze(mo.ur[:,indy,:])
utheta_tran_x = np.squeeze(mo.utheta[:,indy,:])
eta_tran_x = np.squeeze(mo.eta[:,indy,:])

# X-transect Negative infinity
ur_infn_x = ur_tran_x[:,4]
utheta_infn_x = utheta_tran_x[:,4]
eta_infn_x = eta_tran_x[:,4]
mag_infn_x = np.abs(ur_infn_x + 1j* utheta_infn_x)


# Y-transect Negative infinity
ur_infn_y = ur_tran_y[:,4]
utheta_infn_y = utheta_tran_y[:,4]
eta_infn_y = eta_tran_y[:,4]
mag_infn_y = np.abs(ur_infn_y + 1j* utheta_infn_y)

# X-transect Positive infinity
ur_infp_x = ur_tran_x[:,-5]
utheta_infp_x = utheta_tran_x[:,-5]
eta_infp_x = eta_tran_x[:,-5]
mag_infp_x = np.abs(ur_infp_x + 1j* utheta_infp_x)

# Y-transect Positive infinity
ur_infp_y = ur_tran_y[:,-5]
utheta_infp_y = utheta_tran_y[:,-5]
eta_infp_y = eta_tran_y[:,-5]
mag_infp_y = np.abs(ur_infp_y + 1j* utheta_infp_y)




a = g.Xrad


kn_y = (ymath<0)
kp_y = (ymath>0)

kn_x = (xmath<0)
kp_x = (xmath>0)

rstar_y = ymath/a
rstar_x = xmath/a

mask_tran = g.mask_h

ur_ratio_x = np.ones_like(ur_tran_x)
utheta_ratio_x = np.ones_like(utheta_tran_x)
eta_ratio_x = np.ones_like(eta_tran_x)

ur_ratio_y = np.ones_like(ur_tran_y)
utheta_ratio_y = np.ones_like(utheta_tran_y)
eta_ratio_y = np.ones_like(eta_tran_y)


## Normalize by velocity magnitude
ur_ratio_x[:,kn_x] = (ur_tran_x[:,kn_x].T/mag_infn_x).T
utheta_ratio_x[:,kn_x] = (utheta_tran_x[:,kn_x].T/mag_infn_x).T
eta_ratio_x[:,kn_x] = ((eta_tran_x[:,kn_x].T * 9.8) / (a * g.fcor[0] * mag_infn_x)).T

ur_ratio_x[:,kp_x] = (ur_tran_x[:,kp_x].T/mag_infp_x).T
utheta_ratio_x[:,kp_x] = (utheta_tran_x[:,kp_x].T/mag_infp_x).T
eta_ratio_x[:,kp_x] = ((eta_tran_x[:,kp_x].T * 9.8) / (a * g.fcor[0] * mag_infp_x)).T


ur_ratio_y[:,kn_y] = (ur_tran_y[:,kn_y].T/mag_infn_y).T
utheta_ratio_y[:,kn_y] = (utheta_tran_y[:,kn_y].T/mag_infn_y).T
eta_ratio_y[:,kn_y] = ((eta_tran_y[:,kn_y].T * 9.8) / (a * g.fcor[0] * mag_infn_y)).T

ur_ratio_y[:,kp_y] = (ur_tran_y[:,kp_y].T/mag_infp_y).T
utheta_ratio_y[:,kp_y] = (utheta_tran_y[:,kp_y].T/mag_infp_y).T
eta_ratio_y[:,kp_y] = ((eta_tran_y[:,kp_y].T * 9.8) / (a * g.fcor[0] * mag_infp_y)).T






ext = find_extrema(utheta_infp_y)



nt = ext.all[1:]
ntss = ext.all[2]


ymath_g = np.ones_like(utheta_ratio_y.T)*ymath[:,np.newaxis]
t_gy = np.ones_like(utheta_ratio_y.T)* mo.tdays[np.newaxis,:]

xmath_g = np.ones_like(utheta_ratio_x.T)*xmath[:,np.newaxis]
t_gx = np.ones_like(utheta_ratio_x.T)* mo.tdays[np.newaxis,:]




nt = ext.all[:]
utheta_norm = utheta_ratio_y[nt,:].T
eta_norm = eta_ratio_y[nt,:].T
ur_norm = ur_ratio_x[nt,:].T
#
nmid = int(np.floor(Nx/2))

ma = g.mask_h[:,nmid]
ma = ma[np.newaxis,:]
ma_rep = np.repeat(ma,len(nt),axis=0)

utheta_norm.mask = (ma_rep.T == 0)
ur_norm.mask = (ma_rep.T == 0)
eta_norm.mask = (ma_rep.T == 0)



mask = (ma_rep.T == 0)

tnt = mo.tdays[nt]
uu = utheta_infp_y

np.savez(anpath + "normalization_snapshot_FINAL", mask = mask, rstar_x = rstar_x, rstar_y = rstar_y, utheta = utheta_norm, eta=eta_norm,ur=ur_norm,t=tnt,nt=nt,td=mo.tinertial,uu=uu)
