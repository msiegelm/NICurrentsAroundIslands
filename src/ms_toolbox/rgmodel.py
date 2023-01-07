#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  4 10:29:18 2019

@author: msiegelman
"""
import time
import numpy as np
import matplotlib.pyplot as plt
from pycurrents.system import Bunch
from scipy.io import savemat
import netCDF4
from netCDF4 import Dataset

def writeNC(ncFname,h,u,v,t,n,taux,tauy,mm):

    kk = np.shape(h)

    ncout=Dataset(ncFname + '.nc','w',format='NETCDF4')
    ncout.createDimension('ntime',kk[0])
    ncout.createDimension('ny',kk[1])
    ncout.createDimension('nx',kk[2])

    # define variables for output to netcdf
    hh=ncout.createVariable('hh',float,('ntime','ny','nx'),zlib=True)
    uu=ncout.createVariable('uu',float,('ntime','ny','nx'),zlib=True)
    vv=ncout.createVariable('vv',float,('ntime','ny','nx'),zlib=True)
    tt=ncout.createVariable('tt',float,('ntime'),zlib=True)
    nn=ncout.createVariable('nn',float,('ntime'),zlib=True)
    tauxx=ncout.createVariable('taux',float,('ntime'),zlib=True)
    tauyy=ncout.createVariable('tauy',float,('ntime'),zlib=True)




            # define attributes that correspond to each variable
    hh.long_name='Depth of thermocline [m]'
    hh.units='meters'
    hh.standard_name='thermocline depth'
    hh.field='h, scalar'
    hh.coordinates='ntime ny nx'

    uu.long_name='u-velocity (E-W velocity) [ms^{-1}]'
    uu.units='meters per second'
    uu.standard_name='u-velocity'
    uu.field='u, scalar'
    uu.coordinates='ntime ny nx'

    vv.long_name='v-velocity (N-S velocity) [ms^{-1}]'
    vv.units='meters per second'
    vv.standard_name='v-velocity'
    vv.field='v, scalar'
    vv.coordinates='ntime ny nx'

    tauxx.long_name='x-component, Wind Stress  [N m^-2]'
    tauxx.units='seconds'
    tauxx.standard_name='tau-x'
    tauxx.field='t, scalar'
    tauxx.coordinates='ntime'

    tauyy.long_name='y-component, Wind Stress  [N m^-2]'
    tauyy.units='seconds'
    tauyy.standard_name='tau-y'
    tauyy.field='t, scalar'
    tauyy.coordinates='ntime'

    tt.long_name='time [seconds]'
    tt.units='seconds'
    tt.standard_name='time'
    tt.field='t, scalar'
    tt.coordinates='ntime'

    tt.long_name='time [seconds]'
    tt.units='seconds'
    tt.standard_name='time'
    tt.field='t, scalar'
    tt.coordinates='ntime'

    nn.long_name='number of time steps'
    nn.units='number'
    nn.standard_name='number'
    nn.field='t, scalar'
    nn.coordinates='ntime'


    ncout.variables['hh'][:]=np.ma.array(h)
    ncout.variables['uu'][:]=np.ma.array(u)
    ncout.variables['vv'][:]=np.ma.array(v)
    ncout.variables['taux'][:]=np.ma.array(taux)
    ncout.variables['tauy'][:]=np.ma.array(tauy)
    ncout.variables['tt'][:]=np.ma.array(t)
    ncout.variables['nn'][:]=np.ma.array(n)

    ncout.sync()

    ncout.close()


def writeNC_grid(ncFname,h,mm,xrad,yrad):

    kk = np.shape(h)

    ncout=Dataset(ncFname + '.nc','w',format='NETCDF4')

    # Global Attributes

    ncout.grid_bc = mm.grid_bc
    ncout.is_bc = mm.is_bc

    # Create Dimensions
    ncout.createDimension('ny',kk[0])
    ncout.createDimension('nx',kk[1])
    ncout.createDimension('sing',1)

    # define variables for output to netcdf
    Ah=ncout.createVariable('Ah',float,('sing'),zlib=True)
    ff=ncout.createVariable('f',float,('ny'),zlib=True)
    hh=ncout.createVariable('H',float,('sing'),zlib=True)
    xx=ncout.createVariable('x',float,('nx'),zlib=True)
    yy=ncout.createVariable('y',float,('ny'),zlib=True)
    lon=ncout.createVariable('lon',float,('nx'),zlib=True)
    lat=ncout.createVariable('lat',float,('ny'),zlib=True)
    el=ncout.createVariable('el',float,('ny','nx'),zlib=True)

    mask_h=ncout.createVariable('mask_h',float,('ny','nx'),zlib=True)
    mask_u=ncout.createVariable('mask_u',float,('ny','nx'),zlib=True)
    mask_v=ncout.createVariable('mask_v',float,('ny','nx'),zlib=True)
    Xrad=ncout.createVariable('Xrad',float,('sing'),zlib=True)
    Yrad=ncout.createVariable('Yrad',float,('sing'),zlib=True)

            # define attributes that correspond to each variable
    Ah.long_name='Viscosity [m^2 s^-1]'
    Ah.units='meters^2 s^-1'
    Ah.standard_name='Horizontal Eddy Viscosity'
    Ah.field='Ah, scalar'
    Ah.coordinates='sing'

    ff.long_name='Coriolis parameter [radians s^-1]'
    ff.units='radians per second'
    ff.standard_name='Coriolis parameter'
    ff.field='f, scalar'
    ff.coordinates='ny'

    hh.long_name='Layer thickness [m]'
    hh.units='meters'
    hh.standard_name='Layer thickness'
    hh.field='H, scalar'
    hh.coordinates='sing'

    xx.long_name='x [meters]'
    xx.units='meters'
    xx.standard_name='x'
    xx.field='x, scalar'
    xx.coordinates='nx'

    yy.long_name='y [meters]'
    yy.units='meters'
    yy.standard_name='y'
    yy.field='y, scalar'
    yy.coordinates='ny'

    lon.long_name='Longitude [degrees]'
    lon.units='degrees'
    lon.standard_name='Longitude'
    lon.field='lon, scalar'
    lon.coordinates='nx'

    lat.long_name='Latitude [degrees]'
    lat.units='Degrees'
    lat.standard_name='Latitude'
    lat.field='lat, scalar'
    lat.coordinates='ny'

    el.long_name='Topography, Mask for h-grid (Arakawa C-grid)'
    el.units='Binary (0= masked)'
    el.standard_name='Topography'
    el.field='el, scalar'
    el.coordinates='ny nx'

    mask_h.long_name='Mask for h-grid (Arakawa C-grid)'
    mask_h.units='Binary (0= masked)'
    mask_h.standard_name='Mask h-grid'
    mask_h.field='mask_h, scalar'
    mask_h.coordinates='ny nx'

    mask_u.long_name='Mask for u-grid (Arakawa C-grid)'
    mask_u.units='Binary (0= masked)'
    mask_u.standard_name='Mask u-grid'
    mask_u.field='mask_u, scalar'
    mask_u.coordinates='ny nx'

    mask_v.long_name='Mask for v-grid (Arakawa C-grid)'
    mask_v.units='Binary (0= masked)'
    mask_v.standard_name='Mask v-grid'
    mask_v.field='mask_v, scalar'
    mask_v.coordinates='ny nx'

    Xrad.long_name='x-axis length (a)'
    Xrad.units='meters'
    Xrad.standard_name='X Radius a'
    Xrad.field='radius, scalar'
    Xrad.coordinates='sing'

    Yrad.long_name='y-axis length (b)'
    Yrad.units='meters'
    Yrad.standard_name='Y Radius b'
    Yrad.field='radius, scalar'
    Yrad.coordinates='sing'


    ncout.variables['Ah'][:]=mm.Ah
    ncout.variables['f'][:]=np.squeeze(mm.f)
    ncout.variables['H'][:]=mm.H
    ncout.variables['x'][:]=mm.x
    ncout.variables['y'][:]=mm.y
    ncout.variables['lon'][:]=mm.phi
    ncout.variables['lat'][:]=mm.theta
    ncout.variables['el'][:]=mm.el
    ncout.variables['mask_h'][:]=mm.mask_h
    ncout.variables['mask_u'][:]=mm.mask_u
    ncout.variables['mask_v'][:]=mm.mask_v
    ncout.variables['Xrad'][:]=xrad
    ncout.variables['Yrad'][:]=yrad



    ncout.sync()

    ncout.close()


def int_saveNC(m,nint,ncFname,ndiv,xrad,yrad):

    h_f = []
    u_f = []
    v_f = []
    t_f = []
    n_f = []
    taux_f = []
    tauy_f = []

    for n,t in m.integrate():
        if not np.mod(n,ndiv):
            if n < nint+20:
                    h_f.append(m.h[:,:,m.nnext].copy())
                    u_f.append(m.u[:,:].copy())
                    v_f.append(m.v[:,:].copy())
                    t_f.append(t)
                    n_f.append(n)
                    taux_f.append(m.taux[n])
                    tauy_f.append(m.tauy[n])
                    print(n)
                    if n == nint:
                        writeNC(ncFname,h_f,u_f,v_f,t_f,n_f,taux_f,tauy_f,m)

                        gn = ncFname + '_grid'
                        writeNC_grid(gn,h_f,m,xrad,yrad)
            else:
                break


def int_return(m,nint,ndiv):

    h_f = []
    u_f = []
    v_f = []
    t_f = []
    n_f = []

    for n,t in m.integrate():
        if not np.mod(n,ndiv):
            if n < nint+20:
                    h_f.append(m.h[:,:,m.nnext].copy())
                    u_f.append(m.u[:,:].copy())
                    v_f.append(m.v[:,:].copy())
                    t_f.append(t)
                    n_f.append(n)
                    print(n)
                    if n == nint:
                        return h_f , u_f, v_f, t_f, n_f
            else:
                break


def animator(prob1,hh,uu,vv,t,Ny,nr1,nr2,linevar):

    #s=np.s_[rng[0]:rng[-1]+1]
    ss=np.s_[::3]

    nmidy = int(np.floor(Ny/2))

    fig, (ax0,ax1)  = plt.subplots(nrows = 2,ncols = 1, figsize=(10,7.5),gridspec_kw = {'height_ratios':[3, 1]})
    line,=ax1.plot(prob1.phi,prob1.phi*0,lw=2)


    for nt in range(nr1,nr2):
        ax0.cla()
        p=ax0.pcolormesh(prob1.phi,prob1.theta,np.ma.array(hh[nt,:,:],mask=prob1.mask_h==0),vmin=-.1,vmax=.1)
        ax0.quiver(prob1.phi[ss],prob1.theta[ss],
                              np.ma.array(uu[nt,ss,ss],mask=(prob1.mask_u[ss,ss]==0)),
                              np.ma.array(vv[nt,ss,ss],mask=(prob1.mask_v[ss,ss]==0)),scale=1)
        ax0.set_title("[%d]" % (nt))
        ax0.set_aspect("equal")

        if linevar is "h":
            line.set_data(prob1.phi,hh[nt,nmidy,:])

        elif linevar is "u":
            line.set_data(prob1.phi,uu[nt,nmidy,:])

        elif linevar is "v":
            line.set_data(prob1.phi,vv[nt,nmidy,:])

        plt.pause(.2)



def vorticity_calc(prob1,uu,vv,x,y):
    dx = x[1]-x[0]
    dy = y[1]-y[0]

    dudy = (uu[:,1:,:] - uu[:,:-1,:]) / dy

    dvdx = (vv[:,:,1:] - vv[:,:,:-1]) / dx

    dudymn = .5 * ( dudy[:,:,1:] + dudy[:,:,:-1])

    dvdxmn = .5 * ( dvdx[:,1:,:] + dvdx[:,:-1,:])

    vort = dvdxmn - dudymn

    xmn = .5 * (x[1:] + x[:-1])

    ymn = .5 * (y[1:] + y[:-1])


    return vort, xmn , ymn


def circum(Xrad,Yrad,shape):
    if shape == "circle":
        cir = 2 * np.pi * Xrad

    elif shape == "ellipse":

        cir = 4 * (Xrad + Yrad) * (np.pi/4)**((4*Xrad*Yrad)/(Xrad + Yrad)**2)

    return cir

def min_dist(xgrid,ygrid,xcoord,ycoord):
    dismat = np.sqrt((xgrid-xcoord)**2 + (ygrid-ycoord)**2)
    minin = np.unravel_index(np.nanargmin(dismat), dismat.shape)
    return minin



def vel_polar2cart(ur,utheta,th):
    u = (ur * np.cos(th)) - (utheta * np.sin(th))
    v = (ur * np.sin(th)) + (utheta * np.cos(th))
    return u, v

def vel_cart2polar(u,v,xg,yg):
    """
    Converts vectors in cartesian coordinate system to Polar coordinates
    Input: u,v in cartesian coordinates
           xg,yg is cartesian grid to get theta

    Output: ur,utheta polar coordinate vectors
    """

    Z = xg + 1j * yg
    th = np.angle(Z)
    ur = (u * np.cos(th)) + (v * np.sin(th))
    utheta = (-u * np.sin(th)) + (v * np.cos(th))
    return ur, utheta



def writeNC_slow(ncFname,m,nint,ndiv,xrad,yrad):

    ncout=Dataset(ncFname + '.nc','w',format='NETCDF4')
    ncout.createDimension('ntime',size=None)
    ncout.createDimension('ny',m.Ny)
    ncout.createDimension('nx',m.Nx)

    # define variables for output to netcdf
    hh=ncout.createVariable('hh',np.float32,('ntime','ny','nx'))
    uu=ncout.createVariable('uu',np.float32,('ntime','ny','nx'))
    vv=ncout.createVariable('vv',np.float32,('ntime','ny','nx'))
    tt=ncout.createVariable('tt',float,('ntime'))
    nn=ncout.createVariable('nn',np.int16,('ntime'))
    tauxx=ncout.createVariable('taux',np.float32,('ntime'))
    tauyy=ncout.createVariable('tauy',np.float32,('ntime'))


    # define attributes that correspond to each variable
    hh.long_name='Depth of thermocline [m]'
    hh.units='meters'
    hh.standard_name='thermocline depth'
    hh.field='h, scalar'
    hh.coordinates='ntime ny nx'

    uu.long_name='u-velocity (E-W velocity) [ms^{-1}]'
    uu.units='meters per second'
    uu.standard_name='u-velocity'
    uu.field='u, scalar'
    uu.coordinates='ntime ny nx'

    vv.long_name='v-velocity (N-S velocity) [ms^{-1}]'
    vv.units='meters per second'
    vv.standard_name='v-velocity'
    vv.field='v, scalar'
    vv.coordinates='ntime ny nx'

    tt.long_name='time [seconds]'
    tt.units='seconds'
    tt.standard_name='time'
    tt.field='t, scalar'
    tt.coordinates='ntime'

    nn.long_name='number of time steps'
    nn.units='number'
    nn.standard_name='number'
    nn.field='t, scalar'
    nn.coordinates='ntime'

    tauxx.long_name='x-component, Wind Stress  [N m^-2]'
    tauxx.units='Newtons per meters squared'
    tauxx.standard_name='tau-x'
    tauxx.field='t, scalar'
    tauxx.coordinates='ntime'

    tauyy.long_name='y-component, Wind Stress  [N m^-2]'
    tauyy.units='Newtons per meters squared'
    tauyy.standard_name='tau-y'
    tauyy.field='t, scalar'
    tauyy.coordinates='ntime'

    toc = time.time()
    ind = 0
    for n,t in m.integrate():
        if not np.mod(n,ndiv):
            print(t)
            tic = time.time()
            int_time = tic - toc
            ncout.variables['hh'][ind,:,:]=m.h[:,:,m.nnext]
            ncout.variables['uu'][ind,:,:]=m.u[:,:]
            ncout.variables['vv'][ind,:,:]=m.v[:,:]
            ncout.variables['taux'][ind]=m.taux[n]
            ncout.variables['tauy'][ind]=m.tauy[n]
            ncout.variables['tt'][ind]=np.array(t)
            ncout.variables['nn'][ind]=np.array(n)
            toc = time.time()
            print("%d %.3f %.3f" % (n, int_time, toc-tic))
            ind+=1
            print(ind)
            print(nint)
        if n == nint:
            ncout.close()
            gn = ncFname + '_grid'
            h_f = m.h[:,:,m.nnext]
            print(type(h_f))
            writeNC_grid(gn,h_f,m,xrad,yrad)
            break

def make_ellipse(x,y,a,b,xc,yc):
    """
    Input:
     x = x-grid
     y = y-grid
     a = x radius
     b = y radius
     xc = x index of middle of ellipse on grid
     yc = y index of middle of ellipse on grid

     Return:
     el = Masked grid with ellipse

    """
    xx, yy = np.meshgrid(x,y)
    el = ((xx - x[xc]) ** 2 / a**2) + ((yy - y[yc]) ** 2 / b**2)

    kk = (el <= 1)

    el[kk] = 0
    el[~kk] = 1


    # Remove unmask rows and columns that have only one masked element

    nmsk = (el == 0)

    nzero = nmsk.sum(axis = 0)

    kk = (nzero == 1)

    el[:,kk] = 1


    nzero = nmsk.sum(axis = 1)

    kk = (nzero == 1)

    el[kk,:] = 1

    return el


def make_wind(nevent,omega, amp, times,nstart,taper=True,wind_type="CW"):

    """
    nevent = number of cycles for event
    omega = angular freqency of wind event (units of radians s ** -1)
    times = time of model run
    nstart = number of time steps to delay start of wind event
    taper = True = Taper wind event
    wind_type = type of wind event (CW = CW rotating, constant_x, constant_y)
    """

    freq = omega / (2 * np.pi)

    T = 1/freq

    nsec_event = nevent * T

    tstart = times[0]
    dt = times[1]-times[0]
    tevent = np.arange(tstart,nsec_event + dt, dt)


    if taper:
        w = np.blackman(len(tevent))

    else:
        w = np.ones(len(tevent))

    if wind_type is "CW":
        WW = amp * np.exp(-1j * omega * tevent)
        taux_ev = w * np.real(WW)
        tauy_ev = w * np.imag(WW)

    elif wind_type is "constant_x":

        taux_ev = amp * np.ones(len(tevent))
        tauy_ev = np.zeros(len(tevent))

    elif wind_type is "constant_y":
        taux_ev = np.zeros(len(tevent))
        tauy_ev = amp * np.ones(len(tevent))


    sl_ = slice(nstart,nstart+len(taux_ev))

    lent = len(times) + len(taux_ev)

    taux = np.zeros(lent)
    tauy = np.zeros(lent)

    taux[sl_] = taux_ev
    tauy[sl_] = tauy_ev

    taux = taux[:len(times)]
    tauy = tauy[:len(times)]
    return taux, tauy


def read_grid(xrng,yrng,ncFnameG,skip=1):
    """
    Input:
    xrng,yrng = tuple of min and max x,y
    ncFnameG = path + file name to read grid

    output:
    g structure
    """
    dnc_grid = Dataset(ncFnameG)
    g = Bunch()
    xln = len(xrng)
    if xln == 0 and skip == 1:

        g["x"] = dnc_grid['x'][:]
        g["y"] = dnc_grid['y'][:]
        g["f"] = dnc_grid['f'][:]
        g["Ah"] = dnc_grid['Ah'][:]
        g["H"] = dnc_grid['H'][:]
        g["phi"] = dnc_grid['lon'][:]
        g["theta"] = dnc_grid['lat'][:]
        g["mask_h"] = dnc_grid['mask_h'][:,:]
        g["mask_u"] = dnc_grid['mask_u'][:,:]
        g["mask_v"] = dnc_grid['mask_v'][:,:]
        g["Xrad"] = dnc_grid['Xrad'][:]
        g["Yrad"] = dnc_grid['Yrad'][:]
        g["el"] = dnc_grid['el'][:,:]

    elif xln == 0 and skip != 1:
        sl = np.s_[::skip]
        print("hi")
        g["x"] = dnc_grid['x'][sl]
        print("1")
        g["y"] = dnc_grid['y'][sl]
        print("2")
        g["f"] = dnc_grid['f'][sl]
        print("3")
        g["Ah"] = dnc_grid['Ah'][:]
        print("4")
        g["H"] = dnc_grid['H'][:]
        print("5")
        g["phi"] = dnc_grid['lon'][sl]
        print("6")
        g["theta"] = dnc_grid['lat'][sl]
        print("7")
        g["mask_h"] = dnc_grid['mask_h'][sl,sl]
        print("8")
        g["mask_u"] = dnc_grid['mask_u'][sl,sl]
        print("9")
        g["mask_v"] = dnc_grid['mask_v'][sl,sl]
        print("10")
        g["Xrad"] = dnc_grid['Xrad'][:]
        print("11")
        g["Yrad"] = dnc_grid['Yrad'][:]
        print("12")
        g["el"] = dnc_grid['el'][sl,sl]
        print("13")

    else:
        xsl = slice(xrng[0],xrng[1])
        ysl = slice(yrng[0],yrng[1])

        g["x"] = dnc_grid['x'][xsl]
        g["y"] = dnc_grid['y'][ysl]
        g["f"] = dnc_grid['f'][ysl]
        g["Ah"] = dnc_grid['Ah'][:]
        g["H"] = dnc_grid['H'][:]
        g["phi"] = dnc_grid['lon'][xsl]
        g["theta"] = dnc_grid['lat'][ysl]
        g["mask_h"] = dnc_grid['mask_h'][ysl,xsl]
        g["mask_u"] = dnc_grid['mask_u'][ysl,xsl]
        g["mask_v"] = dnc_grid['mask_v'][ysl,xsl]
        g["Xrad"] = dnc_grid['Xrad'][:]
        g["Yrad"] = dnc_grid['Yrad'][:]
        g["el"] = dnc_grid['el'][ysl,xsl]

    return g

def read_output(xrng,yrng,ncFname,skip=1):
    """
    Input:
    xrng,yrng = tuple of min and max x,y
    ncFname = path + file name to read model output

    output:
    g structure
    """
    dnc = Dataset(ncFname)
    opns = Bunch()
    xln = len(xrng)
    if xln == 0 and skip == 1:
        opns["hh"] = dnc["hh"][:,:,:]
        opns["uu"] = dnc["uu"][:,:,:]
        opns["vv"] = dnc["vv"][:,:,:]
        opns["tt"] = dnc["tt"][:]
        opns["taux"] = dnc["taux"][:]
        opns["tauy"] = dnc["tauy"][:]
        opns["tdays"] = opns["tt"]/86400

    elif xln == 0 and skip != 1:

        sl = np.s_[::skip]
        opns["hh"] = dnc["hh"][:,sl,sl]
        opns["uu"] = dnc["uu"][:,sl,sl]
        opns["vv"] = dnc["vv"][:,sl,sl]
        opns["tt"] = dnc["tt"][:]
        opns["taux"] = dnc["taux"][:]
        opns["tauy"] = dnc["tauy"][:]
        opns["tdays"] = opns["tt"]/86400

    else:

        xsl = slice(xrng[0],xrng[1])
        ysl = slice(yrng[0],yrng[1])
        opns["hh"] = dnc["hh"][:,ysl,xsl]
        opns["uu"] = dnc["uu"][:,ysl,xsl]
        opns["vv"] = dnc["vv"][:,ysl,xsl]
        opns["tt"] = dnc["tt"][:]
        opns["taux"] = dnc["taux"][:]
        opns["tauy"] = dnc["tauy"][:]
        opns["tdays"] = opns["tt"]/86400

    return opns


def read_output_ts(trng,xrng,yrng,ncFname):
    """
    Input:
    xrng,yrng = tuple of min and max x,y
    ncFname = path + file name to read model output

    output:
    g structure
    """
    dnc = Dataset(ncFname)
    opns = Bunch()
    xln = len(xrng)
    if xln == 0:
        tsl = slice(trng[0],trng[1])
        opns["hh"] = dnc["hh"][tsl,:,:]
        opns["uu"] = dnc["uu"][tsl,:,:]
        opns["vv"] = dnc["vv"][tsl,:,:]
        opns["tt"] = dnc["tt"][tsl]
        opns["taux"] = dnc["taux"][tsl]
        opns["tauy"] = dnc["tauy"][tsl]
        opns["tdays"] = opns["tt"]/86400

    else:
        tsl = slice(trng[0],trng[1])
        xsl = slice(xrng[0],xrng[1])
        ysl = slice(yrng[0],yrng[1])
        opns["hh"] = dnc["hh"][tsl,ysl,xsl]
        opns["uu"] = dnc["uu"][tsl,ysl,xsl]
        opns["vv"] = dnc["vv"][tsl,ysl,xsl]
        opns["tt"] = dnc["tt"][tsl]
        opns["taux"] = dnc["taux"][tsl]
        opns["tauy"] = dnc["tauy"][tsl]
        opns["tdays"] = opns["tt"]/86400

    return opns

def test_Econserve(mo,save=False,figpath=[]):

    gg = 9.8
    rho = 1000
    dt = mo.tt[1]-mo.tt[0]
    ke = .5 * rho * (mo.uu**2 + mo.vv**2)
    pe = .5 * rho * gg * mo.eta**2
    E = pe + ke

    kets = np.sum(ke,axis=(1,2))
    pets = np.sum(pe,axis=(1,2))
    Ets = kets + pets


    fig,ax = plt.subplots(ncols=1,nrows=3,figsize=(10,10))
    ax[0].plot(mo.tdays,kets)
    ax[0].set_title("Kinetic Energy")
    ax[1].plot(mo.tdays,pets)
    ax[1].set_title("Potential Energy")
    ax[2].plot(mo.tdays,Ets)
    ax[2].set_title("Total Energy")

    if save:
        fig.savefig(figpath + 'total_energy' + ".png", dpi=300, bbox_inches='tight',orientation="horizontal")

    return kets, pets, Ets
