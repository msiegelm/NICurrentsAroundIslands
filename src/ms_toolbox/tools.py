#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 24 09:12:09 2018

@author: msiegelman
"""

import numpy as np
import matplotlib.dates as dates
from pyproj import Proj

def num2dt64(tnum):
    """
    Converts matplotlib datenum to numpy datetime64.

    Input:
    -----
    tnum = datenum

    Output:
    ------
    t    = numpy datetime64
    """
    tdates = dates.num2date(tnum)
    t = [np.datetime64(tdates[ii].strftime("%Y-%m-%dT%H:%M:%S")) for ii in range(len(tdates))]

    return t


def vec_rot(u,v,ang):
    """
    Rotates vector(u,v) by angle (ang); + ang --> clockwise rotation

    Parameters
    ----------
    u: ndarray,
    v: ndarray,
    ang: float, degrees

    Returns
    -------
    urot: ndarray
    vrot: ndarray

    rotated vectors

    """


    uu = (u + v * 1j) * np.exp(-1j * ang * np.pi / 180)
    urot = np.real(uu)
    vrot = np.imag(uu)
    return urot, vrot

def pcm_xy(x,y):
	"""
	Makes x and y for pcolormesh to match with correct locations

	Parameters
	----------
	x: ndarray,
	y: ndarray,

	Returns
	-------
	xpc: ndarray
	ypc: ndarray
	x and y for pcolormesh

	"""
	dx = x[1]-x[0]
	dy = y[1]-y[0]
	xpc = np.hstack((x - (.5*dx),x[-1] + (.5*dx)))
	ypc = np.hstack((y - (.5*dy),y[-1] + (.5*dy)))

	return xpc,ypc

def interpolated_intercept(x, y1, y2):
    """Find the intercept of two curves, given by the same x data"""

    def intercept(point1, point2, point3, point4):
        """find the intersection between two lines
        the first line is defined by the line between point1 and point2
        the first line is defined by the line between point3 and point4
        each point is an (x,y) tuple.

        So, for example, you can find the intersection between
        intercept((0,0), (1,1), (0,1), (1,0)) = (0.5, 0.5)

        Returns: the intercept, in (x,y) format
        """

        def line(p1, p2):
            A = (p1[1] - p2[1])
            B = (p2[0] - p1[0])
            C = (p1[0]*p2[1] - p2[0]*p1[1])
            return A, B, -C

        def intersection(L1, L2):
            D  = L1[0] * L2[1] - L1[1] * L2[0]
            Dx = L1[2] * L2[1] - L1[1] * L2[2]
            Dy = L1[0] * L2[2] - L1[2] * L2[0]

            x = Dx / D
            y = Dy / D
            return x,y

        L1 = line([point1[0],point1[1]], [point2[0],point2[1]])
        L2 = line([point3[0],point3[1]], [point4[0],point4[1]])

        R = intersection(L1, L2)

        return R

    idx = np.squeeze(np.argwhere(np.diff(np.sign(y1 - y2)) != 0))
    # if np.any(np.shape(idx)):
    #     xc, yc = intercept((x[idx[0]], y1[idx[0]]),((x[idx[0]+1], y1[idx[0]+1])), ((x[idx[0]], y2[idx[0]])), ((x[idx[0]+1], y2[idx[0]+1])))
    # else:
    xc, yc = intercept((x[idx], y1[idx]),((x[idx+1], y1[idx+1])), ((x[idx], y2[idx])), ((x[idx+1], y2[idx+1])))
    kk = ~np.isnan(xc)
    if np.any(kk):
        return np.squeeze(xc[kk][0]),np.squeeze(yc[kk][0])
    else:
        return np.nan,np.nan


def convert_utmll(x,y,zone,hemisphere,toLatLon=False):
    """
    Converts between UTM coordinates and Lat/Lon using pyproj
    Input:
    -----
    x          = either UTMx or Lon
    y          = either UTMy or Lat
    zone       = WGS-84 Zone eg. for Palau zone = 53 (string or int)
    hemisphere = "north" or "south"
    toLatLon   = True if going from UTM to Lat/Lon

    Output:
    ------
    UTMx,UTMy or lon,lat if toLatLon = True

    """
    myProj = Proj('+proj=utm +zone=%s +%s +ellps=WGS84'%(zone,hemisphere), preserve_units=False)

    if toLatLon:
        lon,lat = myProj(x,y,inverse=True)
        return lon,lat
    else:
        UTMx, UTMy = myProj(x,y)
        return UTMx,UTMy
