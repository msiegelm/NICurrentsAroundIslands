#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 14 10:33:17 2017

@author: msiegelman
"""
import numpy as np

def coriolis(lat):
	"""
	f = coriolis frequency (rad/sec)
	ft = coriolis frequency(cycle/day)
	T = coriolis period (hours/cycle)
	"""
	omega = (2*np.pi)/86400
	latrad = lat*np.pi/180
	f = 2*omega*np.sin(latrad)
	T = ((1/f)*2*np.pi)/60/60 # hours
	ft = 1/(T/24)
	return f,ft,T
    