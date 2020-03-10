#!j/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue March 10 2020

@author: benhills
"""

import numpy as np

def principal_strain(xc,yc,xs,ys,dxs,dys,dx_errs=None,dy_errs=None,mindist=0.,verbose=True):
    """
    Calculate the principal strain tensor from velocity vectors

    Parameters
    ---------
    xc:     x center point about which to calculate the strain rates (m)
    yc:     y center point about which to calculate the strain rates (m)
    xs:     x points where the velocities are measured (m)
    ys:     y points where the velocities are measured (m)
    dxs:    x velocities (m/yr)
    dys:    y velocities (m/yr)
    dx_errs:    x velocity error (m/yr)
    dy_errs:    y velocity error (m/yr)

    Output
    ---------
    theta_p:    principal angle from horizontal (rad)
    eps_11:     first principal strain (along theta_p) (yr-1)
    eps_22:     second principal strain (along theta_p + np.pi/2) (yr-1)
    eps_err:    strain rate error (yr-1)

    """

    if len(xs) <= 1:
        if verbose:
            print('Only one vector input to the strain rate calculation.')
        return np.nan, np.nan, np.nan, np.nan

    # Distance from center to each point
    xdists = (xs-xc)
    ydists = (ys-yc)

    # Strain rates oriented with x-y
    eps_yy = np.nanmean((dys[abs(ydists)>mindist]-np.nanmean(dys[abs(ydists)>mindist]))/ydists[abs(ydists)>mindist])
    eps_xy = np.nanmean((dxs[abs(ydists)>mindist]-np.nanmean(dxs[abs(ydists)>mindist]))/ydists[abs(ydists)>mindist])
    eps_yx = np.nanmean((dys[abs(xdists)>mindist]-np.nanmean(dys[abs(xdists)>mindist]))/xdists[abs(xdists)>mindist])
    eps_xx = np.nanmean((dxs[abs(xdists)>mindist]-np.nanmean(dxs[abs(xdists)>mindist]))/xdists[abs(xdists)>mindist])

    # Strain rate error oriented with x-y
    if dx_errs is not None and dy_errs is not None:
        eps_yy_err = np.nanmean(dy_errs[abs(ydists)>mindist]/ydists[abs(ydists)>mindist])
        eps_xy_err = np.nanmean(dx_errs[abs(ydists)>mindist]/ydists[abs(ydists)>mindist])
        eps_yx_err = np.nanmean(dy_errs[abs(xdists)>mindist]/xdists[abs(xdists)>mindist])
        eps_xx_err = np.nanmean(dx_errs[abs(xdists)>mindist]/xdists[abs(xdists)>mindist])

    # Principal angle
    theta_p = np.arctan(2.*np.mean([eps_xy,eps_yx])/(eps_xx-eps_yy))/2.
    if (eps_xx-eps_yy) == 0. and np.mean([eps_xy,eps_yx]) == 0.:
        theta_p = 0.

    # Principal strain rates
    M1 = np.array([[np.cos(theta_p),np.sin(theta_p)],[-np.sin(theta_p),np.cos(theta_p)]])
    M2 = np.array([[eps_xx,np.mean([eps_xy,eps_yx])],[np.mean([eps_xy,eps_yx]),eps_yy]])
    M3 = np.array([[np.cos(theta_p),-np.sin(theta_p)],[np.sin(theta_p),np.cos(theta_p)]])
    Mresult = np.matmul(np.matmul(M1,M2),M3)
    eps_11 = Mresult[0,0]
    eps_22 = Mresult[1,1]

    # Strain rate error
    if dx_errs is not None and dy_errs is not None:
        eps_err = (eps_xx_err+eps_yy_err)/2. + np.nanmean([eps_xy_err,eps_yx_err])
    else:
        eps_err = None

    return theta_p, eps_11, eps_22, eps_err
