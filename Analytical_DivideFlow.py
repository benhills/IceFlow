#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 4 2018

@author: benhills
"""

def flow_dansgaard_johnson(ws,Eps,Eps_c):
    """
    Vertical velocity from Kingslake et al., 2014
    """
    w = ws*Eps**2. / (Eps_c*(2.-Eps_c))
    w[Eps >= Eps_c] = ws * (2.*Eps[Eps>=Eps_c]-Eps_c)/(2.-Eps_c)
    return w

def flow_lliboutry(parameters,Eps):
    """
    Vertical velocity from Kingslake et al., 2014
    """
    ws = parameters[0]
    p = parameters[1]
    w = ws * (1. - (p+2.)/(p+1.)*Eps + 1./(p+1.)*Eps**(p+2.))
    return w
