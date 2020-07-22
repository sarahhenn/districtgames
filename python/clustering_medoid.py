#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 20 10:30:27 2016

@author: tsz-xhu
"""

from __future__ import division
import numpy as np
import math
import python.k_medoids as k_medoids

def _distances(values, norm=2):
    """
    """
    d = np.zeros((values.shape[1], values.shape[1]))

    dist = (lambda day1, day2, r: 
            math.pow(np.sum(np.power(np.abs(day1 - day2), r)), 1/r))

    # Remember: This matrix is symmetrical!
    for i in range(values.shape[1]): # loop over first days
        for j in range(i+1, values.shape[1]): # loop second days
            d[i, j] = dist(values[:,i], values[:,j], norm)
    
    # Fill the remaining entries
    d = d + d.T
    
    return d


def cluster(inputs, number_clusters=12, norm=2, time_limit=300, mip_gap=0.0):
    """
    """
    # Determine time steps per day
    len_day = int(inputs.shape[1] / 365)
    
    # Manipulate inputs
    # Initialize arrays
    inputsTransformed = []
    inputsScaled = []
    inputsScaledTransformed = []
    
    # Fill and reshape
    # Scaling to values between 0 and 1, thus all inputs shall have the same
    # weight and will be clustered equally in terms of quality 
    for i in range(inputs.shape[0]):
        vals = inputs[i,:]
        temp = (vals - np.min(vals)) / (np.max(vals) - np.min(vals))
        inputsScaled.append(temp)
        inputsScaledTransformed.append(temp.reshape((len_day, 365), order="F"))
        inputsTransformed.append(vals.reshape((len_day, 365), order="F"))

    # Put the scaled and reshaped inputs together
    L = np.concatenate(tuple(inputsScaledTransformed))

    # Compute distances
    d = _distances(L, norm)

    # Execute optimization model
    (y, z, obj) = k_medoids.k_medoids(d, number_clusters, time_limit, mip_gap)
    
    # Retain typical days
    nc = np.zeros_like(y)
    typicalDays = []

    # nc contains how many days are there in each cluster
    nc = []
    for i in range(len(y)):
        temp = np.sum(z[i,:])
        if temp > 0:
            nc.append(temp)
            typicalDays.append([ins[:,i] for ins in inputsTransformed])

    typicalDays = np.array(typicalDays)
    nc = np.array(nc, dtype="int")
    nc_cumsum = np.cumsum(nc) * len_day

    # Construct (yearly) load curves
    # ub = upper bound, lb = lower bound
    clustered = np.zeros_like(inputs)
    for i in range(len(nc)):
        if i == 0:
            lb = 0
        else:
            lb = nc_cumsum[i-1]
        ub = nc_cumsum[i]
        
        for j in range(len(inputsTransformed)):
            clustered[j, lb:ub] = np.tile(typicalDays[i][j], nc[i])

    # Scaling to preserve original demands
    sums_inputs = [np.sum(inputs[j,:]) for j in range(inputs.shape[0])]
    scaled = np.array([nc[day] * typicalDays[day,:,:] 
                       for day in range(number_clusters)])
    sums_scaled = [np.sum(scaled[:,j,:]) for j in range(inputs.shape[0])]
    scaling_factors = [sums_inputs[j] / sums_scaled[j] 
                       for j in range(inputs.shape[0])]
    scaled_typ_days = [scaling_factors[j] * typicalDays[:,j,:]
                       for j in range(inputs.shape[0])]
    
    return (scaled_typ_days, nc, z)