#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 21 15:38:47 2015

@author: tsz-xhu
"""

from __future__ import division
import numpy as np
import python.model_subproblem_op as model_subproblem_op

class house(object):
    """
    Overview of all methods and usage of this class
    """
    
    def __init__(self, params, eco ,devs, houses):
        """
        """
        
        self.res_x = []
        self.res_y = []
        self.res_energy = []
        self.res_power = []
        self.res_heat = []
        self.res_soc = []
        self.res_p_imp = []
        self.res_p_ch = []
        self.res_p_dch = []
        self.res_p_use = []
        self.res_p_sell = []
        self.res_area = 0
        self.res_cap = 0
        self.res_volume = 0
        self.res_temperature = []
        self.obj = []
        self.res_c_inv = 0
        self.res_c_om = 0
        self.res_c_dem = 0
        self.res_c_met = 0
        self.res_rev = {}
        self.res_chp_sub = []
        self.res_soc_nom = []
        self.res_power_nom = []
        self.res_heat_nom = []
        self.res_soc_init = []
        self.res_p_use = {}
        self.res_p_sell ={}
        self.devs = devs
        self.params = params
        self.eco = eco
        self.houses = houses
        
    def compute_proposal(self, houses, marginals, eco, devs, clustered, params, housedev, iteration=0, t_init=0):
        """
        This function computes a new proposal (P and k).
        Internally, the results of the subproblem are stored.
        If this is the first time in the current optimization period that new
            proposals have to generated, _iteration_ has to be set to 0
        """
        self.iteration = iteration
        if iteration == 0:
            self.opti_res = {}
            self.houses = houses
            opti_res = {}
           
            opti_res = model_subproblem_op.compute(houses, marginals, eco, devs, clustered, params, housedev)
        
        (res_y, res_energy, res_power, res_heat, res_soc, res_p_imp,
         res_p_ch, res_p_dch, res_p_use, res_p_sell, res_cap,
         res_volume, res_temperature, obj, res_c_dem,
         res_c_met, res_soc_nom, res_power_nom,
         res_heat_nom, res_soc_init, devs, costs, proposals, cost, objVal, runtime, mipgap) = opti_res
        
        self.opti_res = opti_res  

        #k = 0  #Heat pump proposals do not generate costs!
        
        return (opti_res)