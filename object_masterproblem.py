# -*- coding: utf-8 -*-
"""

"""


# http://en.wikipedia.org/wiki/Branch_and_price

import numpy as np
import model_masterproblem_test_2

class Master(object):
    """
    This is the Master class that handles all attributes relevant to the masterproblem
    """
    
    def __init__(self, houses_number, params, houses, eco, P_demand, clustered):
        """
        Initialization parameters are:
            number_chp: Number of installed CHP units
            
            number_hp:  Number of installed HP units
        """
        dev = ["bat","pv","eh","stc","hp","chp","boiler","tes"]
        self.count_iteration = 0
        self.houses_number = houses_number
        self.houses = houses
        self.params = params
        self.eco = eco
        self.P_demand = P_demand
        self.costs = []
        self.proposals = {}
        self.bounds = {}
        self.marginals = {}
        self.marginals["sigma"] = {}
        self.weights = clustered[0]["weights"]
        for n in range(len(houses)):
            self.marginals["sigma"][n] = []
        self.marginals["pi"] = []
        self.proposals["chp"] = []
        self.proposals["hp"] = []
        self.proposals["pv"] = []
        self.proposals["boiler"] = []
        self.proposals["eh"] = []
            # self.bounds
        self.marginals["pi"] = []
        self.proposals["house"] = []
        
    def update_timestep(self, P_demand, houses):
        """
        Each new timestep requires a new load profile as well as a new RES profile
        """
        dev = ["bat","pv","eh","stc","hp","chp","boiler","tes"]
        self.count_iteration = 0
        self.P_demand = P_demand
        self.costs = []
        self.proposals = {}
        self.bounds = {}
        self.marginals = {}
        self.marginals["sigma"] = {}
        for n in range(len(houses)):
            self.marginals["sigma"][n] = []
        # for device in ("pv","hp","chp"):
        #    self.proposals[device] = []
            #self.bounds[device] = []
        self.marginals["pi"] = []
        self.proposals["chp"] = []
        self.proposals["hp"] = []
        self.proposals["pv"] = []
        self.proposals["boiler"] = []
        self.proposals["eh"] = []
        self.proposals["house"] = []
        
    def update_proposals(self, costs, proposals, houses):
        """
        Compute new marginals, based on proposed costs and electricity proposals
        
        In case of initializing Branch&Price with the masterproblem, 
            use empty lists as proposals:
        update_proposals({},{})
        
        This function returns the objective value of the masterproblem as 
            well as the marginals of the resource constraint (electricity balance)
        """
        if self.count_iteration > 0:
            self.append_proposals(costs, proposals, houses)
            #if self.count_iteration > 1:
               # self.update_bounds()
            
        (r_obj, r) = model_masterproblem_test_2.optimize(self, False) 
        self.marginals["pi"].append(r["pi"])
        
        for j in range(len(houses)):        
            self.marginals["sigma"][j].append(r["sigma"][j])
        
        self.count_iteration += 1
        
        return (r_obj, r)

    def append_proposals(self, costs, proposals, houses):
        """
        Append the proposals for costs and electricity production/consumption
        """
        self.costs.append(costs)
        self.proposals["chp"].append(proposals["chp"])
        self.proposals["hp"].append(proposals["hp"])
        self.proposals["pv"].append(proposals["pv"])
        self.proposals["eh"].append(proposals["eh"])       
        self.proposals["house"].append(proposals["house"])  
    
    def finalize(self, max_time=10):
        """
        Finalize the Branch&Price process by solving the masterproblem with 
            binaries instead of continuous variables
        """
        (r_obj, r) = model_masterproblem_test_2.optimize(self, True, max_time)
        
        return (r_obj, r)
#    
#    def update_bounds(self):
#        """
#        Find and eliminate all colums that do not improve the masterproblem.
#        Set the bounds in the masterproblem to zero, if the reduced costs are 
#            higher than the corresponding marginals.
#        """
#        dev = ["bat","pv","eh","stc","hp","chp","boiler","tes"]
#        temp = {}
#        costs = {}
#        for device in dev:
#            temp[device] = np.zeros(self.houses_number)
#            costs[device] = np.array(self.houses_number)
#    
#            temp[device][(costs[device][-1,:] < self.marginals["sigma"][device])[0]]=1
#        
#            self.bounds[device].append(temp[device])