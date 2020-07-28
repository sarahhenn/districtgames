# -*- coding: utf-8 -*-
"""
Run file for hierarchical district control

@author: Sarah Henn, based on Thomas Schütz, Xiaolin Hu
"""

from __future__ import division

import python.object_subproblem as object_subproblem
import python.object_masterproblem as object_masterproblem
import python.parse_inputs as parse_inputs
import numpy as np
import datetime
"""
Inputs: heat demand (includes space heating and dhw)
        electricity demand
        solar irradiation for PV and STC
        ambient temperature
"""

#%% INPUT PARAMETERS

# adjust the following values 
obs_period = 10          # Observation period in days
obs_houses = 113        # Number of buildings to be observed (max 113)
iteration = 10          # Determine the number of iterations


#%% 

time = {}
datetime.datetime.now()
time["begin"] = datetime.datetime.now()
print ("This program begin at " + str(datetime.datetime.now()) + ".")

weather = {}
weather["solar_irrad"] = np.loadtxt("raw_inputs/solar_rad_35deg.csv", max_rows=obs_period*24) / 1000
weather["temperature"] = np.loadtxt("raw_inputs/temperature.csv", max_rows=obs_period*24)
    
# load list of buildings, load demands per building
    
nodes,houses = parse_inputs.read_demands("nodes.txt", obs_period, obs_houses)

fullhorizon = len(weather["temperature"])

number_houses = len(houses)

# load devs per building
housedevs = parse_inputs.read_housedevs()

# Read devices, economic date and other parameters
devs = {}

devs = parse_inputs.read_devices(timesteps=fullhorizon,
                         temperature_ambient=weather["temperature"],
                         temperature_flow=35, 
                         temperature_design=-12, 
                         solar_irradiation=weather["solar_irrad"])

(eco, par, devs) = parse_inputs.read_economics(devs)

par = parse_inputs.compute_parameters(par, fullhorizon)

days = 1
dt = par["dt"]
#times = range(fullhorizon)
times = range(120)

# Subproblem Object
house = []

for n in range(len(houses)):
    house.append(object_subproblem.house(par, eco, devs, houses, nodes, housedevs, weather))

# P_demand for all the buildings in each time step
P_demand = {}
for t in times:
    P_demand[t] = sum(nodes[n]["elec"][t] for n in range(number_houses))     

# Create masterproblem
mp = object_masterproblem.Master(len(houses), par, houses, eco, P_demand, nodes)

# Store results of masterproblem
res_obj = []
res_marginals = []
res_costs = {}
res_proposals = {}
    
# Initialize masterproblem
(r_obj, r) = mp.update_proposals({},{}, houses)
res_obj.append(r_obj)
res_marginals.append(r["pi"])

it_counter = 0 # Iteration counter

opti_res = {}

while it_counter < iteration:
    print()
    print ("Begin iteration " + str(it_counter))
    print()
    datetime.datetime.now() # get the time 
    time_begin = datetime.datetime.now()
    print ("*********************************************")
    print ("The " + str(it_counter) + " iteration begins at " + str(datetime.datetime.now())+".")
    print ("*********************************************")
    costs = []
    proposals = {}
    proposals["chp"] = []
    proposals["hp"] = []
    proposals["pv"] = []
    proposals["eh"] = []
    proposals["house"] = []    
    res_costs[it_counter] = []
    opti_res[it_counter] = {}
    
    for n in range(number_houses):
        marginals = {}        
        marginals["sigma"] = r["sigma"][n]
        marginals["pi"] = r["pi"]
        opti_res[it_counter][n] = house[n].compute_proposal(houses, marginals, eco, devs, nodes[n], par, housedevs.iloc[n], weather)      
        costs.append(opti_res[it_counter][n][26])
        res_costs[it_counter].append(opti_res[it_counter][n][21])
        proposals["chp"].append(opti_res[it_counter][n][22]["chp"])
        proposals["hp"].append(opti_res[it_counter][n][22]["hp"])
        proposals["pv"].append(opti_res[it_counter][n][22]["pv"])
        proposals["eh"].append(opti_res[it_counter][n][22]["eh"])
        proposals["house"].append(opti_res[it_counter][n][22]["house"])
    
    res_proposals[it_counter] = proposals   
         
    (r_obj, r) = mp.update_proposals(costs, proposals, houses)

    res_obj.append(r_obj)
    res_marginals.append(r["pi"])
    print()
    print ("End iteration " + str(it_counter))
    print()
    it_counter += 1
    
    datetime.datetime.now()
    time_interval = datetime.datetime.now() - time_begin
    time["time_interval_"+str(it_counter)] = time_interval.seconds


# Solve masterproblem with binary restrictions
datetime.datetime.now()
time_begin_finalize = datetime.datetime.now()
time["begin_finalize"] = time_begin_finalize
print ("The final iteration begins at " + str(datetime.datetime.now()) +".")

#(obj, lambda_house) = mp.finalize(max_time=200)
(r_obj, lambda_house) = mp.finalize(max_time=200)

datetime.datetime.now()
time_interval_finalize = datetime.datetime.now() - time_begin_finalize
time["time_interval_finalize"] = time_interval_finalize.seconds
time["end"] = datetime.datetime.now()
print ("The program ends at " + str(datetime.datetime.now()) + ".")

# Store the results into a pkl-date
import pickle
filename = "results//results_column_generation_"+str(number_houses)+"_buildings_"+"_typtage_"+str(iteration)+"_iteration"+".pkl"
with open(filename, "wb") as f_in:
    pickle.dump(opti_res, f_in, pickle.HIGHEST_PROTOCOL)
    pickle.dump(eco, f_in, pickle.HIGHEST_PROTOCOL)
    pickle.dump(devs, f_in, pickle.HIGHEST_PROTOCOL)
    pickle.dump(nodes, f_in, pickle.HIGHEST_PROTOCOL)
    pickle.dump(par, f_in, pickle.HIGHEST_PROTOCOL)
    pickle.dump(res_marginals, f_in, pickle.HIGHEST_PROTOCOL)
    pickle.dump(res_obj, f_in, pickle.HIGHEST_PROTOCOL)
    pickle.dump(r_obj, f_in, pickle.HIGHEST_PROTOCOL)
    pickle.dump(proposals, f_in, pickle.HIGHEST_PROTOCOL)
    pickle.dump(r, f_in, pickle.HIGHEST_PROTOCOL)
    pickle.dump(time, f_in, pickle.HIGHEST_PROTOCOL)
    pickle.dump(res_costs, f_in, pickle.HIGHEST_PROTOCOL)
    pickle.dump(lambda_house, f_in, pickle.HIGHEST_PROTOCOL)
    pickle.dump(res_proposals, f_in, pickle.HIGHEST_PROTOCOL)
