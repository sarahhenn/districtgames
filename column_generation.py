# -*- coding: utf-8 -*-
"""
Run file for hierarchical district control

@author: Sarah Henn, based on Thomas Schütz, Xiaolin Hu
"""

from __future__ import division

import object_subproblem
import object_masterproblem
import python.parse_inputs as parse_inputs
import python.clustering_medoid as clustering
import numpy as np
import datetime
"""
Inputs: heat demand (includes space heating and dhw)
        electricity demand
        solar irradiation for PV and STC
        design heat load
        weights for typical days
        ambient temperature
"""
time = {}
datetime.datetime.now()
time["begin"] = datetime.datetime.now()
print ("This program begin at " + str(datetime.datetime.now()) + ".")

# Create lists with buildings
    
design_heat_loads = {}

with open("design_heat_load.txt") as f:
    for line in f:
        splitLine = line.split()
design_heat_loads = {splitLine[2*n] : int(splitLine[(2*n+1)]) / 1000 for n in range(int(len(splitLine)/2))} 

houses = []
for n in range(int(len(splitLine)/2)):
    houses.append(splitLine[2*n])

house = []

number_houses = len(houses)

raw_inputs = {}
raw_inputs["solar_irrad"] = np.loadtxt("raw_inputs/solar_rad_35deg.csv") / 1000
raw_inputs["temperature"] = np.loadtxt("raw_inputs/temperature.csv")
dhw = {}
sh = {}

for n in range(number_houses):
    dhw[n] = np.loadtxt("raw_inputs/dhw_"+houses[n]+".csv")
    sh[n] = np.loadtxt("raw_inputs/heat_"+houses[n]+".csv")
    raw_inputs[n] = {"electricity": np.loadtxt("raw_inputs/elec_"+houses[n]+".csv") / 1000,
                     "heat": (dhw[n] + sh[n]) / 1000}

# Clustering

inputs_clustering = []
for n in range(number_houses):
    inputs_clustering.append(raw_inputs[n]["electricity"])
    inputs_clustering.append(raw_inputs[n]["heat"])
inputs_clustering.append(raw_inputs["solar_irrad"])
inputs_clustering.append(raw_inputs["temperature"])

number_clusters = 12 # Determine the number of type days

(inputs, nc, z) = clustering.cluster(np.array(inputs_clustering),
                                     number_clusters=number_clusters,
                                     norm=2,
                                     mip_gap=0.01)

len_day = np.shape(inputs[0])[1] # Determine time steps per day

clustered = {}
for n in range(number_houses):
    clustered[n] = {"electricity":      inputs[2*n],
                    "heat":             inputs[2*n+1],
                    "design_heat_load": design_heat_loads[houses[n]],
                    "temperature":      inputs[-1],
                    "solar_irrad":      inputs[-2],
                    "weights":          nc}
#clustered["z"]           = z

# Read devices, economic date and other parameters
devs = {}

for n in range(number_houses):
    devs[n] = parse_inputs.read_devices(timesteps=len_day, days=number_clusters,
                         temperature_ambient=clustered[n]["temperature"],
                         temperature_flow=35, 
                         temperature_design=-12, 
                         solar_irradiation=clustered[n]["solar_irrad"])
    (eco, par, devs[n]) = parse_inputs.read_economics(devs[n])

par = parse_inputs.compute_parameters(par, number_clusters, len_day)

days = range(number_clusters)
dt = par["dt"]
times = range(24)

# Set economic parameters
#eco["sub_chp"] = 2 * eco["sub_chp"] # Reduce or add the subsidy for chp
#
#dev = ["bat","pv","eh","stc","hp","chp","boiler","tes"]
#
#for n in range(number_houses):
#    for i in dev:
#        devs[n][i]["c_inv_fix"] = devs[n][i]["c_inv_fix"] * 1
#        devs[n][i]["c_inv_var"] = devs[n][i]["c_inv_var"] * 1

# Subproblem Object
for n in range(len(houses)):
    house.append(object_subproblem.house(par, eco, devs, houses))

# P_demand for all the buildings in each time step
P_demand = {}
for d in days:
    for t in times:
        P_demand[d,t] = sum(clustered[n]["electricity"][d,t] for n in range(number_houses))     

# Create masterproblem
mp = object_masterproblem.Master(len(houses), par, houses, eco, P_demand, clustered)

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

iteration = 2 # Determine the number of iterations

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
        opti_res[it_counter][n] = house[n].compute_proposal(houses, marginals, eco, devs[n], clustered[n], par)      
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
filename = "results//results_column_generation_"+str(number_houses)+"_buildings_"+str(number_clusters)+"_typtage_"+str(iteration)+"_iteration"+".pkl"
with open(filename, "wb") as f_in:
    pickle.dump(opti_res, f_in, pickle.HIGHEST_PROTOCOL)
    pickle.dump(eco, f_in, pickle.HIGHEST_PROTOCOL)
    pickle.dump(devs, f_in, pickle.HIGHEST_PROTOCOL)
    pickle.dump(clustered, f_in, pickle.HIGHEST_PROTOCOL)
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

# Retrieve optimal schedules from each house:
#x_hp = np.zeros((timesteps, len(hp_nom)))
#P_hp = np.zeros((timesteps, len(hp_nom)))
#for i in xrange(len(hp_nom)):
#    (temp_x, temp_y, temp_T, temp_P) = hp[i].get_optimal_schedule(lambda_hp[:,i])
#    x_hp[:,i] = temp_x
#    P_hp[:,i] = temp_P

#x_chp = np.zeros((timesteps, len(chp_nom)))    
#P_chp = np.zeros((timesteps, len(chp_nom)))
#for j in xrange(len(chp_nom)):
#    (temp_x, temp_y, temp_T, temp_P, temp_Q) = chp[j].get_optimal_schedule(lambda_chp[:,j])
#    x_chp[:,j] = temp_x
#    P_chp[:,j] = temp_P