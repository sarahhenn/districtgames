#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 02 10:49:10 2015

@author: tsz
"""

import pickle
import numpy as np

houses = ["AB", "LCS68a", "LCS74"]
house = houses[0]
filename = "results_"+house+".pkl"

## Load results
#with open(filename, "rb") as f_out:
#    res = pickle.load(f_out)
#
#(x, y, z, power, heat, energy, soc, p_imp, ch, dch, p_use, p_sell, c_inv, c_om,
# c_dem, c_meter, rev, chp_sub, soc_init, objVal, runtime, mipgap) = res
#
#filename = "results_"+house+"2.pkl"
#
## Load results
#with open(filename, "rb") as f_out:
#    eco = pickle.load(f_out)
#    devs = pickle.load(f_out)
#    clustered = pickle.load(f_out)
#    par = pickle.load(f_out)

# Load results
with open(filename, "rb") as f_out:
    res = pickle.load(f_out)
    eco = pickle.load(f_out)
    devs = pickle.load(f_out)
    clustered = pickle.load(f_out)
    par = pickle.load(f_out)

(res_x, res_y, res_energy, res_power, res_heat, res_soc, res_p_imp,
 res_p_ch, res_p_dch, res_p_use, res_p_sell, res_area, res_cap,
 res_volume, res_temperature, obj, res_c_inv, res_c_om, res_c_dem,
 res_c_met, res_chp_sub, res_soc_nom, res_power_nom,
 res_heat_nom, res_soc_init, devs, costs, proposals, cost, objVal, runtime, mipgap) = res[0]

###################################################

days = range(len(clustered[0]["weights"]))

device = ["bat","pv","eh","stc","hp","chp","boiler","tes"]

# Determine structure
structure = {dev: np.sum(res_x[dev]) for dev in device}
nominals = {dev: np.sum(res_x[dev] * res_cap[dev]) for dev in ("chp","boiler","hp","eh")}
nominals["tes"] = np.sum(res_x["tes"] * res_volume)
nominals["bat"] = np.sum(res_x["bat"] * res_cap["bat"])
nominals["inv"] = np.sum(res_x["inv"] * res_cap["inv"])
nominals["stc"] = np.sum(res_energy["stc"] * res_cap["stc"])
nominals["pv"] = np.sum(res_energy["pv"] * res_cap["pv"])

# Characteristics
# Total heat and electricity demand (building, not energy system)
demand_electrical_building = 0
demand_thermal_building = 0

# Solar fractions

# Total solar generation
solar_thermal_total = 0
solar_electrical_total = 0
for d in days:
    demand_thermal_building += np.sum(clustered["heat"][d,:] * clustered["weights"][d])
    demand_electrical_building += np.sum(clustered["electricity"][d,:] * clustered["weights"][d])
    
    for n in devs["stc"].keys():
        solar_thermal_total += np.sum(res_heat["stc"][n][d,:] * clustered["weights"][d])
    for n in devs["pv"].keys():
        solar_electrical_total += np.sum(res_power["pv"][n][d,:] * clustered["weights"][d])

# Transform to kWh/a
to_kWh = lambda j: j * par["dt"] / (3600 * 1000)
demand_electrical_building = to_kWh(demand_electrical_building)
demand_thermal_building = to_kWh(demand_thermal_building)
solar_thermal_total = to_kWh(solar_thermal_total)
solar_electrical_total = to_kWh(solar_electrical_total)

solar_fraction_electrical_building = solar_electrical_total / demand_electrical_building
solar_fraction_thermal_building = solar_thermal_total / demand_thermal_building

# Operating and full load hours:
full_load_hours = {}
operating_hours = {}
for dev in ("chp", "boiler", "hp", "eh"):
    full_load_hours[dev] = 0
    operating_hours[dev] = 0
    for d in days:
        for n in res_x[dev].keys():
            operating_hours[dev] += np.sum(clustered["weights"][d] * np.sum(res_y[dev][n][d,:]==1))
            full_load_hours[dev] += np.sum(clustered["weights"][d] * np.sum(res_heat[dev][n][d,:]/devs[dev][n]["Q_nom"]))
    
    full_load_hours[dev] *= par["dt"] / 3600
    operating_hours[dev] *= par["dt"] / 3600


# Investments (not annualized)
costs_total_inv = {}
for dev in ("boiler","chp","hp","eh","inv","tes","bat"):
    costs_total_inv[dev] = np.sum(res_x[dev][n] * devs[dev][n]["c_inv"] for n in res_x[dev].keys())
for dev in ("pv", "stc"):
    costs_total_inv[dev] = np.sum(res_energy[dev][n] * devs[dev][n]["c_inv"] for n in res_x[dev].keys())

# Electricity generation/consumption/grid
electricity_generation = {}
electricity_consumption = {}
electricity_consumption["grid"] = 0
electricity_generation["chp", "total"] = 0
electricity_generation["chp", "sell"] = 0
electricity_generation["chp", "use"] = 0
electricity_generation["pv", "total"] = 0
electricity_generation["pv", "sell"] = 0
electricity_generation["pv", "use"] = 0
for d in days:
    electricity_consumption["grid"] += np.sum(clustered["weights"][d] * res_p_imp[d,:])
    for dev in ("chp", "pv"):
        electricity_generation[dev, "total"] += np.sum(clustered["weights"][d] * (res_p_sell[dev][d,:] + res_p_use[dev][d,:]))
        electricity_generation[dev, "sell"] += np.sum(clustered["weights"][d] * res_p_sell[dev][d,:])
        electricity_generation[dev, "use"] += np.sum(clustered["weights"][d] * res_p_use[dev][d,:])

electricity_consumption["grid"] = to_kWh(electricity_consumption["grid"])
for dev in ("chp", "pv"):
    electricity_generation[dev, "total"] = to_kWh(electricity_generation[dev, "total"])
    electricity_generation[dev, "sell"] = to_kWh(electricity_generation[dev, "sell"])
    electricity_generation[dev, "use"] = to_kWh(electricity_generation[dev, "use"])
    
for dev in ("hp", "eh"):
    electricity_consumption[dev] = 0
    for d in days:
        for n in devs[dev].keys():
            electricity_consumption[dev] += np.sum(clustered["weights"][d] * res_power[dev][n][d,:])

    electricity_consumption[dev] = to_kWh(electricity_consumption[dev])

# Heat generation / consumption
heat_total_consumption = np.sum([clustered["weights"][d] * clustered["heat"][d,:] for d in days])
heat_total_consumption = to_kWh(heat_total_consumption)

heat_total_generation = {}
for dev in ("boiler", "chp", "hp", "eh", "stc"):
    heat_total_generation[dev] = 0
    for d in days:
        for n in devs[dev].keys():
            heat_total_generation[dev] += np.sum(clustered["weights"][d] * res_heat[dev][n][d,:])
    
    heat_total_generation[dev] = to_kWh(heat_total_generation[dev])

# Compute average and max. storage usage
soc_average = {}
soc_max = {}
for dev in ("tes", "bat"):
    soc_average[dev] = 0
    soc_max[dev] = 0
    for n in devs[dev].keys():
        soc_average[dev] += np.mean(res_soc[dev][n])
        soc_max[dev] += np.max(res_soc[dev][n])

# Compute storage losses
sto_loss_tes = (sum([heat_total_generation[dev] for dev in heat_total_generation.keys()]) - heat_total_consumption)
sto_loss_bat = (sum(electricity_generation[dev, "use"] for dev in ("chp", "pv"))
                + electricity_consumption["grid"]
                - sum(electricity_consumption[dev] for dev in ("eh", "hp"))
                - demand_electrical_building
                - sum(electricity_generation[dev, "sell"] for dev in ("chp", "pv")))

heat_fractions = {dev: heat_total_generation[dev] / (heat_total_consumption + sto_loss_tes) for dev in heat_total_generation.keys()}
power_building_total = electricity_generation["chp", "use"] + electricity_generation["pv", "use"] + electricity_consumption["grid"]
power_fractions = {}
power_fractions["chp"] = electricity_generation["chp", "use"] / power_building_total
power_fractions["pv"] = electricity_generation["pv", "use"] / power_building_total
power_fractions["grid"] = electricity_consumption["grid"] / power_building_total