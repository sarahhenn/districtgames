#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Subproblem optimization model for hierarchical district control 

@author: Sarah Henn, based on Thomas Schütz, Xiaolin Hu

"""

from __future__ import division

import gurobipy as gp
import numpy as np
import datetime

def compute(pass_house, marginals, eco, devs, demands, params):
    """
    Parameters
    ----------
    eco : dict
        - b :
        - crf :
        - prChange : 
        - q : 
        - rate : 
        - sub_CHP : 
        - t_calc : 
        - tax : 
        - gas, Boiler : 
        - gas, CHP : 
        - gas, c_meter : 
        - pr, el : 
        - sell, CHP : 
        - sell, PV : 
    devs : dict
        - bat :
        - boiler :
        - chp :
        - eh :
        - hp :
        - pv :
        - stc :
        - tes :
    demands : dict
        - electricity : 
        - heat : 
        - solar_irrad : 
        - temperature : 
        - weights : 
    params : dict
        - dt : time step length (s)
        - maximum roof area (m2)
        - mip_gap : Solver setting (-)
        - time_limit : Solver setting (s)
        - days : 
        - time_steps : 
    """
    
    # Extract parameters
    dt = params["dt"]
    sigma = marginals["sigma"]
    
    # Define subsets
    heater = ("boiler", "chp", "eh", "hp")
    storage = ("bat", "tes")
    solar = ("pv", "stc")
    device = heater+storage+solar
    
    time_steps = range(params["time_steps"])
    days = range(params["days"])
    
    try:
        model = gp.Model("Design computation")
        
        # Define variables
        # Costs and Revenues
        c_dem = {dev: model.addVar(vtype="C", name="c_dem_"+dev)
                  for dev in ("boiler", "chp", "grid")}
        c_meter = model.addVar(vtype="C", name="c_meter")
        
        revenue = {dev: model.addVar(vtype="C", name="revenue_"+dev)
                  for dev in ("chp", "pv")}
        chp_subsidy = model.addVar(vtype="C", name="chp_subsidy")
        
        # SOC, power, heat and energy
        soc = {}
        power = {}
        heat = {}
        energy = {}
        soc_nom = {}
        power_nom = {}
        heat_nom = {}
        for d in days: # All days
            for t in time_steps: # All time steps of all days
                timetag = "_"+str(d)+"_"+str(t)
                for dev in storage: # All storage devices
                    soc[dev,d,t] = model.addVar(vtype="C",
                                                name="SOC_"+dev+"_"+timetag)
                
                for dev in (heater+solar):
                    power[dev,d,t] = model.addVar(vtype="C",
                                                  name="P_"+dev+"_"+timetag)
                    heat[dev,d,t] = model.addVar(vtype="C",
                                                 name="Q_"+dev+"_"+timetag)

                for dev in heater:
                    energy[dev,d,t] = model.addVar(vtype="C",
                                                   name="E_"+dev+"_"+timetag)
        
                for dev in heater:
                    heat_nom[dev,d,t] = model.addVar(vtype="C",
                                              name="Q_nom_"+dev+"_"+timetag)
                dev = "hp"
                power_nom[dev,d,t] = model.addVar(vtype="C",
                                           name="P_nom_"+dev+"_"+timetag)

        # Storage initial SOC's
        soc_init = {}
        for dev in storage:
            for d in days:
                tag = dev + "_" + "_" + str(d)
                soc_init[dev,d] = model.addVar(vtype="C", name="SOC_init_"+tag)


        # Storage charging and discharging
        p_ch = {}
        p_dch = {}
        for d in days:
            for t in time_steps:
                timetag = "_" + str(d) + "_" + str(t)

                p_ch[d,t] = model.addVar(vtype="C", name="P_ch"+timetag)
                p_dch[d,t] = model.addVar(vtype="C", name="P_dch"+timetag)
                    
        
        # Electricity imports, sold and self-used electricity
        p_imp = {}
        p_use = {}
        p_sell = {}
        for d in days:
            for t in time_steps:
                timetag = "_"+str(d)+"_"+str(t)
                
                p_imp[d,t] = model.addVar(vtype="C", name="P_imp_"+timetag)
                for dev in ("chp", "pv"):
                    p_use[dev,d,t] = model.addVar(vtype="C", 
                                                  name="P_use_"+dev+timetag)
                    p_sell[dev,d,t] = model.addVar(vtype="C", 
                                                   name="P_sell_"+dev+timetag)

        # Existing parameters
        # TODO: import from extern
        
        # Variables for determining the variable costs (based on capacity)
        capacity = {}
        capacity["bat"] = 10.0
        capacity["boiler"] = 0.0
        capacity["chp"] = 0.0
        capacity["eh"] =  5.0  
        capacity["hp"] =  5.0  
        capacity["pv"] =  30.0    # in m2
        capacity["inv"] =  20.0    
        capacity["stc"] = 0.0     # in m2
        capacity["tes"] = 300.0   # in m3
        
        x = {}
        for dev in device:
            x[dev] = 1
        
        # mapping PV and STC areas
        area = {}
        for dev in solar:
            area[dev] = capacity[dev]
            
        # maping storage sizes
        soc_nom = {}
        dev = "tes"
        volume = capacity[dev]
        soc_nom[dev] = params["rho_w"] * params["c_w"] * devs[dev]["dT_max"] / 3600000 * volume
        dev = "bat"
        soc_nom[dev] = capacity[dev]      
        
        # Activation decision variables
        # TODO: check if needed
        y = {}  # Acitivation (heaters)       
        for d in days:
            for t in time_steps:
                timetag = "_"+str(d)+"_"+str(t)
                for dev in heater: # All heating devices
                    y[dev,d,t] = model.addVar(vtype="B", lb=0.0, ub=1.0,
                                              name="y_"+dev+"_"+timetag)
        
        # Update
        model.update()

        # Define constraints
        # Objective
        model.setObjective(sum(c_dem[key] for key in c_dem.keys())
                         + c_meter
                         - chp_subsidy
                         - sum(revenue[key] for key in revenue.keys())
                         - sigma, 
                         gp.GRB.MINIMIZE)


        # Determine niminal heat at every timestep
        for d in days:
            for t in time_steps:
                for dev in heater:
                    model.addConstr(capacity[dev] >= heat_nom[dev,d,t],
                                    name="Capacity_"+dev+"_"+str(d)+"_"+str(t))
        
        # Economic constraints

        # Demand related costs (gas)
        for dev in ("boiler", "chp"):
            model.addConstr(c_dem[dev] == eco["b"]["gas"] * eco["crf"] *
                eco["gas",dev] * dt * 
                sum(demands["weights"][d] * 
                    sum(energy[dev,d,t] for t in time_steps)
                    for d in days),
                name="Demand_costs_"+dev)
                
        # Demand related costs (electricity)    # eco["b"]["el"] * eco["crf"]
        dev = "grid"
        model.addConstr(c_dem[dev] == sum(#demands["weights"][d] * 
#                sum((demands["electricity"][d,t] + power["hp",d,t] +power["eh",d,t] -power["chp",d,t] -power["pv",d,t]) * marginals["pi"][d,t] for t in time_steps)
                 sum((p_imp[d,t]) * marginals["pi"][d,t] for t in time_steps)
                 for d in days),
            name="Demand_costs_"+dev)
        
        # Revenues for selling electricity to the grid / neighborhood
        for dev in ("chp", "pv"):
            model.addConstr(revenue[dev] == sum(#demands["weights"][d] * 
                    sum(p_sell[dev,d,t] * marginals["pi"][d,t] for t in time_steps) 
                    for d in days),
                name="Feed_in_rev_"+dev)
                
#        # Revenues for selling electricity to the grid
#        for dev in ("chp", "pv"):
#            model.addConstr(revenue[dev] == eco["b"]["eex"] * eco["crf"] * dt *
#                sum(demands["weights"][d] *
#                    sum(p_sell[dev,d,t] * eco["sell",dev] 
#                    for t in time_steps)
#                for d in days),
#                name="Feed_in_rev_"+dev)

        # Metering costs
        for dev in ("boiler", "chp"):
            model.addConstr(c_meter >= eco["b"]["infl"] * eco["crf"] 
                            * eco["gas","c_meter"] * x[dev],
                            name="Metering_costs_"+dev)

        # CHP subsidies
        model.addConstr(chp_subsidy == eco["b"]["eex"] * eco["crf"] * dt *
                sum(demands["weights"][d] *
                    sum(power["chp",d,t] * eco["sub_chp"] 
                    for t in time_steps)
                for d in days),
                name="Subsidies_chp")       

        # Technical constraints
        
        # Devices can be switched on only if they exist
        for dev in heater:
            model.addConstr(sum( sum(y[dev,d,t] for t in time_steps)
                            for d in days) 
                  <= params["time_steps"] * params["days"] * x[dev],
                  name="Activation_"+dev)

        # Compute nominal power consumption of HP:
        dev = "hp"
        for d in days:
            for t in time_steps:
                model.addConstr(power_nom[dev,d,t] * devs[dev]["cop_a2w35"]
                                == heat_nom[dev,d,t],
                                name="Power_nom_"+dev+"_"+str(d)+"_"+str(t))

        # Devices operation
        # Heat output between mod_lvl*Q_nom and Q_nom (P_nom for heat pumps)
        # Power and Energy directly result from Heat output
        for dev in heater:
            for d in days:
                for t in time_steps:
                    # Abbreviations
                    timetag = "_" + str(d) + "_" + str(t)
                    
                    mod_lvl = devs[dev]["mod_lvl"]
                    eta     = devs[dev]["eta"][d,t]
                    omega   = devs[dev]["omega"][d,t]
                    
                    if dev == "hp":
                        model.addConstr(power[dev,d,t] <= power_nom[dev,d,t],
                            name="Max_pow_operation_"+dev+"_"+timetag)
                        model.addConstr(power[dev,d,t] >= 
                            power_nom[dev,d,t] * mod_lvl,
                            name="Min_pow_operation_"+dev+"_"+timetag)
                    else:
                        model.addConstr(heat[dev,d,t] <= heat_nom[dev,d,t],
                            name="Max_heat_operation_"+dev+"_"+timetag)
                        model.addConstr(heat[dev,d,t] >= 
                            heat_nom[dev,d,t] * mod_lvl,
                            name="Min_heat_operation_"+dev+"_"+timetag)
                    
                    model.addConstr(power[dev,d,t] == 1/eta * heat[dev,d,t],
                          name="Power_equation_"+dev+"_"+timetag)
                            
                    model.addConstr(energy[dev,d,t] == 
                          1/omega * (heat[dev,d,t]+power[dev,d,t]),
                          name="Energy_equation_"+dev+"_"+timetag)

        # Solar components
        for dev in solar:
            for d in days:
                for t in time_steps:
                    timetag = "_" + str(d) + "_" + str(t)
                    eta_th = devs[dev]["eta_th"][d][t]
                    eta_el = devs["inv"]["eta"] * devs[dev]["eta_el"][d][t]
                    solar_irrad = demands["solar_irrad"][d][t]

                    model.addConstr(heat[dev,d,t] == 
                          eta_th * area[dev] * solar_irrad,
                          name="Solar_thermal_"+dev+timetag)
                    model.addConstr(power[dev,d,t] == 
                          eta_el * area[dev] * solar_irrad,
                          name="Solar_electrical_"+dev+timetag)
        
        # Storage equations
        # TES soc[dev,d,t] soc_init[dev,d] soc_nom[dev]
        dev = "tes"
        k_loss = devs[dev]["k_loss"]
        eta_ch = devs[dev]["eta_ch"]
        eta_dch = devs[dev]["eta_dch"]
        for d in days:
            for t in time_steps:
                if t == 0:
                    if np.max(demands["weights"]) == 1:
                        if d == 0:
                           soc_prev = soc_init[dev,d]
                        else:
                           soc_prev = soc[dev,d-1,params["time_steps"]-1]
                    else:
                        soc_prev = soc_init[dev,d]
                else:
                    soc_prev = soc[dev,d,t-1]
                
                timetag = "_" + str(d) + "_" + str(t)
                
                charge = eta_ch * sum(heat[dv,d,t] for dv in (heater+solar))
                discharge = 1 / eta_dch * demands["heat"][d,t]
                
                model.addConstr(soc[dev,d,t] == (1-k_loss) * soc_prev + 
                                dt * (charge - discharge),
                                name="Storage_bal_"+dev+timetag)
        
        # Nominal storage content (SOC)
        for dev in storage:
            for d in days:
                # Inits
                model.addConstr(soc_nom[dev] >= soc_init[dev,d], 
                                name="SOC_nom_inits_"+dev+"_"+str(d))
                for t in time_steps:
                    # Regular storage loads
                    model.addConstr(soc_nom[dev] >= soc[dev,d,t],
                                    name="SOC_nom_"+dev+"_"+str(d)+"_"+str(t))

        dev = "bat"
        k_loss = devs[dev]["k_loss"]
        for d in days:
            for t in time_steps:
                if t == 0:
                    if np.max(demands["weights"]) == 1:
                        if d == 0:
                           soc_prev = soc_init[dev,d]
                        else:
                           soc_prev = soc[dev,d-1,params["time_steps"]-1]
                    else:
                        soc_prev = soc_init[dev,d]
                else:
                    soc_prev = soc[dev,d,t-1]

                timetag = "_" + str(d) + "_" + str(t)
    
                model.addConstr(soc[dev,d,t] == (1-k_loss) * soc_prev + 
                    dt * (devs[dev]["eta"] * p_ch[d,t] - 
                          1 / devs[dev]["eta"] * p_dch[d,t]),
                    name="Storage_balance_"+dev+timetag)
    
                model.addConstr(p_ch[d,t] <= x[dev] * devs[dev]["P_ch_fix"] + 
                                capacity[dev] * devs[dev]["P_ch_var"],
                                name="P_ch_max"+timetag)
    
                model.addConstr(p_dch[d,t] <= x[dev] * devs[dev]["P_dch_fix"] + 
                                capacity[dev] * devs[dev]["P_dch_var"],
                                name="P_dch_max"+timetag)
                
        # SOC repetitions
        for dev in storage:
            for d in range(params["days"]):
                if np.max(demands["weights"]) > 1:
                    model.addConstr(soc[dev,d,params["time_steps"]-1] ==
                                    soc_init[dev,d],
                                    name="repetitions_"+dev+"_"+str(d))

        # Electricity balance (house)
        for d in days:
            for t in time_steps:
                model.addConstr(demands["electricity"][d,t]
                    + p_ch[d,t] - p_dch[d,t]
                    + power["hp",d,t] + power["eh",d,t]
                    - p_use["chp",d,t] - p_use["pv",d,t]
                    == p_imp[d,t],
                    name="Electricity_balance_"+str(d)+"_"+str(t))
        
        # Split CHP and PV generation into self-consumed and sold powers
        for dev in ("chp", "pv"):
            for d in days:
                for t in time_steps:
                    model.addConstr(p_sell[dev,d,t] + p_use[dev,d,t] 
                            == power[dev,d,t],
                            name="power=sell+use_"+dev+"_"+str(d)+"_"+str(t))
                            
        # Heat pump's operation depends on storage temperature
        for d in days:
            for t in time_steps:
                # Abbreviations
                dT_relative = devs["hp"]["dT_max"] / devs["tes"]["dT_max"]
                # Residual storage content
                resSC = (devs["tes"]["volume_max"] * devs["tes"]["dT_max"]
                         * params["rho_w"] * params["c_w"] * (1 - dT_relative)
                         / 3600000)
                
                model.addConstr(soc["tes",d,t] <= soc_nom["tes"] * dT_relative 
                      + (1 - y["hp",d,t]) * resSC,
                      name="Heat_pump_act_"+str(d)+"_"+str(t))
        
        
        # Set solver parameters
        model.Params.TimeLimit = params["time_limit"]
        model.Params.MIPGap = params["mip_gap"]
        model.Params.MIPFocus = 3
        
        # Execute calculation
        model.optimize()
#        model.computeIIS()
#        model.write("model.ilp")
        
        # Retrieve results
        res_y = {}
        res_energy = {}
        res_power = {}
        res_heat = {}
        res_soc = {}
        
        for dev in heater:
            res_y[dev] = {(d,t): y[dev,d,t].X
                                 for d in days
                                 for t in time_steps}
            res_energy[dev] = {(d,t) : energy[dev,d,t].X  
                                       for d in days
                                       for t in time_steps}
        for dev in (heater + solar):
            res_power[dev] = {(d,t) : power[dev,d,t].X  
                                      for d in days
                                      for t in time_steps}
        
            res_heat[dev] = {(d,t) : heat[dev,d,t].X  
                                     for d in days
                                     for t in time_steps}
        
        for dev in storage:
            res_soc[dev] = {(d,t): soc[dev,d,t].X 
                                   for d in days
                                   for t in time_steps}
    
        res_p_imp = {(d,t) : p_imp[d,t].X for d in days for t in time_steps}
        res_p_ch = {(d,t) : p_ch[d,t].X for d in days for t in time_steps}
        res_p_dch = {(d,t) : p_dch[d,t].X for d in days for t in time_steps}
        
        res_cap = {dev : capacity[dev] for dev in capacity.keys()}
        
        res_volume = max(0.001, volume)
        res_temperature = {(d,t): res_soc["tes"][d,t]/(res_volume*1000*4180) for d in days for t in time_steps}
        
        obj = model.ObjVal
        
        
        res_c_dem = {dev: c_dem[dev].X for dev in c_dem.keys()}
        res_c_met = c_meter.X
        
        # res_rev = {dev: revenue[dev].X for dev in revenue.keys()}
#        res_chp_sub = chp_subsidy.X
        
        res_soc_nom = {dev: soc_nom[dev] for dev in storage}
        res_power_nom = {}
        res_heat_nom = {}
        for dev in heater:
            res_heat_nom[dev] = {(d,t): heat_nom[dev,d,t].X for d in days for t in time_steps}
        dev = "hp"
        res_power_nom[dev] = {(d,t): power_nom[dev,d,t].X for d in days for t in time_steps}
        
        res_soc_init = {}
        for dev in storage:
            res_soc_init[dev] = {d: soc_init[dev,d].X for d in days}
        
        res_p_use = {}
        res_p_sell = {}
        for dev in ("chp", "pv"):
            res_p_use[dev] = {(d,t): p_use[dev,d,t].X for d in days for t in time_steps}
            res_p_sell[dev] = {(d,t): p_sell[dev,d,t].X for d in days for t in time_steps}
        
        print
        for dev in devs:
            print ("Capacity " + dev + ": " + str(res_cap[dev]))
                
        print ()
        print ("Obj: " + str(model.ObjVal))
        
        # Compute "costs" of the proposal:
        cost = {}
        cost["inv"] = []
        cost["om"] = []
        cost["dem"] = []
        cost["grid"] = []
        cost["dem"] = res_c_dem["boiler"]+res_c_dem["chp"]
        cost["grid"] = res_c_dem["grid"]
#        costs = (res_c_met - res_chp_sub + cost["inv"] + cost["om"] +cost["dem"])
        costs = (cost["dem"])
            
        # Compute "proposals" of the devices:
        proposals = {}
        proposals["house"] = {}
        for dev in ("chp", "pv"):
            proposals[dev] = {(d,t): power[dev,d,t].X  for d in days for t in time_steps}
        proposals["hp"] = res_power["hp"]
        proposals["eh"] = res_power["eh"]
        for d in days:
            for t in time_steps:
                proposals["house"][d,t] = res_p_imp[d,t] - res_p_sell["chp"][d,t] - res_p_sell["pv"][d,t]
    
        #res_p_imp
        #res_p_imp = {}
        
        objVal = obj
        runtime = model.getAttr("Runtime")
#        mipgap = model.getAttr("MIPgap") 
        mipgap = 0
        datetime.datetime.now()   
#        model.computeIIS()
#        model.write("model.ilp")
#        print('\nConstraints:')        
#        for c in model.getConstrs():
#            if c.IISConstr:
#                print('%s' % c.constrName)
#        print('\nBounds:')
#        for v in model.getVars():
#            if v.IISLB > 0 :
#                print('Lower bound: %s' % v.VarName)
#            elif v.IISUB > 0:
#                print('Upper bound: %s' % v.VarName) 
        
        # Return results
        return (res_y, res_energy, res_power, res_heat, res_soc, res_p_imp,
                res_p_ch, res_p_dch, res_p_use, res_p_sell, res_cap,
                res_volume, res_temperature, obj, res_c_dem,
                res_c_met, res_soc_nom, res_power_nom,
                res_heat_nom, res_soc_init, devs, costs, proposals, cost, objVal, runtime, mipgap)
        
    except gp.GurobiError:
        print("Error")
    