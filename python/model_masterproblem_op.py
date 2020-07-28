# -*- coding: utf-8 -*-
"""
Masterproblem optimization model for hierarchical district control 

@author: Sarah Henn, based on Thomas Schütz, Xiaolin Hu

"""
from __future__ import division

import gurobipy as gp
import numpy as np

def optimize(mp, final_iteration=False, max_time=10):

    times = range(mp.params["time_steps"])
#    days = range(mp.params["days"])
#    dev = ["bat","pv","eh","stc","hp","chp","boiler","tes"]
#    days_number = len(days)
    times_number = len(times)
    
    houses_number = len(mp.houses)  
    prop = {}
    number_props = {}
    k = {}
    
    prop["house"] = np.array(mp.proposals["house"])
    number_props["house"] = np.array(np.shape(prop["house"]))[0]

    if number_props["house"] == 0:             
        number_props["house"] = 1
        prop["house"] = np.zeros((1,houses_number,times_number))
        k = np.zeros((1,houses_number))
        
    else:
        k = np.array(mp.costs)
                        
        prop["house"] = np.array(mp.proposals["house"])
        number_props["house"] = (np.shape(prop["house"]))[0]
            
#    P_dem = mp.P_demand  
    dt = mp.params["dt"]
    k_el = mp.eco["pr","el"] # Electricity price
    r_el = mp.eco["sell","chp"] # Feed-in compensation for PV units
      
    # Gurobi optimization model
    try:
        # Create a new model
        model = gp.Model("masterproblem")

        # Create variables with one or more sets:
        l = {} # Weighting variables of all proposals
        l["house"] = {} # Weighting variables of each device
        P_imp = {} # Imported electricity
        P_exp = {} # Exported electricity
   
        if final_iteration:   # Final iteration calculate with binary variables       
            for p in range(number_props["house"]):
                for h in range(houses_number):
                    l["house"][p,h] = model.addVar(vtype="B", name="l_house_"+str(p)+"_"+str(h), lb=0.0, ub=1.0)

        else:
            for p in range(number_props["house"]):
                for h in range(houses_number):
                    l["house"][p,h] = model.addVar(vtype="C", name="l_house_"+str(p)+"_"+str(h), lb=0.0, ub=1.0)            
     
        for t in times:
            P_imp[t] = model.addVar(vtype="C", name="P_imp_"+str(t), lb=0.0)
            P_exp[t] = model.addVar(vtype="C", name="P_exp_"+str(t), lb=0.0)
       
        # Integrate new variables into the model
        model.update()    
    
        # Set objective
        costs_electricity = (mp.eco["b"]["el"] * mp.eco["crf"] * k_el * dt * 
                            sum(P_imp[t] for t in times) - 
                            mp.eco["b"]["eex"] * mp.eco["crf"] * r_el * dt * 
                             sum(P_exp[t] for t in times))
        costs_others = sum(sum((k[p,h] * l["house"][p,h]) for p in range(number_props["house"])) for h in range(houses_number))        

        model.setObjective(costs_electricity + costs_others, gp.GRB.MINIMIZE)
        
        # Add constraints
        # Electricity balance:
        for t in times:
            house_proposal = sum(sum(prop["house"][p][h][t] * l["house"][p,h]
                                 for p in range(number_props["house"])) for h in range(houses_number))
            model.addConstr(P_imp[t] - P_exp[t] - house_proposal == 0, "ElectricityBalance_"+str(t)) # - P_dem[d,t]

        # Convexity constraints
        for h in range(houses_number):
            model.addConstr(sum(l["house"][p,h] for p in range(number_props["house"])) == 1, "Convex_house_"+str(h))
        
        # Set Gurobi parameters
        model.Params.Presolve = 0
        if final_iteration:
            model.Params.Cuts = 3
            model.Params.MIPFocus = 1
#        model.Params.MIPGap = 0.01
        model.Params.TimeLimit = max_time
        
        # Run model
        model.optimize()
    
        # Print final solution
        if model.status == gp.GRB.OPTIMAL or model.status == gp.GRB.TIME_LIMIT:
            r = {}
            r["sigma"] = {}
            r["imp"] = {}
            r["exp"] = {}
            r_obj = model.ObjVal 
            print("Current objective of the master problem: " + str(r_obj))
            
            if final_iteration:
                r["imp"] = {(t) : P_imp[t].X for t in times}
                r["exp"] = {(t) : P_exp[t].X for t in times}
                r["house"] = np.zeros((number_props["house"], houses_number))        
                for p in range(number_props["house"]):
                    for h in range(houses_number):
                        r["house"][p,h] = l["house"][p,h].X                      
            else:
                r["pi"] = np.zeros(times_number)
                r["sigma"] = np.zeros(houses_number)
                r["imp"] = {(t) : P_imp[t].X for t in times}
                r["exp"] = {(t) : P_exp[t].X for t in times}
                                           
                for t in range(times_number):
                    r["pi"][t] = (model.getConstrByName("ElectricityBalance_"+str(t))).Pi #.getAttr("Pi")
                for h in range(houses_number):
                    r["sigma"][h]  = (model.getConstrByName("Convex_house_"+str(h))).Pi #.getAttr("Pi")
            
        else: 
            model.computeIIS()
            print('\nConstraints:')        
            for c in model.getConstrs():
                if c.IISConstr:
                    print('%s' % c.constrName)
            print('\nBounds:')
            for v in model.getVars():
                if v.IISLB > 0 :
                    print('Lower bound: %s' % v.VarName)
                elif v.IISUB > 0:
                    print('Upper bound: %s' % v.VarName)
    except gp.GurobiError:
        print('Error in masterproblem')
        
    if final_iteration:
        return (r_obj, r)
    else:
        return (r_obj, r)