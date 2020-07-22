#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 20 18:42:33 2015

@author: Thomas
"""
from __future__ import division

import xlrd
import numpy as np
import scipy.stats as stats

def read_economics(devices, filename="further_parameters.xlsx"):
    """
    Read in economic parameters and update residual values of devices.
    
    Parameters
    ----------
    devices : dictionary
        All device specific characteristics.
    filename : string, optional
        Excel-file with the economic and other parameters.
    
    Returns
    -------
    eco : dictionary
        Information on economic parameters.
    par : dictionary
        All non-economic and non-technical parameters.
    devices : dictionary
        All device specific characteristics.
    """
    book = xlrd.open_workbook(filename)
    
    sheet_eco = book.sheet_by_name("economics")
    sheet_par = book.sheet_by_name("further_parameters")
    
    eco = {}
    par = {}
    
#    kWh2J = 3600 * 1000 # Factor from kWh to Joule
    
    # Economics
    t_calc = sheet_eco.cell_value(1,1)
    eco["t_calc"] = t_calc
    eco["tax"]    = sheet_eco.cell_value(2,1)
    eco["rate"]   = sheet_eco.cell_value(3,1)
    eco["q"]      = 1 + eco["rate"]
    eco["crf"]    = ((eco["q"] ** eco["t_calc"] * eco["rate"]) / 
                     (eco["q"] ** eco["t_calc"] - 1))
    
    eco["prChange"] = {}
    eco["prChange"]["el"]   = sheet_eco.cell_value(6,1)
    eco["prChange"]["gas"]  = sheet_eco.cell_value(7,1)
    eco["prChange"]["eex"]  = sheet_eco.cell_value(8,1)
    eco["prChange"]["infl"] = sheet_eco.cell_value(9,1)

    pC = eco["prChange"]
    eco["b"] = {key: ((1 - (pC[key] / eco["q"]) ** eco["t_calc"]) / 
                      (eco["q"] - pC[key]))
                for key in pC.keys()}
    
    # Always EUR per kWh (meter per anno)
    eco["sub_chp"]        = sheet_eco.cell_value(12,1)
    eco["pr",   "el"]     = sheet_eco.cell_value(13,1)
    eco["sell", "chp"]    = sheet_eco.cell_value(14,1)
    eco["sell", "pv"]     = sheet_eco.cell_value(15,1)
    eco["gas", "chp"]     = sheet_eco.cell_value(16,1)
    eco["gas", "boiler"]  = sheet_eco.cell_value(17,1)
    eco["gas", "c_meter"] = sheet_eco.cell_value(18,1)
    
    # Determine residual values
    for dev in devices.keys():
        t_life = devices[dev]["T_op"]
        rval = (t_life - t_calc) / t_life / (eco["q"] ** t_calc)
        devices[dev]["rval"] = rval
    
    # Further parameters
    par["A_max"]           = sheet_par.cell_value(2,1)
    par["mip_gap"]         = sheet_par.cell_value(3,1)
    par["time_limit"]      = sheet_par.cell_value(4,1)
    par["rho_w"]           = sheet_par.cell_value(7,1)
    par["c_w"]             = sheet_par.cell_value(8,1)

    return (eco, par, devices)
    
def compute_parameters(par, number_clusters, len_day):
    """
    Add number of days, time steps per day and temporal discretization to par.
    
    Parameters
    ----------
    par : dictionary
        Dictionary which holds non-device-characteristic and non-economic 
        parameters.
    number_clusters : integer
        Number of allowed clusters.
    len_day : integer
        Time steps per day
    """
    par["days"] = number_clusters
    par["time_steps"] = len_day
    par["dt"] = 24 / len_day
    
    return par
    
    

def read_devices(timesteps, days, 
                 temperature_ambient, temperature_flow, temperature_design,
                 solar_irradiation,
                 filename="devices.xlsx"):
    """
    Read all devices from a given file.
    
    Parameters
    ----------
    timesteps : integer
        Number of time steps per typical day
    days : integer
        Number of typical days
    temperature_ambient : array_like
        2-dimensional array [days, timesteps] with the ambient temperature in 
        degree Celsius
    temperature_flow : float or array_like
        Required flow temperature in degree Celsius. Either as float value or
        as 2-dimensional array [days, timesteps]
    temperature_design : float
        Nominal design temperature in degree Celsius (-12 for Aachen, Germany)
    solar_irradiation : array_like
        Solar irradiation in Watt per square meter on the tilted areas on 
        which STC or PV will be installed.
    filename : string, optional
        Path to the *.xlsx file containing all available devices
    
    Return
    ------
    results : dictionary
        Dictionary containing the information for each device specified in 
        the given input file.
    """
    # Initialize results
    results = {}
    
    # Open work book
    book = xlrd.open_workbook(filename)
    
    # Get all available sheets
    available_sheets = book.sheet_names()
    
    # Iterate over all sheets
    for dev in available_sheets:
        # Read each sheet
        current_sheet = _read_sheet(book.sheet_by_name(dev), dev, timesteps)
        
        results[dev] = _handle_sheet(current_sheet, dev, timesteps, days, 
                                     temperature_ambient, temperature_flow,
                                     temperature_design,
                                     solar_irradiation)
    
    return results

def _handle_sheet(sheet, dev, timesteps, days,
                  temperature_ambient, temperature_flow, temperature_design,
                  solar_irradiation):
    """
    Parameters
    ----------
    sheet : dictionary
        Read device characteristics
    dev : string
        - `"boiler"`    : Boiler
        - `"chp"`       : CHP unit
        - `"hp"`        : Heat pump
        - `"eh"`        : Electrical heater
        - `"pv"`        : Photovoltaic modules
        - `"stc"`       : Solar thermal collectors
        - `"tes"`       : Thermal energy storage units
        - `"bat"`       : Battery units
        - `"inv"`       : Inverters    
    timesteps : integer
        Number of time steps per typical day
    days : integer
        Number of typical days
    temperature_ambient : array_like
        2-dimensional array [days, timesteps] with the ambient temperature in 
        degree Celsius
    temperature_flow : float or array_like
        Required flow temperature in degree Celsius. Either as float value or
        as 2-dimensional array [days, timesteps]
    temperature_design : float
        Nominal design temperature in degree Celsius (-12 for Aachen, Germany)
    filename : string, optional
        Path to the *.xlsx file containing all available devices
    solar_irradiation : array_like
        Solar irradiation in Watt per square meter on the tilted areas on 
        which STC or PV will be installed.
        
    Implemented characteristics
    ---------------------------
    - eta = Q/P
    - omega = (Q+P) / E
    """
    results = {}
    
    # Define infinity
    infinity = np.inf
    
    ones = np.ones((days, timesteps))
    
    keys = sheet.keys()
    
    if dev == "bat":
        capacity = np.array([sheet[i]["cap"] for i in keys])
        c_inv    = np.array([sheet[i]["c_inv"] for i in keys])
        c_om     = np.array([sheet[i]["c_om"] for i in keys])
        p_ch     = np.array([sheet[i]["P_ch_max"] for i in keys])
        p_dch    = np.array([sheet[i]["P_dch_max"] for i in keys])

        results["T_op"]   = np.mean([sheet[i]["T_op"] for i in keys])
        results["k_loss"] = np.mean([sheet[i]["k_loss"] for i in keys])
        results["eta"]    = np.mean([sheet[i]["eta"] for i in keys])

        results["cap_min"]  = np.min(capacity) # kWh
        results["cap_max"]  = np.max(capacity) # kWh
        results["c_om_rel"] = np.mean(c_om / c_inv)
        
        # Regression: p_ch = slope * capacity + intercept
        lin_reg = stats.linregress(x=capacity, y=p_ch)
        results["P_ch_fix"] = lin_reg[1]
        results["P_ch_var"] = lin_reg[0] # kW/kWh
        
        # Regression: p_dch = slope * capacity + intercept
        lin_reg = stats.linregress(x=capacity, y=p_dch)
        results["P_dch_fix"] = lin_reg[1]
        results["P_dch_var"] = lin_reg[0] # kW/kWh
        
        # Regression: c_inv = slope * capacity + intercept
        lin_reg = stats.linregress(x=capacity, y=c_inv)
        results["c_inv_fix"] = lin_reg[1]
        results["c_inv_var"] = lin_reg[0] # Euro/kWh
        
    elif dev == "boiler":
        c_inv       = np.array([sheet[i]["c_inv"] for i in keys])
        c_om        = np.array([sheet[i]["c_om"] for i in keys])
        heat_output = np.array([sheet[i]["Q_nom"] for i in keys])
        eta         = np.mean([sheet[i]["eta"] for i in keys])
                
        results["T_op"]    = np.mean([sheet[i]["T_op"] for i in keys])
        results["mod_lvl"] = np.mean([sheet[i]["mod_lvl"] for i in keys])

        results["c_om_rel"]  = np.mean(c_om / c_inv)
        results["Q_nom_min"] = np.min(heat_output)
        results["Q_nom_max"] = np.max(heat_output)
        
        # Boilers overall efficiency is equal to its thermal performance
        # and assumed constant for all times
        results["omega"] = ones * eta
        
        # The electrical efficiency (Q/P) is zero for infinite, as boilers
        # do not require a significant amount of electricity
        results["eta"] = ones * infinity
        
        # Regression: c_inv = slope * heat_output + intercept
        lin_reg = stats.linregress(x=heat_output, y=c_inv)
        results["c_inv_fix"] = lin_reg[1]
        results["c_inv_var"] = lin_reg[0]   # Euro/Watt
        
    elif dev == "chp":
        c_inv       = np.array([sheet[i]["c_inv"] for i in keys])
        c_om        = np.array([sheet[i]["c_om"] for i in keys])
        heat_output = np.array([sheet[i]["Q_nom"] for i in keys])
        eta         = np.mean([sheet[i]["eta"] for i in keys])
        omega       = np.mean([sheet[i]["omega"] for i in keys])
        
        results["T_op"]    = np.mean([sheet[i]["T_op"] for i in keys])
        results["mod_lvl"] = np.mean([sheet[i]["mod_lvl"] for i in keys])
        
        results["c_om_rel"]  = np.mean(c_om / c_inv)
        results["Q_nom_min"] = np.min(heat_output)
        results["Q_nom_max"] = np.max(heat_output)
                
        results["omega"] = ones * omega
        results["eta"]   = ones * eta
        
        # Regression: c_inv = slope * heat_output + intercept
        lin_reg = stats.linregress(x=heat_output, y=c_inv)
        results["c_inv_fix"] = lin_reg[1]
        results["c_inv_var"] = lin_reg[0]   # Euro/Watt
        
    elif dev == "eh":
        c_inv       = np.array([sheet[i]["c_inv"] for i in keys])
        c_om        = np.array([sheet[i]["c_om"] for i in keys])
        heat_output = np.array([sheet[i]["Q_nom"] for i in keys])
        eta         = np.mean([sheet[i]["eta"] for i in keys])

        results["T_op"]    = np.mean([sheet[i]["T_op"] for i in keys])
        results["mod_lvl"] = np.mean([sheet[i]["mod_lvl"] for i in keys])
        
        results["c_om_rel"]  = np.mean(c_om / c_inv)
        results["Q_nom_min"] = np.min(heat_output)
        results["Q_nom_max"] = np.max(heat_output)
        
        # Electrical heaters do not require any gas energy (E), therefore
        # omega has to be large
        results["omega"] = ones * infinity
        
        # The effiency is defined like eta (Q/P)
        results["eta"]   = ones * eta
        
        # Regression: c_inv = slope * heat_output + intercept
        lin_reg = stats.linregress(x=heat_output, y=c_inv)
        results["c_inv_fix"] = lin_reg[1]
        results["c_inv_var"] = lin_reg[0]   # Euro/Watt
        
    elif dev == "hp":
        c_inv       = np.array([sheet[i]["c_inv"] for i in keys])
        c_om        = np.array([sheet[i]["c_om"] for i in keys])
        heat_output = np.array([sheet[i]["Q_nom"] for i in keys])
        
        results["c_om_rel"]  = np.mean(c_om / c_inv)
        results["Q_nom_min"] = np.min(heat_output)
        results["Q_nom_max"] = np.max(heat_output)
        
        results["T_op"]    = np.mean([sheet[i]["T_op"] for i in keys])
        results["mod_lvl"] = np.mean([sheet[i]["mod_lvl"] for i in keys])
        results["dT_max"]  = np.mean([sheet[i]["dT_max"] for i in keys])
        
        results["cop_a2w35"] = np.mean([sheet[i]["cop_a2w35"] for i in keys])
        results["P_nom_min"] = results["Q_nom_min"] / results["cop_a2w35"]
        results["P_nom_max"] = results["Q_nom_max"] / results["cop_a2w35"]
        
        # Define function to compute the reversible coefficient of 
        # performance (COP) based on ambient temperature and required flow
        # temperature (both given in degree Celsius).
        cop_rev = lambda t_amb, t_flow: (t_flow+273.15) / (t_flow - t_amb)

        # Read nominal COP
        cop_nom = results["cop_a2w35"]
        
        # Compute efficiency of the given device: zeta = cop_nom / cop_rev
        zeta = cop_nom / cop_rev(2, 35)
        
        # Assume that zeta is constant for all temperatures
        cop = cop_rev(temperature_ambient, temperature_flow) * zeta

        # Compute COP at design conditions
        cop_design = cop_rev(temperature_design, temperature_flow) * zeta
        results["cop_design"] = cop_design
        
        # Heat pumps do not require any gas energy (E), therefore omega 
        # has to be a large value
        results["omega"] = ones * infinity
        
        # The COP value is defined like eta (Q/P)
        results["eta"]   = ones * cop
        
        # Regression: c_inv = slope * heat_output + intercept
        lin_reg = stats.linregress(x=heat_output, y=c_inv)
        results["c_inv_fix"] = lin_reg[1]
        results["c_inv_var"] = lin_reg[0]   # Euro/Watt
        
    elif dev == "inv":
        c_inv    = np.array([sheet[i]["c_inv"] for i in keys])
        c_om     = np.array([sheet[i]["c_om"] for i in keys])
        power_dc = np.array([sheet[i]["P_nom_DC"] for i in keys])
        
        results["c_om_rel"]  = np.mean(c_om / c_inv)
        results["power_min"] = np.min(power_dc)
        results["power_max"] = np.max(power_dc)
        
        results["T_op"] = np.mean([sheet[i]["T_op"] for i in keys])
        results["eta"]  = np.mean([sheet[i]["eta"] for i in keys])
        
        # Regression: c_inv = slope * power_dc + intercept
        lin_reg = stats.linregress(x=power_dc, y=c_inv)
        results["c_inv_fix"] = lin_reg[1]
        results["c_inv_var"] = lin_reg[0]   # Euro/Watt
        
    elif dev == "pv":
        c_inv = np.array([sheet[i]["c_inv"] for i in keys])
        c_om  = np.array([sheet[i]["c_om"] for i in keys])
        area  = np.array([sheet[i]["area"] for i in keys])
        
        results["c_om_rel"]  = np.mean(c_om / c_inv)
        results["area_mean"] = np.mean(area)
        results["area_min"]  = np.min(area)
        
        results["T_op"]   = np.mean([sheet[i]["T_op"] for i in keys])
        results["p_NOCT"] = np.mean([sheet[i]["p_NOCT"] for i in keys])
        results["t_NOCT"] = np.mean([sheet[i]["t_NOCT"] for i in keys])
        results["gamma"]  = np.mean([sheet[i]["gamma"] for i in keys])
                
        i_NOCT = 0.8 # kW / m2
        
        # Interpolate cell temperature.
        # Without solar irradiation, the cell temperature has to be equal
        # to the ambient temperature. At NOCT irradiation, the cell's 
        # temperature has to be equal to t_NOCT
        t_cell = (temperature_ambient + solar_irradiation / i_NOCT * 
                                  (results["t_NOCT"] - temperature_ambient))
        eta_NOCT = results["p_NOCT"] / (results["area_mean"] * i_NOCT)
        # Compute electrical efficiency of the cell
        eta_el   = eta_NOCT * (1 + results["gamma"] / 100 * 
                               (t_cell - results["t_NOCT"]))
        
        results["eta_th"] = np.zeros((days, timesteps))
        results["eta_el"] = eta_el
        
        results["c_inv_fix"] = 0
        results["c_inv_var"] = np.mean(c_inv) / results["area_mean"]  # Euro/m2
        
    elif dev == "stc":
        c_inv = np.array([sheet[i]["c_inv"] for i in keys])
        c_om  = np.array([sheet[i]["c_om"] for i in keys])
        area  = np.array([sheet[i]["area"] for i in keys])

        results["c_om_rel"]  = np.mean(c_om / c_inv)
        results["area_mean"] = np.mean(area)
        results["area_min"]  = np.min(area)
        
        results["T_op"] = np.mean([sheet[i]["T_op"] for i in keys])

        zero_loss    = np.mean([sheet[i]["zero_loss"] for i in keys])
        first_order  = np.mean([sheet[i]["first_order"] for i in keys])
        second_order = np.mean([sheet[i]["second_order"] for i in keys])
        
        results["eta_el"] = np.zeros((days, timesteps))
        temp_diff = temperature_flow - temperature_ambient
        eta_th = (zero_loss - 
                  first_order / solar_irradiation * temp_diff - 
                  second_order / solar_irradiation * (temp_diff**2))
        results["eta_th"] = np.maximum(eta_th, 0)
        
        results["c_inv_fix"] = 0
        results["c_inv_var"] = np.mean(c_inv) / results["area_mean"]  # Euro/m2
        
    elif dev == "tes":
        c_inv  = np.array([sheet[i]["c_inv"] for i in keys])
        c_om   = np.array([sheet[i]["c_om"] for i in keys])
        volume = np.array([sheet[i]["volume"] for i in keys])
        
        results["c_om_rel"]   = np.mean(c_om / c_inv)
        results["volume_min"] = np.min(volume)
        results["volume_max"] = np.max(volume)
        
        results["T_op"]    = np.mean([sheet[i]["T_op"] for i in keys])
        results["k_loss"]  = np.mean([sheet[i]["k_loss"] for i in keys])
        results["eta_ch"]  = np.mean([sheet[i]["eta_ch"] for i in keys])
        results["eta_dch"] = np.mean([sheet[i]["eta_dch"] for i in keys])
                
        results["dT_max"] = 40  # Max. temperature spread in Kelvin
        
        # Regression: c_inv = slope * volume + intercept
        lin_reg = stats.linregress(x=volume, y=c_inv)
        results["c_inv_fix"] = lin_reg[1]
        results["c_inv_var"] = lin_reg[0]   # Euro/m3
    return results

def _read_sheet(sheet, device, timesteps):
    """
    sheet : sheet-object
        Sheet of the workbook containing all available devices
    device : string
        - `"boiler"`    : Boiler
        - `"chp"`       : CHP unit
        - `"hp"`        : Heat pump
        - `"eh"`        : Electrical heater
        - `"pv"`        : Photovoltaic modules
        - `"stc"`       : Solar thermal collectors
        - `"tes"`       : Thermal energy storage units
        - `"bat"`       : Battery units
        - `"inv"`       : Inverters
    timesteps : integer
        Number of time steps per typical day
    
    Implemented characteristics
    ---------------------------
    - eta = Q/P
    - omega = (Q+P) / E
    """
    
    # Initialize results
    results = {}
    
    # Read all rows but the headers:
    for row in range(1, sheet.nrows):
        # Create new dictionary for current entry. Add common inputs.
        current_results = {}
        # Handle each device separately
        if device == "bat":
            current_results["c_inv"]     = sheet.cell_value(row, 1)
            current_results["c_om"]      = sheet.cell_value(row, 2)
            current_results["T_op"]      = sheet.cell_value(row, 3)
            current_results["cap"]       = sheet.cell_value(row, 4)
            current_results["eta"]       = sheet.cell_value(row, 5)
            current_results["P_ch_max"]  = sheet.cell_value(row, 6)
            current_results["P_dch_max"] = sheet.cell_value(row, 7)
            current_results["k_loss"]    = 0

        elif device == "boiler":
            current_results["Q_nom"]   = sheet.cell_value(row, 1)
            current_results["mod_lvl"] = sheet.cell_value(row, 2)
            current_results["c_inv"]   = sheet.cell_value(row, 3)
            current_results["c_om"]    = sheet.cell_value(row, 4)
            current_results["T_op"]    = sheet.cell_value(row, 5)
            current_results["eta"]     = sheet.cell_value(row, 6)

        elif device == "chp":
            current_results["Q_nom"]   = sheet.cell_value(row, 1)
            current_results["mod_lvl"] = sheet.cell_value(row, 2)
            current_results["c_inv"]   = sheet.cell_value(row, 3)
            current_results["c_om"]    = sheet.cell_value(row, 4)
            current_results["T_op"]    = sheet.cell_value(row, 5)
            current_results["eta"]     = 1 / sheet.cell_value(row, 6)
            current_results["omega"]   = sheet.cell_value(row, 7)
            
        elif device == "eh":
            current_results["Q_nom"]   = sheet.cell_value(row, 1)
            current_results["mod_lvl"] = sheet.cell_value(row, 2)
            current_results["c_inv"]   = sheet.cell_value(row, 3)
            current_results["c_om"]    = sheet.cell_value(row, 4)
            current_results["T_op"]    = sheet.cell_value(row, 5)
            current_results["eta"]     = sheet.cell_value(row, 6)
            
        elif device == "hp":
            current_results["Q_nom"]   = sheet.cell_value(row, 1)
            current_results["mod_lvl"] = sheet.cell_value(row, 2)
            current_results["c_inv"]   = sheet.cell_value(row, 3)
            current_results["c_om"]    = sheet.cell_value(row, 4)
            current_results["T_op"]    = sheet.cell_value(row, 5)
            current_results["dT_max"]  = sheet.cell_value(row, 7)
            current_results["cop_a2w35"] = sheet.cell_value(row, 6)
            
        elif device == "inv":
            current_results["P_nom_DC"]  = sheet.cell_value(row, 1)
            current_results["c_inv"]     = sheet.cell_value(row, 2)
            current_results["c_om"]      = sheet.cell_value(row, 3)
            current_results["T_op"]      = sheet.cell_value(row, 4)
            current_results["eta"]       = sheet.cell_value(row, 5)
            
        elif device == "pv":
            current_results["c_inv"] = sheet.cell_value(row, 2)
            current_results["c_om"]  = sheet.cell_value(row, 3)
            current_results["T_op"]  = sheet.cell_value(row, 4)
            current_results["area"]  = sheet.cell_value(row, 5)
            
            current_results["p_NOCT"]  = sheet.cell_value(row, 1)
            current_results["t_NOCT"]  = sheet.cell_value(row, 6)
            current_results["gamma"]   = sheet.cell_value(row, 7)

        elif device == "stc":
            current_results["c_inv"] = sheet.cell_value(row, 1)
            current_results["c_om"]  = sheet.cell_value(row, 2)
            current_results["T_op"]  = sheet.cell_value(row, 3)
            current_results["area"]  = sheet.cell_value(row, 4)

            current_results["zero_loss"]    = sheet.cell_value(row, 5)
            current_results["first_order"]  = sheet.cell_value(row, 6)
            current_results["second_order"] = sheet.cell_value(row, 7)

        elif device == "tes":
            current_results["c_inv"]    = sheet.cell_value(row, 1)
            current_results["c_om"]     = sheet.cell_value(row, 2)
            current_results["T_op"]     = sheet.cell_value(row, 3)
            current_results["eta_ch"]   = sheet.cell_value(row, 6)
            current_results["eta_dch"]  = sheet.cell_value(row, 7)
            
            standby_losses = sheet.cell_value(row, 4) # in kWh / 24h
            volume = sheet.cell_value(row, 5)         # in m3
            
            temp_diff_norm = 45 # K
            heat_cap       = 4180 # J/(kgK)
            density        = 1000 # kg/m3
            energy_content = volume * temp_diff_norm * heat_cap * density # J
            
            k_loss_day = 1 - standby_losses / energy_content * 3600*1000
            current_results["k_loss"] = 1 - (k_loss_day ** (1 / timesteps))
            current_results["volume"] = volume
            
        results[row] = current_results
        
    return results

if __name__ == "__main__":
    timesteps = 24
    days = 5
    # Random temperatures between -10 and +20 degC:
    temperature_ambient = np.random.rand(days, timesteps) * 30 - 10
    
    temperature_design = -12 # Aachen
    
    solar_irradiation = np.random.rand(days, timesteps) * 800
    solar_irradiation
    
    devs = read_devices(timesteps, days, temperature_ambient, 
                        temperature_flow=35,
                        temperature_design=temperature_design,
                        solar_irradiation=solar_irradiation)
                        
    (eco, par, devs) = read_economics(devs)
