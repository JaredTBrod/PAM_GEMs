import cobra
from cobra.core.metabolite import elements_and_molecular_weights
elements_and_molecular_weights['R']=0.0
elements_and_molecular_weights['Z']=0.0
import pandas as pd
import numpy as np
import csv

#### Change Biomass composition
# define a function change a biomass reaction in the model
def update_biomass(model, rxn, stoich, metabolite):
    r = model.reactions.get_by_id(rxn)

    new_stoich = stoich

    # you now have a dictionary of new stoichs for your model
    for m,s in r.metabolites.items():
        stoich = s*-1
        temp_dict = {m:stoich}
        r.add_metabolites(temp_dict)
    r.add_metabolites(new_stoich)

    # Then get the total to equal 1 mg biomass DW
    total = 0
    for m,s in r.metabolites.items():
        gfw = model.metabolites.get_by_id(m.id).formula_weight
        mass = gfw*s*-1
        total = total+mass
    correction = total/1000 # this will get it to 1000 ug total mass

    # Then adjust the stoichiometry as appropriate
    for m,s in r.metabolites.items(): # now change the stoich
        to_add = s/correction-s
        r.add_metabolites({m:to_add})

    # Finally build the biomass_c metabolite
        imbal = r.check_mass_balance()
        if 'charge' in imbal.keys():
            met_charge = imbal['charge']*-1
            del imbal['charge']
        met_mass = 0
        formula_string = ''
        for e,v in imbal.items():
            if v > 1e-10 or v < -1e-10:
                mass = elements_and_molecular_weights[e]
                met_mass = met_mass+(mass*-1*v)
                form_str = e+str(-1*v)
                formula_string = formula_string + form_str

    met = model.metabolites.get_by_id(metabolite)
    met.formula = formula_string
    met.charge = met_charge
    r.add_metabolites({met:1})

    # Add GAM constraint
    if rxn == 'bof_c':
        gam_met = model.metabolites.GAM_const_c
        r.add_metabolites({gam_met:1})
    
    model.repair()

#     print(model.reactions.get_by_id(rxn).reaction)
#     print('')
#     print(model.metabolites.get_by_id(met.id).formula)
#     print('')
#     print(model.metabolites.get_by_id(met.id).formula_weight)
#     print('')
#     print(model.reactions.get_by_id(rxn).check_mass_balance())
    return model

#############################################
#### Simulate growth with all constraints####
#############################################
def figure_2LL(model):
    ## Run all 27 parameter combos and capture:
    #### - growth rate
    #### - biomass
    import math as m
    cols = ['ID','GR','mgDW','Cells']
    data_out = pd.DataFrame(columns=cols)

    Po_vals = ['mean','ub','lb']
    a_star_vals = ['mean','ub','lb']
    cell_mass = ['mean','ub','lb'] 

    a_star_all = pd.read_csv('pt_a_star.csv',index_col=0,header=0)
    rel = pd.read_csv('fluor_lamp.csv',header=0)

    xsec = 80. #culture cross sectional area in cm2
    path_len = 4.7 #cm
    volume = 375.
    time_interval = 20

    for p in Po_vals:
        for a in a_star_vals:
            for c in cell_mass:
                sample_id = p+'_'+a+'_'+c
                model = model
                if a == 'mean':
                    a_star = a_star_all['LL'] 
                elif a == 'ub':
                    a_star = a_star_all['LL_UB']
                else: 
                    a_star = a_star_all['LL_LB']

                
                if c == 'mean':
                    mgCell = 19.1/1e9
                elif c == 'ub':
                    mgCell = 21.1/1e9
                else: 
                    mgCell = 17.1/1e9
                
                innoc = 2.8e6*volume # cells/mL * total mL 
                iDW = mgCell*innoc # initial culture biomass
                gDW = iDW
                cells = (gDW/mgCell)

                biomass = pd.DataFrame(data=[gDW,cells],index = ['Biomass','Cells'],columns=['0']) # add it up at the end of each simulation

                for t in range(time_interval,1440+time_interval,time_interval): #One simulation every 20 minutes from t=0 to t=24 hrs (1440 min)
                    import math as m
                    interval_bm = 0 #initializes the total biomass
                    photon_count = 0
                    atten = np.zeros(len(a_star)) #captures light attenuation

                    tmp_o2_evo = 0

                    gDW = biomass[str(t-time_interval)]['Biomass'] #mg DW
                    cells = (gDW/mgCell)

                    if p == 'mean':
                        Ps = 2.94e-11 # umol O2 cell-1 s-1
                        alpha = 9.82e-2
                        beta = 0.0
                        resp = 1.83e-12
                    
                    elif p == 'ub':
                        Ps = 3.06e-11 # umol O2 cell-1 s-1
                        alpha = 9.75e-2
                        beta = 0.0
                        resp = 1.75e-12
      
                    else: 
                        Ps = 2.83e-11 # umol O2 cell-1 s-1
                        alpha = 9.91e-2
                        beta = 0.0
                        resp = 1.92e-12

                    irrad = 60.                       

                    # Photon_flux is the initial amount of light delivered to the culture at each wavelength from 400-700 nm
                    photon_flux = (rel['rel_height'].values)*irrad*xsec/10000*time_interval*60. #umol/(m2*s) * cm2 * 1m2/10000cm2 * 60s/min * min = umol photons/time interval
                    total_photons = sum(photon_flux)

                    total_absorbed = 0.

                    ti_o2 = 0.
                    for nm in range(len(photon_flux)):

                        abs_coeff = a_star[400+nm]              # a* value for the given nm (cm2/cell)
                        Io = photon_flux[nm]    # incident photon flux at this nm at this slice (umol_photon/time_interval)
                        Ia = Io-Io*(m.exp(-1*abs_coeff*cells/volume*path_len))
                        nm_abs = Ia

                        total_absorbed = total_absorbed+nm_abs
                    conv_abs = total_absorbed/time_interval/60./(cells) # converts abs to O2 evo curve units  umol/TI * TI/min * min/s * 1/cells => umol/(mgchla*s)
                    slice_o2 = (Ps*(1-m.exp(-1*alpha*conv_abs/Ps))*(m.exp(-1*beta*conv_abs/Ps)))-resp #umol O2 cell-1 s-1
                    ti_o2 = ti_o2+(slice_o2*(cells)*60.*time_interval) #umol O2 cell-1 s-1 * cells * s/min * min/TI = umol O2/TI

                    o2evo = ti_o2
                    
                    model.reactions.EX_o2_e.bounds = (0.99*o2evo,o2evo)#                            <----Po Constraint
                    
                    ngam = resp*60.*time_interval*cells
                    model.reactions.NGAM.lower_bound = ngam
                    cef_rate = 1.1*gDW/60.*time_interval #CEF sets CEF_h upper bound.  umol /(mgDW*h) * mgDW * h/60min * min/TI = umol/TI
                    model.reactions.CEF_h.upper_bound = cef_rate
                                        
                    model.reactions.EX_photon_e.bounds = (total_absorbed*-1.,total_absorbed*-0.9999)

                ###### Parameters for PSII fluorescence

                ## Fv/Fm
                    FvFm_LL = 0.69

                ## Calculate Y(II) based on absorbed 
                    abs_conv = total_absorbed/time_interval/60./(cells)

                    yII_LL = 0.7016*np.exp(-8.535e8*abs_conv)
                #     yII_LL2 = (2.77e16*(m.pow(abs_conv,2))-(2.51e8*abs_conv)+6.97e-1)

                ## Y(NPQ)
                    yNPQ_LL = 1./(1.+(np.power(2.3e-9/abs_conv,3.)))
                    if yNPQ_LL < 0:
                        yNPQ_LL = 0.     

                ### Constraints

                    phoYII = round(yII_LL,2) # Y(II)

                    regNPQ = round(yNPQ_LL,2) # Y(NPQ)

                    unrNPQ = round((1-phoYII-regNPQ),2) # Y(NO)

                ### Edit model constraints with the appropriate values

                ## PHO_PSIIt_u

                    rxn = model.reactions.PHO_PSIIt_u

                    ## reset the stoich
                    for met,s in rxn.metabolites.items():
                        stoich = s*-1
                        temp_dict = {met:stoich}
                        rxn.add_metabolites(temp_dict)

                    m1 = model.metabolites.photon_YII_u
                    m2 = model.metabolites.photon_YNPQ_u
                    m4 = model.metabolites.photon_YNO_u
                    m3 = model.metabolites.photon_h

                    rxn.add_metabolites({m1:phoYII,
                                        m2:regNPQ,
                                        m4:unrNPQ,
                                        m3:-1.})

                    model.reactions.DM_photon_c.upper_bound = 0. # constrained                                    <----hv Constraint

                    # Add D1 damage cost uncoupled to PSII
                    ## If uncoupled, set the lower and upper bounds to the experimentally determined values
                    # Damage rates
                    D1_rate = 7e-6# # LL: umol D1/mgDW h-1         <---- D1 Constraint
                    D1_rate = D1_rate * gDW/60.# ugD1 ugDW-1 min-1 * mgDW * 1h/60 min = umolD1 min-1
                    D1_rate = D1_rate * time_interval # umolD1 min-1 * min/TI = umol D1/TI
                    
                    model.reactions.NGAM_D1_u.bounds = (D1_rate, 1.001*(D1_rate))                    

                    ## Solve the model
                    model.objective = 'bof_c'
                    solution = model.optimize() 
                    if solution.status == 'optimal':
                        obj_rxn = model.reactions.bof_c
                        biomass[str(t)]=(gDW+obj_rxn.x,(gDW+obj_rxn.x)/(mgCell))

                    ### collect data
                    if t == 1440:
                        dry_weight = biomass['1420']['Biomass']
                        cell_count = biomass['1420']['Cells']
                        mu = np.log(biomass['1440']['Cells']/biomass['0']['Cells'])/(1440/60)
                        data_out = data_out.append({'ID':sample_id,'GR':mu,'mgDW':dry_weight,'Cells':cell_count},ignore_index=True)
    return data_out

def figure_2HL(model):
    ## Run all 27 parameter combos and capture:
    #### - growth rate
    #### - biomass
    import math as m
    cols = ['ID','GR','mgDW','Cells']
    data_out = pd.DataFrame(columns=cols)

    Po_vals = ['mean','ub','lb']
    a_star_vals = ['mean','ub','lb']
    cell_mass = ['mean','ub','lb'] 

    a_star_all = pd.read_csv('pt_a_star.csv',index_col=0,header=0)
    rel = pd.read_csv('fluor_lamp.csv',header=0)

    xsec = 80. #culture cross sectional area in cm2
    path_len = 4.7 #cm
    volume = 375.
    time_interval = 20

    for p in Po_vals:
        for a in a_star_vals:
            for c in cell_mass:
                sample_id = p+'_'+a+'_'+c
                model = model
                if a == 'mean':
                    a_star = a_star_all['HL'] 
                elif a == 'ub':
                    a_star = a_star_all['HL_UB']
                else: 
                    a_star = a_star_all['HL_LB']

                    
                if c == 'mean':
                    mgCell = 20.4/1e9
                elif c == 'ub':
                    mgCell = 21.8/1e9
                else: 
                    mgCell = 19.0/1e9
                
                innoc = 3.5e6*volume # cells/mL * total mL 
                iDW = mgCell*innoc # initial culture biomass
                gDW = iDW
                cells = (gDW/mgCell)

                biomass = pd.DataFrame(data=[gDW,cells],index = ['Biomass','Cells'],columns=['0']) # add it up at the end of each simulation

                for t in range(time_interval,720+time_interval,time_interval): #One simulation every 20 minutes from t=0 to t=24 hrs (1440 min)
                    import math as m
                    interval_bm = 0 #initializes the total biomass
                    photon_count = 0
                    atten = np.zeros(len(a_star)) #captures light attenuation

                    tmp_o2_evo = 0

                    gDW = biomass[str(t-time_interval)]['Biomass'] #mg DW
                    cells = (gDW/mgCell)

                    if p == 'mean':
                        Ps = 2.46e-11 # umol O2 cell-1 s-1
                        alpha = 9.47e-2
                        beta = 0.0
                        resp = 4.72e-12
                    
                    elif p == 'ub':
                        Ps = (2.54e-11) # umol O2 cell-1 s-1
                        alpha = (1.05e-1)
                        beta = 3.67e-5
                        resp = 4.28e-12 
      
                    else: 
                        Ps = (2.35e-11) # umol O2 cell-1 s-1
                        alpha = (8.5e-2)#
                        beta = 2.5e-10
                        resp = 5.15e-12

                 # Calculate photons and O2 uptake
                    irrad =  600. #uE 

                    # Photon_flux is the initial amount of light delivered to the culture at each wavelength from 400-700 nm
                    photon_flux = (rel['rel_height'].values)*irrad*xsec/10000*time_interval*60. #umol/(m2*s) * cm2 * 1m2/10000cm2 * 60s/min * min = umol photons/time interval
                    total_photons = sum(photon_flux)

                    total_absorbed = 0.

                    ti_o2 = 0.
                    for nm in range(len(photon_flux)):

                        abs_coeff = a_star[400+nm]              # a* value for the given nm (cm2/cell)
                        Io = photon_flux[nm]    # incident photon flux at this nm at this slice (umol_photon/time_interval)
                        Ia = Io-Io*(m.exp(-1*abs_coeff*cells/volume*path_len))
                        nm_abs = Ia

                        total_absorbed = total_absorbed+nm_abs
                    conv_abs = total_absorbed/time_interval/60./(cells) # converts abs to O2 evo curve units  umol/TI * TI/min * min/s * 1/cells => umol/(mgchla*s)
                    slice_o2 = (Ps*(1-m.exp(-1*alpha*conv_abs/Ps))*(m.exp(-1*beta*conv_abs/Ps)))-resp #umol O2 cell-1 s-1
                    ti_o2 = ti_o2+(slice_o2*(cells)*60.*time_interval) #umol O2 cell-1 s-1 * cells * s/min * min/TI = umol O2/TI

                    o2evo = ti_o2

                    model.reactions.EX_o2_e.bounds = (0.99*o2evo,o2evo)
                    
                    ngam = resp*60.*time_interval*cells
                    model.reactions.AOX_m.lower_bound = ngam
                    cef_rate = 0.7*gDW/60.*time_interval #CEF sets CEF_h upper bound.  umol /(mgDW*h) * mgDW * h/60min * min/TI = umol/TI
                    model.reactions.CEF_h.upper_bound = cef_rate
                    
                    model.reactions.EX_photon_e.bounds = (total_absorbed*-1.,total_absorbed*-0.9999)

                ###### Parameters for PSII fluorescence

                ## Fv/Fm
                    FvFm_HL = 0.63

                ## Calculate Y(II) based on absorbed 
                    abs_conv = total_absorbed/time_interval/60./(cells)

                    yII_HL = 0.6398*np.exp(-1.169e9*abs_conv)

                ## Y(NPQ)
                    yNPQ_HL = 1./(1.+(np.power(1.37e-9/abs_conv,5.5)))
                    if yNPQ_HL < 0:
                        yNPQ_HL = 0.                                            


                ### Constraints 

                    phoYII = round(yII_HL,2) # Y(II)

                    regNPQ = round(yNPQ_HL,2) # Y(NPQ)

                    unrNPQ = round((1-phoYII-regNPQ),2) # Y(NO)

                ### Edit model constraints with the appropriate values

                ## Set Y(II) constraints   
                    rxn = model.reactions.PHO_PSIIt_u

                    ## reset the stoich
                    for met,s in rxn.metabolites.items():
                        stoich = s*-1
                        temp_dict = {met:stoich}
                        rxn.add_metabolites(temp_dict)

                    m1 = model.metabolites.photon_YII_u
                    m2 = model.metabolites.photon_YNPQ_u
                    m4 = model.metabolites.photon_YNO_u
                    m3 = model.metabolites.photon_h

                    rxn.add_metabolites({m1:phoYII,
                                        m2:regNPQ,
                                        m4:unrNPQ,
                                        m3:-1.})
                ## Set photon constraints
                    model.reactions.DM_photon_c.upper_bound = 0. # constrained

                    # Add D1 damage cost uncoupled to PSII
                   # Damage rates
                    D1_rate = 2.52e-4# # HL: umol D1/mgDW h-1         <---- D1 Constraint
                    D1_rate = D1_rate * gDW/60.# ugD1 ugDW-1 min-1 * mgDW * 1h/60 min = umolD1 min-1
                    D1_rate = D1_rate * time_interval # umolD1 min-1 * min/TI = umol D1/TI
                    
                    model.reactions.NGAM_D1_u.bounds = (D1_rate, 1.0001*(D1_rate))
                    
                    ## Solve the model
                    model.objective = 'bof_c'
                    solution = model.optimize() 
                    if solution.status == 'optimal':
                        obj_rxn = model.reactions.bof_c

                        biomass[str(t)]=(gDW+obj_rxn.x,(gDW+obj_rxn.x)/(mgCell))

                       
                                        ### collect data
                        if t == 720:
                            dry_weight = biomass['700']['Biomass']
                            cell_count = biomass['700']['Cells']
                            mu = np.log(biomass['720']['Cells']/biomass['0']['Cells'])/(720/60)
                           
                            data_out = data_out.append({'ID':sample_id,'GR':mu,'mgDW':dry_weight,'Cells':cell_count},ignore_index=True)
    return data_out

def figure_2HL_HCO3(model):
    ## Run all 27 parameter combos and capture:
    #### - growth rate
    #### - biomass
    import math as m
    cols = ['ID','GR','mgDW','Cells']
    data_out = pd.DataFrame(columns=cols)

    Po_vals = ['mean','ub','lb']
    a_star_vals = ['mean','ub','lb']
    cell_mass = ['mean','ub','lb'] 

    a_star_all = pd.read_csv('pt_a_star.csv',index_col=0,header=0)
    rel = pd.read_csv('fluor_lamp.csv',header=0)

    xsec = 80. #culture cross sectional area in cm2
    path_len = 4.7 #cm
    volume = 375.
    time_interval = 20

    for p in Po_vals:
        for a in a_star_vals:
            for c in cell_mass:
                sample_id = p+'_'+a+'_'+c
                model = model
                if a == 'mean':
                    a_star = a_star_all['HL'] 
                elif a == 'ub':
                    a_star = a_star_all['HL_UB']
                else: 
                    a_star = a_star_all['HL_LB']

                    
                if c == 'mean':
                    mgCell = 20.4/1e9
                elif c == 'ub':
                    mgCell = 21.8/1e9
                else: 
                    mgCell = 19.0/1e9
                
                innoc = 3.5e6*volume # cells/mL * total mL 
                iDW = mgCell*innoc # initial culture biomass
                gDW = iDW
                cells = (gDW/mgCell)

                biomass = pd.DataFrame(data=[gDW,cells],index = ['Biomass','Cells'],columns=['0']) # add it up at the end of each simulation

                for t in range(time_interval,720+time_interval,time_interval): #One simulation every 20 minutes from t=0 to t=24 hrs (1440 min)
                    import math as m
                    interval_bm = 0 #initializes the total biomass
                    photon_count = 0
                    atten = np.zeros(len(a_star)) #captures light attenuation

                    tmp_o2_evo = 0

                    gDW = biomass[str(t-time_interval)]['Biomass'] #mg DW
                    cells = (gDW/mgCell)

                    if p == 'mean':
                        Ps = (2.46e-11)*1.15 # umol O2 cell-1 s-1
                        alpha = (9.47e-2)
                        beta = 0.0
                        resp = 4.72e-12
                    
                    elif p == 'ub':
                        Ps = (2.54e-11)*1.15 # umol O2 cell-1 s-1
                        alpha = (1.05e-1)
                        beta = 3.67e-5
                        resp = 4.28e-12 
      
                    else: 
                        Ps = (2.35e-11)*1.15 # umol O2 cell-1 s-1
                        alpha = (8.5e-2)#
                        beta = 2.5e-10
                        resp = 5.15e-12

                 # Calculate photons and O2 uptake
                    irrad =  600. #uE 

                    # Photon_flux is the initial amount of light delivered to the culture at each wavelength from 400-700 nm
                    photon_flux = (rel['rel_height'].values)*irrad*xsec/10000*time_interval*60. #umol/(m2*s) * cm2 * 1m2/10000cm2 * 60s/min * min = umol photons/time interval
                    total_photons = sum(photon_flux)

                    total_absorbed = 0.

                    ti_o2 = 0.
                    for nm in range(len(photon_flux)):

                        abs_coeff = a_star[400+nm]              # a* value for the given nm (cm2/cell)
                        Io = photon_flux[nm]    # incident photon flux at this nm at this slice (umol_photon/time_interval)
                        Ia = Io-Io*(m.exp(-1*abs_coeff*cells/volume*path_len))
                        nm_abs = Ia

                        total_absorbed = total_absorbed+nm_abs
                    conv_abs = total_absorbed/time_interval/60./(cells) # converts abs to O2 evo curve units  umol/TI * TI/min * min/s * 1/cells => umol/(mgchla*s)
                    slice_o2 = (Ps*(1-m.exp(-1*alpha*conv_abs/Ps))*(m.exp(-1*beta*conv_abs/Ps)))-resp #umol O2 cell-1 s-1
                    ti_o2 = ti_o2+(slice_o2*(cells)*60.*time_interval) #umol O2 cell-1 s-1 * cells * s/min * min/TI = umol O2/TI

                    o2evo = ti_o2

                    model.reactions.EX_o2_e.bounds = (0.99*o2evo, o2evo)

                    ngam = resp*60.*time_interval*cells
                    model.reactions.AOX_m.lower_bound = ngam
                    cef_rate = 0.73*gDW/60.*time_interval #CEF sets CEF_h upper bound.  umol /(mgDW*h) * mgDW * h/60min * min/TI = umol/TI
                    model.reactions.CEF_h.upper_bound = cef_rate
                    
                    model.reactions.EX_photon_e.bounds = (total_absorbed*-1.,total_absorbed*-0.9999)

                ###### Parameters for PSII fluorescence

                ## Fv/Fm
                    FvFm_HL = 0.63

                ## Calculate Y(II) based on absorbed 
                    abs_conv = total_absorbed/time_interval/60./(cells)

                    yII_HL = 0.6398*np.exp(-1.169e9*abs_conv)

                ## Y(NPQ)
                    yNPQ_HL = 1./(1.+(np.power(1.37e-9/abs_conv,5.5)))
                    if yNPQ_HL < 0:
                        yNPQ_HL = 0.                                            


                ### Constraints 

                    phoYII = round(yII_HL,2) # Y(II)

                    regNPQ = round(yNPQ_HL,2) # Y(NPQ)

                    unrNPQ = round((1-phoYII-regNPQ),2) # Y(NO)

                ### Edit model constraints with the appropriate values

                ## Set Y(II) constraints   
                    rxn = model.reactions.PHO_PSIIt_u

                    ## reset the stoich
                    for met,s in rxn.metabolites.items():
                        stoich = s*-1
                        temp_dict = {met:stoich}
                        rxn.add_metabolites(temp_dict)

                    m1 = model.metabolites.photon_YII_u
                    m2 = model.metabolites.photon_YNPQ_u
                    m4 = model.metabolites.photon_YNO_u
                    m3 = model.metabolites.photon_h

                    rxn.add_metabolites({m1:phoYII,
                                        m2:regNPQ,
                                        m4:unrNPQ,
                                        m3:-1.})
                ## Set photon constraints
                    model.reactions.DM_photon_c.upper_bound = 0. # constrained

                    # Add D1 damage cost uncoupled to PSII
                   # Damage rates
                    D1_rate = 2.52e-4# # LL: umol D1/mgDW h-1         <---- D1 Constraint
                    D1_rate = D1_rate * gDW/60.# ugD1 ugDW-1 min-1 * mgDW * 1h/60 min = umolD1 min-1
                    D1_rate = D1_rate * time_interval # umolD1 min-1 * min/TI = umol D1/TI
                    
                    model.reactions.NGAM_D1_u.bounds = (D1_rate,1.0001*(D1_rate))

                    ## Solve the model
                    model.objective = 'bof_c'
                    solution = model.optimize()
                    if solution.status == 'optimal':
                        obj_rxn = model.reactions.bof_c

                        biomass[str(t)]=(gDW+obj_rxn.x,(gDW+obj_rxn.x)/(mgCell))

                       
                                        ### collect data
                        if t == 720:
                            dry_weight = biomass['700']['Biomass']
                            cell_count = biomass['700']['Cells']
                            mu = np.log(biomass['720']['Cells']/biomass['0']['Cells'])/(720/60)
                            
                            data_out = data_out.append({'ID':sample_id,'GR':mu,'mgDW':dry_weight,'Cells':cell_count},ignore_index=True)
    return data_out

###################################################
#### Simulate growth with variable constraints ####
###################################################
def simulate(model,light,photon_const,Po_const,YII_const,D1_const,DM20):
    a_star = pd.read_csv('pt_a_star.csv',index_col=0,header=0)
    a_star = a_star[light]
    # read in the cool white lamp spectral density
    rel = pd.read_csv('fluor_lamp.csv',header=0)
    volume = 375.
    if light == 'LL':
        mgCell = 19.1/1e9 # mgDW per cell <-- LL +/- 2.0
        innoc = 2.8e6*volume # cells/mL * total mL
        duration = 1440 # total length of the simulation
        
    if light == 'HL':
        mgCell = 20.4/1e9 # mgDW per cell <-- HL Bounds 20.4+/- 1.4
        innoc = 3.5e6*volume # cells/mL * total mL
        duration = 720 # total length of the simulation
    
    iDW = mgCell*innoc # initial culture biomass

    ### Cell count is mgDW/mgCell
    xsec = 80. #culture cross sectional area in cm2
    path_len = 4.7 #cm
    time_interval = 20
    gDW = iDW
    cells = (gDW/mgCell)
    #__________________________Initialize variables to collect output data___________________________________
    # Initialize an unused photon count
    biomass = pd.DataFrame(data=[gDW,cells],index = ['Biomass','Cells'],columns=['0']) # add it up at the end of each simulation
    # Initialize an O2 evolution rate tracker

    for t in range(time_interval,duration+time_interval,time_interval): #One simulation every 20 minutes
        import math as m
        interval_bm = 0 #initializes the total biomass
        photon_count = 0
        atten = np.zeros(len(a_star)) #captures light attenuation

        tmp_o2_evo = 0
        
        gDW = biomass[str(t-time_interval)]['Biomass'] #mg DW
        cells = (gDW/mgCell)

        if light == 'LL':
            #O2 rate equation  
            # input: umol photons cell-1 s-1
            # Output: umol O2 cell-1 s-1
            # Mean
            Ps = 2.94e-11 # umol O2 cell-1 s-1
            alpha = 9.82e-2
            beta = 0.0
            resp = 1.83e-12
            #     # upper bound    
            #     Ps = 3.06e-11 # umol O2 cell-1 s-1
            #     alpha = 9.75e-2
            #     beta = 0.0
            #     resp = 1.75e-12
            #     # lower bound
            #     Ps = 2.83e-11 # umol O2 cell-1 s-1
            #     alpha = 9.91e-2
            #     beta = 0.0
            #     resp = 1.92e-12
            irrad = 60.  
               
        if light == 'HL':
            #     # HL_bicarb
            Ps = (2.46e-11)*1.15 # umol O2 cell-1 s-1
            alpha = (9.47e-2)
            beta = 0.0
            resp = 4.72e-12
            #    #upper bound
                # HL
            #     Ps = (2.54e-11) # umol O2 cell-1 s-1
            #     alpha = (1.05e-1)
            #     beta = 3.67e-5
            #     resp = 4.28e-12    
            #   #lower bound
                ## HL
            #     Ps = (2.35e-11) # umol O2 cell-1 s-1
            #     alpha = (8.5e-2)#
            #     beta = 2.5e-10
            #     resp = 5.15e-12     
            irrad =  600. #uE 


        # Calculate photons and O2 uptake
        # Photon_flux is the initial amount of light delivered to the culture at each wavelength from 400-700 nm
        photon_flux = (rel['rel_height'].values)*irrad*xsec/10000*time_interval*60. #umol/(m2*s) * cm2 * 1m2/10000cm2 * 60s/min * min = umol photons/time interval
        total_photons = sum(photon_flux)
        
       # For each slice:
        #     - Walk through the spectrum calculating absorbed photons and updating the initial value (delivered)
        #     - At the end of the slice sum the totals to get photons absorbed
        #     - Update a running tally of photons absorbed
        total_absorbed = 0.

        ti_o2 = 0.
        for nm in range(len(photon_flux)):
            
            abs_coeff = a_star[400+nm]              # a* value for the given nm (cm2/cell)
            Io = photon_flux[nm]    # incident photon flux at this nm at this slice (umol_photon/time_interval)
            Ia = Io-Io*(m.exp(-1*abs_coeff*cells/volume*path_len))
            nm_abs = Ia
            
            total_absorbed = total_absorbed+nm_abs
        
        conv_abs = total_absorbed/time_interval/60./(cells) # converts abs to O2 evo curve units  umol/TI * TI/min * min/s * 1/cells => umol/(mgchla*s)
        slice_o2 = (Ps*(1-m.exp(-1*alpha*conv_abs/Ps))*(m.exp(-1*beta*conv_abs/Ps)))-resp #umol O2 cell-1 s-1
        ti_o2 = ti_o2+(slice_o2*(cells)*60.*time_interval) #umol O2 cell-1 s-1 * cells * s/min * min/TI = umol O2/TI

        if Po_const == True:
            # Constrained
            o2evo = ti_o2
            
            model.reactions.EX_o2_e.bounds = (0.99*o2evo,o2evo)

        else:
            #unconstrained
            model.reactions.EX_o2_e.bounds = (-1000.,1000.)   

            # ________________________________Time point constraints_________________________________________#

        ngam = resp*60.*time_interval*cells
        model.reactions.NGAM.lower_bound = ngam
        
        if light == 'LL':
            cef_val = 1.1 #5% of max LET
        if light == 'HL': 
            cef_val = 0.7 #5% of max LET
        
        cef_rate = cef_val*gDW/60.*time_interval #CEF sets CEF_h upper bound.  umol /(mgDW*h) * mgDW * h/60min * min/TI = umol/TI
        model.reactions.CEF_h.upper_bound = cef_rate
        
        
            #Photon constraints
            
        model.reactions.EX_photon_e.bounds = (total_absorbed*-1.,total_absorbed*-0.999)
            
        if YII_const == True:
            if light == 'LL':
            ###### Parameters for PSII fluorescence
            ## Fv/Fm
                FvFm_LL = 0.69  
            ## Calculate Y(II) based on absorbed 
                abs_conv = total_absorbed/time_interval/60./(cells)
                yII_LL = 0.7016*np.exp(-8.535e8*abs_conv)
            ## Y(NPQ)
                yNPQ_LL = 1./(1.+(np.power(2.3e-9/abs_conv,3.)))
                if yNPQ_LL < 0:
                    yNPQ_LL = 0.     
                phoYII = round(yII_LL,2) # Y(II)        
                regNPQ = round(yNPQ_LL,2) # Y(NPQ)
                unrNPQ = round((1-phoYII-regNPQ),2) # Y(NO)

            if light == 'HL':
            ###### Parameters for PSII fluorescence
            ## Fv/Fm
                FvFm_HL = 0.63
            ## Calculate Y(II) based on absorbed 
                abs_conv = total_absorbed/time_interval/60./(cells)
                yII_HL = 0.6398*np.exp(-1.169e9*abs_conv)
            ## Y(NPQ)
                yNPQ_HL = 1./(1.+(np.power(1.37e-9/abs_conv,5.5)))
                if yNPQ_HL < 0:
                    yNPQ_HL = 0.                                            
                phoYII = round(yII_HL,2) # Y(II)
                regNPQ = round(yNPQ_HL,2) # Y(NPQ)
                unrNPQ = round((1-phoYII-regNPQ),2) # Y(NO)
        ### Edit model constraints with the appropriate values
        ## PHO_PSIIt_u
            rxn = model.reactions.PHO_PSIIt_u
            
            ## reset the stoich
            for met,s in rxn.metabolites.items():
                stoich = s*-1
                temp_dict = {met:stoich}
                rxn.add_metabolites(temp_dict)
                
            m1 = model.metabolites.photon_YII_u
            m2 = model.metabolites.photon_YNPQ_u
            m4 = model.metabolites.photon_YNO_u
            m3 = model.metabolites.photon_h
            
            rxn.add_metabolites({m1:phoYII,
                                m2:regNPQ,
                                m4:unrNPQ,
                                m3:-1.})
        else:
        ### Edit model constraints with the appropriate values
        ## PHO_PSIIt_u
            rxn = model.reactions.PHO_PSIIt_u
            
            ## reset the stoich
            for met,s in rxn.metabolites.items():
                stoich = s*-1
                temp_dict = {met:stoich}
                rxn.add_metabolites(temp_dict)
                
            m1 = model.metabolites.photon_YII_u
            m2 = model.metabolites.photon_YNPQ_u
            m3 = model.metabolites.photon_h

        # No Y(II) fraction at PSII        
            rxn.add_metabolites({m1:(1.),
                                 m2:(0.),
                                 m3:-1.})   
        
    ## Photon constraint
        if DM20 == True:
            model.reactions.DM_photon_c.bounds = (-0.9999*0.2*(model.reactions.EX_photon_e.lower_bound),-0.2*(model.reactions.EX_photon_e.lower_bound))# 20% assumption 

        elif photon_const == True:
            model.reactions.DM_photon_c.bounds = (0.,0.) # constrained
            
        else:
            model.reactions.DM_photon_c.bounds = (0.,10000.) # Unconstrained

        if D1_const == True:
            # Add D1 damage cost uncoupled to PSII
            # Damage rates
            if light == 'HL':
                D1_rate = 2.52e-4# # HL: umol D1/mgDW h-1         <---- D1 Constraint
            if light == 'LL':
                D1_rate = 7e-6# # LL: umol D1/mgDW h-1         <---- D1 Constraint
            D1_rate = D1_rate * gDW/60.# ugD1 ugDW-1 min-1 * mgDW * 1h/60 min = umolD1 min-1
            D1_rate = D1_rate * time_interval # umolD1 min-1 * min/TI = umol D1/TI

            
            model.reactions.NGAM_D1_u.bounds = (D1_rate, 1.0001*(D1_rate))
        else:
            model.reactions.NGAM_D1_u.bounds = (0.,0.) 
            
        ## Solve the model
        model.objective = 'bof_c'
        solution = cobra.flux_analysis.parsimonious.optimize_minimal_flux(model) 
        if solution.status == 'optimal':
            obj_rxn = model.reactions.bof_c
            
            biomass[str(t)]=(gDW+obj_rxn.x,(gDW+obj_rxn.x)/(mgCell))
            
    mu = np.log(biomass[str(duration)]['Cells']/biomass['0']['Cells'])/(duration/60)
    save_model = model
    print('***Simulation values (mmol gDW-1 h-1)***\n',
          '**Values are for the final simulation interval**\n',
          '*Growth rate is for the entire culture duration*\n',
          'Total photons absorbed: ', (np.around(model.reactions.EX_photon_e.lower_bound/gDW*-3.,1)),'\n',
          'Photons allowed to be lost upstream of the photosystems: ', (np.around(model.reactions.DM_photon_c.upper_bound/gDW*3.,1)),'\n',
          'Oxygen exchange constraint: ', (np.around(model.reactions.EX_o2_e.upper_bound/gDW*3.,1)),'\n',
          'Y(II) constraint: ', model.reactions.PHO_PSIIt_u.reaction,'\n',
          'D1 damage constraint: ',(np.around(model.reactions.NGAM_D1_u.upper_bound/gDW*3.,1)),'\n',
          'Growth rate: ',np.around(mu,3),' h-1\n')
    return biomass, save_model



















