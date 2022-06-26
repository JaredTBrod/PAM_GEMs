import cobra
from cobra.core.metabolite import elements_and_molecular_weights
elements_and_molecular_weights['R']=0.0
elements_and_molecular_weights['Z']=0.0
import pandas as pd
import numpy as np
import csv
from cobra.flux_analysis.parsimonious import pfba
###################################################
#### Simulate growth with variable constraints ####
###################################################
def simulate(model,light,photon_const,Po_const,YII_const,D1_const,DM20,ngam_comp):
    a_star = pd.read_csv('pt_a_star.csv',index_col=0,header=0)
    a_star = a_star[light]
    # read in the cool white lamp spectral density
    rel = pd.read_csv('fluor_lamp.csv',header=0)
    volume = 375.
    if light == 'LL':
        mgCell = 19.1/1e9 # mgDW per cell <-- LL +/- 2.0
        innoc = 2.8e6*volume # cells/mL * total mL
        duration = 1440 # total length of the simulation
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
        mgCell = 20.4/1e9 # mgDW per cell <-- HL Bounds 20.4+/- 1.4
        innoc = 3.5e6*volume # cells/mL * total mL
        duration = 720 # total length of the simulation
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
            model.reactions.EX_o2_e.bounds = (-10000.,10000.)

            # ________________________________Time point constraints_________________________________________#

        ngam = resp*60.*time_interval*cells
        if ngam_comp == 'mito':
            model.reactions.NGAM.lower_bound = ngam
            model.reactions.PTOX_h.bounds = (0.,0.)
            
        if ngam_comp == 'chloro':
            model.reactions.PTOX_h.bounds = (ngam*0.9999,ngam)
            
            model.reactions.NGAM.bounds = (0.,0.)
        
        if light == 'LL':
            cef_val = 1.1
        if light == 'HL':
            cef_val = 0.7
        
        cef_rate = cef_val*gDW/60.*time_interval #CEF sets CEF_h upper bound.  umol /(mgDW*h) * mgDW * h/60min * min/TI = umol/TI
        model.reactions.CEF_h.upper_bound = cef_rate
        
        
            #Photon constraints
        model.reactions.EX_photon_e.bounds = (total_absorbed*-1.,total_absorbed*-0.9999)
            
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
            model.reactions.DM_photon_c.bounds = (-0.9999*0.2*(model.reactions.EX_photon_e.lower_bound),-0.2*(model.reactions.EX_photon_e.lower_bound)) # 20% assumption 
    
        elif photon_const == True:
            model.reactions.DM_photon_c.bounds = (0.,0.) # constrained


        else:
            model.reactions.DM_photon_c.bounds = (0.,100000.) # Unconstrained

        if D1_const == True:
            # Add D1 damage cost uncoupled to PSII
            # Damage rates

            if light == 'HL':
                D1_rate = 2.52e-4# # HL: umol D1/mgDW h-1         <---- D1 Constraint
            if light == 'LL':
                D1_rate = 7e-6# # LL: umol D1/mgDW h-1         <---- D1 Constraint
            D1_rate = D1_rate * gDW/60.# ugD1 ugDW-1 min-1 * mgDW * 1h/60 min = umolD1 min-1
            D1_rate = D1_rate * time_interval # umolD1 min-1 * min/TI = umol D1/TI
            model.reactions.NGAM_D1_u.bounds = (D1_rate,1.0001*(D1_rate)) #
        else:
            model.reactions.NGAM_D1_u.bounds = (0.,0.) #

        ## Solve the model
        model.objective = 'bof_c'
        solution = cobra.flux_analysis.parsimonious.optimize_minimal_flux(model) 
        if solution.status == 'optimal':
            obj_rxn = model.reactions.bof_c
            
            biomass[str(t)]=(gDW+obj_rxn.x,(gDW+obj_rxn.x)/(mgCell))
    mu = np.log(biomass[str(duration)]['Cells']/biomass['0']['Cells'])/(duration/60)
    print('Growth rate: ',np.around(mu,3))
    out_model = model.copy()
    return biomass, out_model


#A function for calculating excess electron transport (EET)
####################################
## Inputs: 
### Model: the GEM
### ref_gDW: the biomass in the culture at the time of EET determination
### TI: the time interval for the simulation 
#####
## Returns:
### Linear electron transport, Excess electron transport

def calc_eet(model, ref_gDW, TI):# 
    from cobra.flux_analysis.parsimonious import pfba
    #1 maximize growth
    model2 = model.copy()
    model2.reactions.EX_o2_e.lower_bound = -1000.
    max_mu = model2.slim_optimize() ## maximize the biomass production
    
    #2 limit mu max
    model2.reactions.bof_c.upper_bound = max_mu 
    model2.reactions.bof_c.lower_bound = max_mu*0.9999
    sol = pfba(model2) # Get the parsimonious solution
    
    # Calculate Linear Electron Transport and Excess Electron Transport
    bof_let1 = (sol.fluxes['PSI_u'] - sol.fluxes['CEF_h'])/ref_gDW/TI
        
    #3 determine LET to generate biomass and maintenance only
    #3a: Open photon demand and maximize so that the model uses the minimum LET
    model2.reactions.DM_photon_c.upper_bound = 1e5
    model2.objective = 'DM_photon_c'
    #3b: get the parsimonious solution to this
    sol = pfba(model2)
    bof_let2 = (sol.fluxes['PSI_u'] - sol.fluxes['CEF_h'])/ref_gDW/TI
    
    ## EET is the excess electrons, in LET, above biomass and mx requirements 
    
    return bof_let1,bof_let2

# A function for constraining and retrieving the EET for HL, LL or DM20%
#########################################

def sim_eet(base_model, light_level, DM_flag, ngam_comp):
    from cobra.flux_analysis.parsimonious import pfba
    import warnings
    import math as m
    
    # get solution for given light level
    if DM_flag == True:
        YII_flag = False
    else:
        YII_flag = True
        
    (biomass_out, model_out) = simulate(base_model,light=light_level,photon_const=True,
                       Po_const=True,YII_const=YII_flag,D1_const=True,DM20=DM_flag, ngam_comp = ngam_comp)

    ## This will return a model and the biomass DF for the simulation at either HL, LL, and +/- DM20% assumption

    
    ## Function for calulating AEF
    model = model_out

    if light_level == 'HL':
        end_point = '700'
        
    elif light_level == 'LL':
        end_point = '1420'

    
    ref_gDW = biomass_out[end_point]['Biomass'] # mg DW at the time of the simulation
    cells = biomass_out[end_point]['Cells'] ## number of cells at the time of the sim
    TI = 20./60. # simulation time interval in hours
    ref_qf = model.reactions.EX_photon_e.lower_bound # model QF at the time of the sim
    ## Run it at fractions of current QF
    QF_range = np.linspace(1e-10,1.35e-9,21)/(ref_gDW/cells)*3600*TI*ref_gDW
    frac = -1*QF_range/ref_qf

    QF_out = []
    YII_eetLET_out = []
    LET_out = []
    for f in frac:
        model = model_out.copy()
        qf = model.reactions.EX_photon_e.lower_bound
        model.reactions.EX_photon_e.lower_bound = (qf * f, qf * f * 0.999)
        abs_conv = -1.*qf*f/20./60./(cells)
        ## recalc O2
        if light_level == 'HL':
            Ps = (2.46e-11)*1.15 # umol O2 cell-1 s-1
            alpha = (9.47e-2)
            beta = 0.0
            resp = 4.72e-12
            
        if light_level == 'LL':
            Ps = 2.94e-11 # umol O2 cell-1 s-1
            alpha = 9.82e-2
            beta = 0.0
            resp = 1.83e-12
            
        ti_o2 = 0.
        slice_o2 = (Ps*(1-m.exp(-1*alpha*abs_conv/Ps))*(m.exp(-1*beta*abs_conv/Ps)))-resp #umol O2 cell-1 s-1
        ti_o2 = ti_o2+(slice_o2*(cells)*60.*20.) #umol O2 cell-1 s-1 * cells * s/min * min/TI = umol O2/TI
        o2evo = ti_o2
        model.reactions.EX_o2_e.bounds = (0.,o2evo) # o2evo
       # for Y(II), need to recalc Y(II) based on QF

        if YII_flag == True:

        ## Fv/Fm
            if light_level == 'HL':
                FvFm_HL = 0.63
                yII_HL = 0.6398*np.exp(-1.169e9*abs_conv)
                if abs_conv > 0:
                    yNPQ_HL = 1./(1.+(np.power(1.37e-9/abs_conv,5.5)))
                else:
                    yNPQ_HL = 0.                                           
            ## Constraints 
                phoYII = round(yII_HL,2) # Y(II)
                regNPQ = round(yNPQ_HL,2) # Y(NPQ)
                unrNPQ = round((1-phoYII-regNPQ),2) # Y(NO)

            if light_level == 'LL':
                FvFm_LL = 0.68
                yII_LL = 0.7016*np.exp(-8.535e8*abs_conv)
                if abs_conv > 0:
                    yNPQ_LL = 1./(1.+(np.power(2.3e-9/abs_conv,3.)))
                else:
                    yNPQ_LL = 0.                                           
            ## Constraints 
                phoYII = round(yII_LL,2) # Y(II)
                regNPQ = round(yNPQ_LL,2) # Y(NPQ)
                unrNPQ = round((1-phoYII-regNPQ),2) # Y(NO)
                
        # Set Y(II) constraints   
            rxn = model.reactions.PHO_PSIIt_u
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

        if DM_flag == True:  #DM20 assumption, no Y(II), no D1 damage, 20% of photons can leave upstream of the photosystems
            rxn = model.reactions.PHO_PSIIt_u
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


            model.reactions.DM_photon_c.bounds = (-0.9999*0.2*(model.reactions.EX_photon_e.lower_bound), -0.2*(model.reactions.EX_photon_e.lower_bound))# 20% assumption 
            
            model.reactions.NGAM_D1_u.bounds = (0.,0.) # no D1 NGAM for the DM20 assumption

        (a,b) = calc_eet(model, ref_gDW, TI)
        YII_eetLET_out.append(a-b)
        QF_out.append(f*ref_qf/ref_gDW/TI)
        LET_out.append(a)

    return YII_eetLET_out, QF_out, LET_out 


### A function for calculating cross-compartment fluxes for a given solution
def calc_fluxes(model,xcomp_rxns,gDW,TI):
    from cobra.flux_analysis.parsimonious import pfba
    import numpy as np
    rxn4e = ['PTOX_h','NGAM','AOX_m'] # reactions with 4 e- per reaction, all others are 2 e-
    rev_rxns = ['GAPDH_m','SDH_m','MDH_m']
    rxns_output = dict()
    soln = pfba(model)
    total_red = 0
    for k,v in xcomp_rxns.items():
        total_flux = 0
        rxn_list = v
        for rx in v:
            rx_flux = soln.fluxes[rx]        
            if rx in rxn4e:
                rx_flux = 4*rx_flux
            elif rx == 'GCS_m':
                rx_flux = rx_flux
           
            elif rx in rev_rxns:
                rx_flux = rx_flux*-2
                
            else:
                rx_flux = 2*rx_flux

            total_flux = total_flux + rx_flux/gDW/TI*60
            
        if total_flux < 0:
            total_flux = total_flux*-1
        rxns_output[k]=np.around(total_flux,2)
    
    for_total = ['Glyco_m','ORN','PRred','BCAA','TCA','Lysine']
    for k,v in rxns_output.items():
        if k in for_total:
            total_red = total_red+v   
    rxns_output['Total']=np.around(total_red,2)
    return rxns_output

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

    print(model.reactions.get_by_id(rxn).reaction)
    print('')
    print(model.metabolites.get_by_id(met.id).formula)
    print('')
    print(model.metabolites.get_by_id(met.id).formula_weight)
    print('')
    print(model.reactions.get_by_id(rxn).check_mass_balance())
    return model


def get_comp_let(model,DW_accum):
    out_var = {}
    labels = {'biomass_carb_c':'Structural carb',
          'biomass_dna_c': 'DNA',
          'biomass_mem_lipids_c': 'Membrane lipids',
          'biomass_pigm_h':'Pigments',
          'biomass_plastid_lip_h':'Plastid lipids',
          'biomass_pro_c':'Protein',
          'biomass_rna_c':'RNA',
          'carbon_storage_c':'Storage'}
    
    ## Cycle through the components
    for k,v in labels.items():
        met = k
        label = v

        ## get the current bof
        temp_model = model.copy()
        bof_ratio = {}
        for m in temp_model.reactions.bof_c.metabolites:
            if temp_model.reactions.bof_c.metabolites[m] < 0.:
                bof_ratio[m.id]=np.around(-1*temp_model.reactions.bof_c.metabolites[m],3) 

        ## Reduce the maximum biomass production by the component percentage
        max_bof = (1-bof_ratio[met])*DW_accum
        temp_model.reactions.bof_c.upper_bound = max_bof
        temp_model.reactions.bof_c.lower_bound = max_bof*0.999
        temp_model.reactions.EX_o2_e.lower_bound = 0.
        ## remove the selected macromolecule
        bof_ratio[met]=0.

        ## update the BOF
        bof_data = dict()

        for bofcmp,percent in bof_ratio.items():
            met_obj = model.metabolites.get_by_id(bofcmp)
            bof_data[met_obj]=percent*-1

        temp_model = update_biomass(temp_model,'bof_c',bof_ratio,'biomass_c')

        ## simulate
        temp_model.reactions.DM_photon_c.upper_bound = 1e5
        temp_model.objective = 'DM_photon_c'
        soln = cobra.flux_analysis.parsimonious.pfba(temp_model)

        ## 
        LET = np.around(2*(soln['PSI_u']-soln['CEF_h']),1)

        out_var[v]=LET
        
    return out_var


def target_met(model,met):
    reaction = model.reactions.DM_biomass_c.copy()
    reaction.id = 'DM_'+str(met)
    met_obj = model.metabolites.get_by_id(met)
    reaction.name = 'Demand Reaction: '+met_obj.name
    m1 = model.metabolites.biomass_c
    m2 = met_obj

    reaction.add_metabolites({m1:1,
                             m2:-1})
    reaction.upper_bound = 0.
    reaction.lower_bound = 0.

    model.add_reaction(reaction)
    model.repair()
    return model



def calc_eet(model, frac_max, met, ref_gDW, TI, get_fluxes, get_model):# 
    from cobra.flux_analysis.parsimonious import pfba
    #1 maximize growth
    model2 = model.copy()
    model2.reactions.EX_o2_e.lower_bound = model2.reactions.EX_o2_e.upper_bound*0.999
    max_mu = model2.slim_optimize()
    
    #2 limit mu to the desired value
    model2.reactions.bof_c.upper_bound = max_mu*frac_max
    model2.reactions.bof_c.lower_bound = max_mu*frac_max*0.99999
    sol = pfba(model2)
    bof_let1 = (sol.fluxes['PSI_u'] - sol.fluxes['CEF_h'])/ref_gDW/TI*2.
    bof_eet1 = (sol.fluxes['O2t_m'])/ref_gDW/TI*4.
#     print("LET baseline: ",np.around(bof_let1,3), "mmol e- gDW-1 h-1")
#     ref_rubiso = sol.fluxes['RUBISO_h']
    
    
    #3 determine LET to achieve this growth rate
    #3a: open photon demand and maximize so that the model uses the minimum LEF
    model2.reactions.DM_photon_c.upper_bound = 1e5
    model2.objective = 'DM_photon_c'
#     model2.reactions.RUBISO_h.lower_bound = ref_rubiso
    #3b: get the parsimonious solution to this
    sol = pfba(model2)
    bof_let2 = (sol.fluxes['PSI_u'] - sol.fluxes['CEF_h'])/ref_gDW/TI*2.
    bof_eet2 = (sol.fluxes['O2t_m'])/ref_gDW/TI*4.
#     print("LET for biomass: ",np.around(bof_let2,3), "mmol e- gDW-1 h-1")
#     print("Baseline EET: ", np.around(bof_let1-bof_let2,3), "mmol e- gDW-1 h-1")
#     print("")
#     4 Open the exchange reaction for the metabolite of interest
    if met:
        met_demand = 'DM_'+met
        model2.reactions.get_by_id(met_demand).upper_bound = 1e5
        model2.reactions.DM_photon_c.upper_bound = 0.
        model2.objective = met_demand
#         model2.reactions.RUBISO_h.lower_bound = ref_rubiso
        sol = pfba(model2)
        bof_let3 = (sol.fluxes['PSI_u'] - sol.fluxes['CEF_h'])/ref_gDW/TI*2.
        bof_eet3 = (sol.fluxes['O2t_m'])/ref_gDW/TI*4.
#         print("Target met yield:", np.around(sol.fluxes[met_demand]/ref_gDW/TI,3),"mmol gDW-1 h-1")
#         print("LET with target met production: ",np.around(bof_let3,3), "mmol e- gDW-1 h-1")
#         print("")
        if get_fluxes == True:
            fluxes = sol.fluxes
            if get_model == True:
                return (model2,fluxes)
            else:
                return fluxes
        
        #5 new EET with production
        model2.reactions.DM_photon_c.upper_bound = 1e5
#         model2.reactions.RUBISO_h.lower_bound = ref_rubiso
        model2.reactions.get_by_id(met_demand).lower_bound = sol.fluxes[met_demand]*0.999
        model2.objective = 'DM_photon_c'
        sol = pfba(model2)
        bof_let4 = (sol.fluxes['PSI_u'] - sol.fluxes['CEF_h'])/ref_gDW/TI*2.
        bof_eet4 = (sol.fluxes['O2t_m'])/ref_gDW/TI*4.
#         print("LET for met production: ",np.around(bof_let4,3), "mmol e- gDW-1 h-1")
#         print("Met production EET: ", np.around(bof_eet1-bof_eet4,3), "mmol e- gDW-1 h-1")
#         print("*****",frac_max)



def get_prod_env(model,met_eng,bof_frac,inoc,light,time_interval):
    a_star = pd.read_csv('pt_a_star.csv',index_col=0,header=0)
    a_star = a_star[light]
    # read in the cool white lamp spectral density
    rel = pd.read_csv('fluor_lamp.csv',header=0)
    volume = 375.
    innoc = inoc*volume
    met_demand = 'DM_'+met_eng
    if light == 'LL':
        mgCell = 19.1/1e9 # mgDW per cell <-- LL +/- 2.0
        duration = 7200 # total length of the simulation
        #O2 rate equation  
        # input: umol photons cell-1 s-1
        # Output: umol O2 cell-1 s-1
        # Mean
        Ps = 2.94e-11 # umol O2 cell-1 s-1
        alpha = 9.82e-2
        beta = 0.0
        resp = 1.83e-12
        irrad = 60.  

    if light == 'HL':
        mgCell = 20.4/1e9 # mgDW per cell <-- HL Bounds 20.4+/- 1.4
        duration = 7200 # total length of the simulation
    #     # HL_bicarb
        Ps = (2.46e-11)*1.15 # umol O2 cell-1 s-1
        alpha = (9.47e-2)
        beta = 0.0
        resp = 4.72e-12
        irrad =  600. #uE 

    iDW = mgCell*innoc # initial culture biomass

    ### Cell count is mgDW/mgCell
    xsec = 80. #culture cross sectional area in cm2
    path_len = 4.7 #cm
    time_interval = time_interval
    gDW = iDW
    cells = (gDW/mgCell)
    #__________________________Initialize variables to collect output data___________________________________
    # Initialize an unused photon count
    biomass = pd.DataFrame(data=[gDW,cells],index = ['Biomass','Cells'],columns=['0']) # add it up at the end of each simulation
    met_prod = pd.DataFrame(data=[0],index = ['Met Production'],columns=['0'])
    # Initialize an O2 evolution rate tracker
    for t in range(time_interval,duration+time_interval,time_interval): #One simulation every 20 minutes
        import math as m
        interval_bm = 0 #initializes the total biomass
        photon_count = 0
        atten = np.zeros(len(a_star)) #captures light attenuation

        tmp_o2_evo = 0

        gDW = biomass[str(t-time_interval)]['Biomass'] #mg DW
        met_val = met_prod[str(t-time_interval)]['Met Production']
        cells = (gDW/mgCell)

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

        o2evo = ti_o2
        model.reactions.EX_o2_e.bounds = (0.9*o2evo, o2evo)

            # ________________________________Time point constraints_________________________________________#

        ngam = resp*60.*time_interval*cells
        model.reactions.NGAM.lower_bound = ngam
        model.reactions.PTOX_h.bounds = (0.,0.)

        if light == 'LL':
            cef_val = 1.1
        if light == 'HL':
            cef_val = 0.7

        cef_rate = cef_val*gDW/60.*time_interval #CEF sets CEF_h upper bound.  umol /(mgDW*h) * mgDW * h/60min * min/TI = umol/TI
        model.reactions.CEF_h.upper_bound = cef_rate


            #Photon constraints
        model.reactions.EX_photon_e.bounds = (total_absorbed*-1.,total_absorbed*-0.9999)

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
    ## Photon constraint
        model.reactions.DM_photon_c.bounds = (0.,0.) # constrained

        # Add D1 damage cost uncoupled to PSII
        # Damage rates

        if light == 'HL':
            D1_rate = 2.52e-4# # HL: umol D1/mgDW h-1         <---- D1 Constraint
        if light == 'LL':
            D1_rate = 7e-6# # LL: umol D1/mgDW h-1         <---- D1 Constraint
        D1_rate = D1_rate * gDW/60.# ugD1 ugDW-1 min-1 * mgDW * 1h/60 min = umolD1 min-1
        D1_rate = D1_rate * time_interval # umolD1 min-1 * min/TI = umol D1/TI
        
        model.reactions.NGAM_D1_u.bounds = (D1_rate, 1.0001*(D1_rate))

        ## Solve the model
        model.objective = 'bof_c'
        solution = cobra.flux_analysis.parsimonious.optimize_minimal_flux(model) 
        if solution.status == 'optimal':
            obj_rxn = model.reactions.bof_c

            max_mu = obj_rxn.x*(1-bof_frac)
            model2 = model.copy()
            model2.reactions.bof_c.bounds = (max_mu*0.9999,max_mu)

            model2.objective = met_demand
            model2.reactions.get_by_id(met_demand).upper_bound = 1e5
            solution2 = model2.optimize() 
            if solution2.status == 'optimal':
                biomass[str(t)]=(gDW+max_mu,(gDW+max_mu)/(mgCell))
                met_prod[str(t)]=(met_val+solution2.fluxes[met_demand])
    return met_prod

def plot_prod(data_DF, time_interval,metgfw,bof_frac):
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import axes3d
    time_interval = time_interval

    data_DF = data_DF.sort_index()
    x = [float(c)/(60*24) for c in data_DF.columns]
    y = data_DF.index
    X,Y = np.meshgrid(x,y)
    Z = data_DF*metgfw/1000./(bof_frac*100.) # umol -> mg/%biomass rerouted

    fig = plt.figure(figsize=(10,18))
    ax = fig.add_subplot(2,1,1, projection='3d')

    # Plot a 3D surface
    mycmap = plt.get_cmap('inferno')

    surf1 = ax.plot_surface(X, Y, Z, linewidth = 0, cmap = mycmap,vmax=27.1)


    fig.colorbar(surf1, ax=ax, shrink = 0.5, aspect = 10)
    ax.set_zlim(0, 30)
    ax.yaxis.set_ticks(np.arange(0., 2.01e8, 0.4e8))
    ax.view_init(20,210)

    plt.show()

