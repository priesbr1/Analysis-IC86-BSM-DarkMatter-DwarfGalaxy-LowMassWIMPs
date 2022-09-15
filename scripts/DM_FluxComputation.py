#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py3-v4.1.0/icetray-start
#METAPROJECT /data/user/niovine/software/combo_py3/build/
import os
import math
import pickle
import numpy as np
import scipy
from scipy import interpolate

#---------------------------------------------------------------------
# NOTE
# The computation of the J-factor can be find in 
# /home/niovine/projects/DarkMatter_OscNext/PythonCode/DarkMatter_Signal/Jfactor/
# Same goes for the spectra which was extracted from PPPC4/charon tables 
# and oscillated with scripts located in
# /home/niovine/projects/DarkMatter_OscNext/PythonCode/DarkMatter_Signal/Spectra/
#---------------------------------------------------------------------


#---------------------------------------------------------------------
# J-factor: dJ/dOmega(psi) [GeV^2 cm^{-5} sr^{-1}]
# NOTE: Psi has to be in degree
#---------------------------------------------------------------------
def interpolate_Jfactor(psi_value, profile):
        
    # Open file
    Jfactor = pickle.load(open("/data/user/niovine/projects/DarkMatter_OscNext/SignalWeight/Jfactor/Clumpy_Jfactor_{}_GeV2_cm5sr1.pkl".format(profile),"rb"))
    
    y_interp = scipy.interpolate.splrep(Jfactor["psi"], Jfactor["Jpsi"])
    interp_Jpsi = scipy.interpolate.splev(psi_value, y_interp, der=0)
    
    return interp_Jpsi



#---------------------------------------------------------------------
# Oscillate spectra
#---------------------------------------------------------------------
# Computing the PNMS matrices Uij for the oscillation matrix
def PMNS_matrix(t12, t13, t23, d):
    s12 = np.sin(t12)
    c12 = np.cos(t12)
    s23 = np.sin(t23)
    c23 = np.cos(t23)
    s13 = np.sin(t13)
    c13 = np.cos(t13)
    cp  = np.exp(1j*d)
    
    return np.array([[ c12*c13, s12*c13, s13*np.conj(cp)],
                  [-s12*c23 - c12*s23*s13*cp, c12*c23 - s12*s23*s13*cp, s23*c13],
                  [ s12*s23 - c12*s23*s13*cp,-c12*s23 - s12*c23*s13*cp, c23*c13]])

# Probability of flavor to change when L->inf
def prob(a, b, U):
    # Gives the oscillation probability for nu(a) -> nu(b)
    # for PMNS matrix U, and L in km and E in GeV
    s = 0
    for i in range(3):
            s += (np.conj(U[a,i])*U[b,i]*U[a,i]*np.conj(U[b,i])).real
    return s

# Define Oscillation Matrix
def osc_matrix(U):
    return np.array([[prob(0, 0, U), prob(0, 1, U), prob(0,2,U)],
                     [prob(1, 0, U), prob(1, 1, U), prob(1,2,U)],
                     [prob(2, 0, U), prob(2, 1, U), prob(2,2,U)]])

# Oscillate the spectra
def oscillate_spectra(spectra, nutypes=["nu_mu", "nu_e", "nu_tau"]):
    
    oscillated = dict()  
    
    # Define mixing angles theta
    t12 = 0.57596 # Old value: 0.5934
    t13 = 0.1296  # Old value: 0.1514
    t23 = 0.8552  # Old value: 0.785398
    U = PMNS_matrix(t12, t13, t23, 0)
    P = osc_matrix(U)
    
    #print("Oscillation parameters: theta12={}, theta13={}, theta23={}".format(str(t12), str(t13), str(t23)))
    #print("Oscillation matrix: ", P)

    for var in spectra.keys():

        if "dNdE" in var:
            dNdE_nue_osc = []
            dNdE_numu_osc = []
            dNdE_nutau_osc = []

            #Apply Oscillation Matrix
            for i in range(len(spectra[var])):
                dNdE_nue_osc.append(np.dot(P,np.array([spectra[var][i], 
                                                       spectra[var][i], 
                                                       spectra[var][i]]))[0])
                dNdE_numu_osc.append(np.dot(P,np.array([spectra[var][i], 
                                                        spectra[var][i], 
                                                        spectra[var][i]]))[1])
                dNdE_nutau_osc.append(np.dot(P,np.array([spectra[var][i], 
                                                         spectra[var][i], 
                                                         spectra[var][i]]))[2])

                # Spectra after oscillation
                oscillated[var+"_"+nutypes[0]+"_osc"] = dNdE_nue_osc
                oscillated[var+"_"+nutypes[1]+"_osc"] = dNdE_numu_osc
                oscillated[var+"_"+nutypes[2]+"_osc"] = dNdE_nutau_osc

        else:
            oscillated[var] = np.array(spectra[var])
    
    return oscillated



#---------------------------------------------------------------------
# Spectra interpolation (either Charon or PPPC4 tables)
#---------------------------------------------------------------------
def interpolate_spectra(spectra_type, E, pdg, channel, mass):
    
    path = "/data/user/niovine/projects/DarkMatter_OscNext/SignalWeight/Spectra/{}/".format(spectra_type)
    spectra_file = pickle.load(open(path+"Spectra_ann_{}_atEarth_OscillationParametersFromArXiv1811_05487.pkl".format(spectra_type),"rb"))
    
    # Different treatment of nu and anti-nu spectra
    if spectra_type == "PPPC4":
        # No distinction is made between neutrinos and anti-neutrinos in PPPC4 tables
        nu_types = ["nu_e", "nu_mu", "nu_tau"]
        pdg_encoding = {"nu_e":12, "nu_mu":14, "nu_tau":16}
        pdg = np.abs(pdg)
    elif spectra_type == "Charon":
        nu_types = ["nu_e", "nu_mu", "nu_tau", "nu_e_bar", "nu_mu_bar", "nu_tau_bar"]
        pdg_encoding = {"nu_e":12, "nu_mu":14, "nu_tau":16, "nu_e_bar":-12, "nu_mu_bar":-14, "nu_tau_bar":-16}
        
    
    # Define array holding interpolated values of spectra
    interp_dNdE = np.zeros(len(pdg))
    
    for nu_type in nu_types:
        
        # Energy
        spectra_E = np.array(spectra_file[channel][str(mass)][nu_type]["E"])
        # Spectra
        spectra_dNdE = np.array(spectra_file[channel][str(mass)][nu_type]["dNdE"])
         
        # print("Energy to interpolate:", spectra_E)
        # print("Spectra to interpolate:", min(spectra_dNdE), max(spectra_dNdE))
        
        # Define low energy and high energy cut for interpolation
        # Check if there is zero in spectra (for nunu) and define interpolation cuts accordingly
        HE_cut = mass
        zeroes = np.where(spectra_dNdE == 0.)
        if zeroes != np.array([]):
            LE_cut = min(spectra_E[np.where(spectra_dNdE!=0.)])+(0.01*HE_cut)
        else:
            LE_cut = min(spectra_E)
        
        # print( "Energy range for interpolation:[{},{}]".format(str(LE_cut),str(HE_cut)) )
        # print( "Real interpolation range: [{},{}]".format(str(min(spectra_E)),str(max(spectra_E))) )
        
        # Define fct from which interpolate from
        y_interp = scipy.interpolate.splrep(spectra_E, spectra_dNdE)
        
        # Actually interpolate the spectra for our energy array
        # Only interpolate for the proper nu_type & above the energy threshold of the spectra
        interp_dNdE = np.where((pdg==pdg_encoding[nu_type]) & (E>=LE_cut) & (E<=HE_cut), 
                               scipy.interpolate.splev(E, y_interp, der=0), interp_dNdE)
    
    
    # Get rid of negative spectra values due to poor interpolation
    while len(E[np.where(interp_dNdE<0)]) != 0:
        neg_v = np.where(interp_dNdE<0)[0]
        interp_dNdE[neg_v] = interp_dNdE[neg_v-1]
    
    return interp_dNdE, E



#---------------------------------------------------------------------
# Define cut on weight
#---------------------------------------------------------------------
def define_weightcut(weight, cut):
    
    H, edges = np.histogram(weight, bins=1000)
    zeroes = np.where(H==0.)
    
    # print(zeroes)
    
    i = 0
    n = 0

    while (i<zeroes[0].shape[0]-1) and (n<cut+1):
        
        # Check if consecutive zeroes
        if zeroes[0][i]+1 == zeroes[0][i+1]:
            n += 1
        # Reset n to zero if encounter non-null value
        elif zeroes[0][i]+1 != zeroes[0][i+1]:
            # print("Reset to zero")
            n = 0

        if n == cut:
            loc = zeroes[0][i]
            # print("Location:", loc)
    
        i+=1
        
    if n >= cut:
        w_lim = edges[loc]
    else:
        w_lim = max(weight)
        
    return w_lim
