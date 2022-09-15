#!/usr/bin/env python
import numpy as np
import glob
from I3Tray import *
from icecube import icetray, dataclasses, dataio
from math import *
from icecube.dataclasses import *
import argparse
from event_cuts_l5_plus import *
import time
import pickle

parser = argparse.ArgumentParser(description="Compute signal and bg weights for determining optimal energy range")
parser.add_argument("-d", "--datafolder", type=str, default="/mnt/research/IceCube/jpandre/Matt/level5p/",
                    dest="datafolder", help="folder containing i3 MC datafiles")
parser.add_argument("-s", "--spectra", type=str, default="/mnt/home/priesbr1/DM_Search/data/annihilation_spectra/DM_CirelliSpectrum_dictNoEW_neutrinos_oscillated.npy",
                    dest="spectra", help="file containing WIMP annihilation spectra")
parser.add_argument("-c", "--channel", type=str, default="b",
                    dest="channel", help="annihilation channel to optimize energy range for")
parser.add_argument("-m", "--mass", type=int, default=10,
                    dest="mass", help="WIMP mass to optimize energy range for")
parser.add_argument("-p", "--part", type=int, default=0,
                    dest="part", help="subset of MC files to use")
parser.add_argument("-o", "--outfolder", type=str, default="/mnt/home/priesbr1/DM_Search/data/energy_optimization/",
                    dest="outfolder", help="output directory for energy range optimization results")
args = parser.parse_args()

datafolder = args.datafolder
dataset = "674"
mass = args.mass
channel = args.channel
part = str(args.part)
outfolder = args.outfolder

nue_list = glob.glob(datafolder + "nue/12" + dataset + "/*" + part + ".i3.bz2")
numu_list = glob.glob(datafolder + "numu/14" + dataset + "/*" + part + ".i3.bz2")
nutau_list = glob.glob(datafolder + "nutau/16" + dataset + "/*" + part + ".i3.bz2")
nue_files = len(nue_list)
numu_files = len(numu_list)
nutau_files = len(nutau_list)
Infile_List = nue_list + numu_list + nutau_list

tray = I3Tray()

signal_weight_array = [0.0 for x in range(0,500)]
background_weight_array = [0.0 for x in range(0,500)]

signal_weight_range_array = [0.0 for x in range(0,int(500*(500-1)/2))]
background_weight_range_array = [0.0 for x in range(0,int(500*(500-1)/2))]

spectra = np.load(args.spectra, allow_pickle=True)
spectra = spectra.item()
fluxes = spectra[channel][mass]["EdNdE"]/spectra[channel][mass]["Energy"]
energies = spectra[channel][mass]["Energy"]
ntype_map = {0:"numu", 1:"numubar", 2:"nue", 3:"nuebar", 4:"nutau", 5:"nutabar", 6:"numu", 7:"numubar", 8:"nue", 9:"nuebar", 10:"nutau", 11:"nutaubar"}

def find_nearest(array,value):
    return (np.abs(array-value)).argmin()

def calculate_angle(recozen,recoazi,mczen,mcazi):
    return np.arccos(np.sin(recozen)*np.sin(mczen)*np.cos(recoazi-mcazi) + np.cos(recozen)*np.cos(mczen))

def calculate_signal_pdfs(frame):
    mctree = frame["I3MCTree"]
    mcparticle = get_most_energetic_primary(mctree)
    mcweightdict = frame["I3MCWeightDict"]
    if abs(mcparticle.type) == 14:
        infile_list_length = numu_files
        ftype = 0
    if abs(mcparticle.type) == 12:
        infile_list_length = nue_files
        ftype = 1
    if abs(mcparticle.type) == 16:
        infile_list_length = nutau_files
        ftype = 2
    if mcparticle.type > 0:
        typeweight = 0.7
        qtype = 0
    if mcparticle.type < 0:
        typeweight = 0.3
        qtype = 1
    if mcweightdict["InteractionType"] == 1.0:
        itype = 0
    if mcweightdict["InteractionType"] == 2.0:
        itype = 1
    if mcweightdict["InteractionType"] == 0.0:
        return False
    ntype = 6*itype + 2*ftype + qtype
    flv = ntype_map[ntype]
    
    pegleg_neutrino = frame["IC86_Dunkman_L6_PegLeg_MultiNest8D_NumuCC"]
    reco_zenith = pegleg_neutrino.dir.zenith
    reco_azimuth = pegleg_neutrino.dir.azimuth
    pid = frame["IC86_Dunkman_L6"]["delta_LLH"]
    mc_zenith = mcparticle.dir.zenith
    mc_azimuth = mcparticle.dir.azimuth
    total_space_angle = calculate_angle(reco_zenith,reco_azimuth,mc_zenith,mc_azimuth)
    mc_energy = mcparticle.energy
    z = mc_energy/mass
    
    for e in range(0,500):
        if pegleg_neutrino.energy <= e+1:
            background_weight_array[e] += frame["NeutrinoWeights"]["OscillatedRate"]/(infile_list_length*typeweight/0.5)
    counter = 0
    for lower_limit in range(0,500):
        for upper_limit in range(lower_limit+1,500):
            if pegleg_neutrino.energy >= lower_limit and pegleg_neutrino.energy <= upper_limit:
                background_weight_range_array[counter] += frame["NeutrinoWeights"]["OscillatedRate"]/(infile_list_length*typeweight/0.5)
            counter += 1
    
    if z <= 1:
        ebin = energies[find_nearest(energies,mc_energy)]
        flux = fluxes[energies == ebin]
        
        weight = flux*mcweightdict["OneWeight"]/(mass*mcweightdict["NEvents"]*infile_list_length*typeweight)
        for e in range(0,500):
            if pegleg_neutrino.energy <= e+1:
                signal_weight_array[e] += weight
        
        counter = 0
        for lower_limit in range(0,500):
            for upper_limit in range(lower_limit+1,500):
                if pegleg_neutrino.energy >= lower_limit and pegleg_neutrino.energy <= upper_limit:
                    signal_weight_range_array[counter] += weight
                counter += 1

tray.AddModule("I3Reader", "reader", filenamelist=Infile_List)
tray.AddModule(event_cuts, "event_cuts")
tray.AddModule(calculate_signal_pdfs, "calculate_signal_pdfs")
tray.AddModule("TrashCan","Done")

t1 = time.time()
tray.Execute()
t2 = time.time()
tray.Finish()

filename = outfolder + "weight_vs_cutoff_" + channel + "_" + str(mass) + "_part" + part + ".txt"
f = open(filename,"wb")
pickle.dump(signal_weight_array,f)
pickle.dump(background_weight_array,f)
pickle.dump(signal_weight_range_array,f)
pickle.dump(background_weight_range_array,f)
f.close()

print("Finished processing in %s seconds"%(t2-t1))
print("Saving output to %s"%filename)
