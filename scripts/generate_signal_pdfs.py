# Adapted from IceCube Low-Energy Solar WIMP analysis
# https://github.com/IceCubeOpenSource/IC86_LE_solarDM/blob/main/scripts/solarWIMP_pdfs_ch5.py

import numpy as np
import glob
import argparse
from scipy.interpolate import interp1d

parser = argparse.ArgumentParser()
parser.add_argument("-d", "--data", type=str, default=None,
                    dest="data", help="MC file(s) containing data to generate background pdf from")
parser.add_argument("-p", "--spectra", type=str, default="/mnt/home/priesbr1/DM_Search/data/annihilation_spectra/DM_CirelliSpectrum_dict_neutrinos_oscillated.npy",
                    dest="spectra", help="file containing WIMP anniliation spectra")
parser.add_argument("-c", "--channel", type=str, default="b",
                    dest="channel", help="annihilation channel to generate background PDF for")
parser.add_argument("-o", "--output", type=str, default="/mnt/home/priesbr1/DM_Search/data/pdfs/DRAGON_signal_pdfs.npy",
                    dest="outfile", help="file to save signal PDFs to (will also try to load previously-calculated PDFs)")
args = parser.parse_args()

def calculate_angle(dec_source, ra_source, dec_nu, ra_nu):
    return np.arccos(np.cos(dec_source)*np.cos(dec_nu)*np.cos(ra_source - ra_nu) + np.sin(dec_source)*np.sin(dec_nu))

def calculate_signal_pdf(datafiles, spectra, channel, mass, energy_ranges, nbins):
    signal_angles = [[] for i in range(len(energy_ranges))]
    signal_weights = [[] for i in range(len(energy_ranges))]
    spline = interp1d(spectra[channel][mass]["Energy"], spectra[channel][mass]["EdNdE"]/spectra[channel][mass]["Energy"])
    
    for filename in datafiles:
        data = np.load(filename)
        
        for i in range(len(data["run"])):
            nu_time = data["time"][i]
            year = int(filename[filename.find("20"):filename.find("20")+4])
            nu_dec = data["dec"][i]
            nu_ra = data["ra"][i]
            nu_energy = 10**data["logE"][i]
            mc_dec = data["trueDec"][i]
            mc_ra = data["trueRa"][i]
            mc_energy = data["trueE"][i]
            if year > 2015:
                nu_energy *= 1.04
            
            for j, pair in enumerate(energy_ranges):
                if (pair[0] <= nu_energy) and (nu_energy <= pair[1]):
                    angular_separation = calculate_angle(nu_dec, nu_ra, mc_dec, mc_ra)
                    signal_angles[j].append(angular_separation)
                    if mc_energy <= mass:
                        signal_weights[j].append(spline(mc_energy))
                    else:
                        signal_weights[j].append(0.0)
        
        print("Finished data from %i"%year)
    
    signal_hists = []
    for i in range(len(signal_angles)):
        hist, bins = np.histogram(signal_angles[i], bins=nbins, weights=signal_weights[i], range=(0, np.pi))
        signal_hists.append((hist, bins))
    
    print("Finished creating signal PDFs for (%s,%i) WIMPs"%(channel,mass))
    return signal_hists

def get_energy_ranges(channel):
    if channel not in ["b", "Tau", "Mu", "W", "Nu"]:
        raise TypeError("Channel (%s) not recognized -- must be one of b/Mu/Tau/W/Nu"%channel)

    if channel == "b":
        energy_list = [10, 20, 30, 50, 80, 100, 150, 200, 300]
        energy_ranges = [(0,12), (0,18), (0,25), (0,38), (0,57), (0,73), (0,108), (3,152), (6,256)]
        nbins = 36
    if channel == "Mu":
        energy_list = [5, 10, 20, 30, 50, 80, 100, 150, 200, 300]
        energy_ranges = [(0,9), (0,16), (3,30), (6,43), (16,49), (24,107), (31,129), (48,199), (66,266), (100,389)]
        nbins = 36
    if channel == "Tau":
        energy_list = [5, 10, 20, 30, 50, 80, 100, 150, 200, 300]
        energy_ranges = [(0,9), (0,16), (0,30), (4,44), (12,71), (20,108), (27,138), (40,205), (57,272), (92,389)]
        nbins = 36
    if channel == "W":
        energy_list = [90, 100, 150, 200, 300]
        energy_ranges = [(27,108), (31,128), (54,207), (82,298), (120,448)]
        nbins = 36
    if channel == "Nu":
        energy_list = [5, 10, 20, 30, 50, 80, 100, 150, 200, 300]
        energy_ranges = [(1,12), (0,21), (12,39), (19,53), (34,83), (56,128), (70,159), (113,234), (148,298), (234,441)]
        nbins = 36

    return [energy_list, energy_ranges, nbins]

if "*" in args.data or "?" in args.data:
    filenames = sorted(glob.glob(args.data))
else:
    filenames = [args.data]
print("Using data from:")
for filename in filenames:
    print(" - "+str(filename))

spectra = np.load(args.spectra, allow_pickle=True)
spectra = spectra.item()

try:
    pdfs = np.load(args.outfile, allow_pickle=True)
    pdfs = pdfs.item()
    print("Found existing signal PDFs at %s"%args.outfile)
except:
    print("No existing signal PDFs found")
    pdfs = dict()

energy_list, energy_ranges, nbins = get_energy_ranges(args.channel)
if args.channel not in pdfs.keys():
    pdfs[args.channel] = dict()
    
for j, mass in enumerate(energy_list):
    if mass not in pdfs[args.channel].keys():
        pdfs[args.channel][mass] = dict()
    
    signal_hists = calculate_signal_pdf(filenames, spectra, args.channel, mass, energy_ranges, nbins)
    pdfs[args.channel][mass]["Counts"] = signal_hists[j][0]
    pdfs[args.channel][mass]["Bins"] = signal_hists[j][1]

np.save(args.outfile, pdfs)
