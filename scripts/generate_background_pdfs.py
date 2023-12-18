# Adapted from IceCube Low-Energy Solar WIMP analysis
# https://github.com/IceCubeOpenSource/IC86_LE_solarDM/blob/main/scripts/generate_bkg_pdf.py

import numpy as np
import glob
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-d", "--data", type=str, default=None,
                    dest="data", help="file(s) containing data to generate background pdf from")
parser.add_argument("-s", "--sources", type=str, default="/mnt/home/priesbr1/DM_Search/data/analysis_sources_ra_dec_jfactors.txt",
                    dest="sources", help="file containing RA, dec, and J-factor information for each source")
parser.add_argument("-c", "--channel", type=str, default="b",
                    dest="channel", help="annihilation channel to generate background PDFs for")
parser.add_argument("-m", "--mass", type=int, default=10,
                    dest="mass", help="WIMP mass to generate background PDFs for")
parser.add_argument("-b", "--band_size", type=int, default=5,
                    dest="band_size", help="size of declination band in degrees")
parser.add_argument("-e", "--seed", type=int, default=None,
                    dest="seed", help="RNG seed for reproducibility")
parser.add_argument("-o", "--output", type=str, default="/mnt/home/priesbr1/DM_Search/data/pdfs/DRAGON_background_pdfs_5deg.npy",
                    dest="outfile", help="file to save background PDFs to (will also try to load previously-calculated PDFs)")
args = parser.parse_args()

if (args.seed != None):
    np.random.seed(args.seed)

def calculate_angle(dec_source, ra_source, dec_nu, ra_nu):
    return np.arccos(np.cos(dec_source)*np.cos(dec_nu)*np.cos(ra_source - ra_nu) + np.sin(dec_source)*np.sin(dec_nu))

def calculate_background_pdf(datafiles, source_ra, source_dec, source_name, energy_range, nbins):
    bkg_angles = []
    
    for filename in datafiles:
        data = np.load(filename)
        
        for i in range(len(data["run"])):
            nu_time = data["time"][i]
            year = int(filename[filename.find("20"):filename.find("20")+4])
            nu_dec = data["dec"][i]
            nu_energy = 10**data["logE"][i]
            if year > 2015:
                nu_energy *= 1.04
            
            if (nu_dec >= source_dec - args.band_size/2*np.pi/180.) and (nu_dec <= source_dec + args.band_size/2*np.pi/180.):
                if (energy_range[0] <= nu_energy) and (nu_energy <= energy_range[1]):
                    for k in range(30):
                        scrambled_ra = np.random.uniform(0, 2*np.pi)
                        source_angle = calculate_angle(source_dec, source_ra, nu_dec, scrambled_ra)
                        bkg_angles.append(source_angle)
    
        print("Finished data from %i"%year)
    
    bkg_hist = []
    hist, bins = np.histogram(bkg_angles, bins=nbins, range=(0, np.pi))
    bkg_hist.append((hist, bins))
    
    print("Finished creating background PDFs for %s"%source_name)
    return bkg_hist

def get_energy_ranges(channel, mass):
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
    
    if mass not in energy_list:
        raise TypeError("Energy range not available for mass %i -- mass must be one of %s"%(mass, energy_list))
    return [energy_list, energy_ranges, nbins]

if "*" in args.data or "?" in args.data:
    filenames = sorted(glob.glob(args.data))
else:
    filenames = [args.data]
print("Using data from:")
for filename in filenames:
    print(" - "+str(filename))

# Source encoding - 0:name; 1:RA; 2:Dec; 3,4,5:r_h; 6,7,8:J0.1; 9,10,11:J0.2; 12,13,14:J0.5; 15,16,17:J10
sources = np.genfromtxt(args.sources, delimiter=",", dtype=("U20",float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float))
names = []
ra  = np.array([])
dec = np.array([])
for i in range(len(sources)):
    names.append(sources[i][0])
    ra = np.append(ra, sources[i][1]*np.pi/180.)
    dec = np.append(dec, sources[i][2]*np.pi/180.)

try:
    pdfs = np.load(args.outfile, allow_pickle=True)
    pdfs = pdfs.item()
    print("Found existing background PDFs at %s"%args.outfile)
except:
    print("No existing background PDFs found")
    pdfs = dict()

for i, source in enumerate(names):
    if names[i] not in pdfs.keys():
        pdfs[source] = dict()

    energy_list, energy_ranges, nbins = get_energy_ranges(args.channel, args.mass)
    if args.channel not in pdfs[source].keys():
        pdfs[source][args.channel] = dict()

    bkg_hist = calculate_background_pdf(filenames, ra[i], dec[i], source, energy_ranges[energy_list.index(args.mass)], nbins)
    if args.mass not in pdfs[source][args.channel].keys():
        pdfs[source][args.channel][args.mass] = dict()
    pdfs[source][args.channel][args.mass]["Counts"] = bkg_hist[0][0]
    pdfs[source][args.channel][args.mass]["Bins"] = bkg_hist[0][1]

np.save(args.outfile, pdfs)
