import numpy as np
import argparse
import time
import math

from skylab.ps_llh import PointSourceLLH, MultiPointSourceLLH
from skylab.ps_injector import PointSourceInjector
from skylab.llh_models import ClassicLLH, EnergyLLH
from skylab.spectral_models import SplinedSpectrum, PowerLaw
from skylab.datasets import Datasets
from skylab.sensitivity_utils import DeltaChiSquare, estimate_sensitivity
from astropy import units as u

GeV = 1
TeV = 1000 * GeV

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument("-p", "--spectra", type=str, default="/mnt/home/priesbr1/DM_Search/data/annihilation_spectra/DM_CirelliSpectrum_dict_neutrinos_oscillated.npy",
                    dest="spectra", help="file containing dark matter spectra for each channel and mass")
parser.add_argument("-m", "--mass", type=int, default=10,
                    dest="mass", help="WIMP mass to use")
parser.add_argument("-c", "--channel", type=str, default="b",
                    dest="channel", help="annihilation channel to use")
parser.add_argument("-u", "--sources", type=str, default="/mnt/home/priesbr1/DM_Search/data/analysis_sources_ra_dec_jfactors.txt",
                    dest="sources", help="file containing RA, dec, and J-factor information for each source")
parser.add_argument("-j", "--J_type", type=str, choices=["0.1","0.2","0.5","10","min","max"], default="max",
                    dest="J_type", help="J-factor values to use based on opening half-angle")
parser.add_argument("-b", "--background_trials", type=str, default="/mnt/home/priesbr1/DM_Search/data/trials_results/trials_bkg_360deg_Jmax/",
                    dest="background_trials", help="folder containing results from background trials")
parser.add_argument("-s", "--signal_trials", type=str, default="/mnt/home/priesbr1/DM_Search/data/trials_results/trials_sig_360deg_Jmax_distributed/",
                    dest="signal_trials", help="folder containing results from signal trials")
parser.add_argument("-r", "--repo_path", type=str, default="/mnt/research/IceCube/datasets/",
                    dest="repo_path", help="repository path for SkyLab dataset")
parser.add_argument("-o", "--outfile", type=str, default="/mnt/home/priesbr1/DM_Search/data/cross_section_results/cross_section_results_360deg_Jmax_distributed_erange.npy",
                    dest="outfile", help="outfile to save results to")
args = parser.parse_args()

# Set up data
sample = 'PointSourceDRAGON_v001p00'  # PointSourceTracks_v00?p0?, GFU, PointSourceDRAGON_v001p00, etc.
seasons = ['IC86, 2011', 'IC86, 2012', 'IC86, 2013', 'IC86, 2014', 'IC86, 2015', 'IC86, 2016', 'IC86, 2017']  # years

def extract_source_info(sources, J_type, J_indices_map):
    names = []
    ra = np.array([])
    dec = np.array([])
    half_angles = []
    J_factors = np.array([])

    J_indices = J_indices_map[J_type]

    for i in range(len(sources)):
        for idx in J_indices:
            if not np.isnan(sources[i][idx]):
                names.append(sources[i][0])
                ra = np.append(ra, sources[i][1]*np.pi/180.)
                dec = np.append(dec, sources[i][2]*np.pi/180.)
                half_angles.append(float(list(J_indices_map.keys())[list(J_indices_map.values()).index([idx])]))
                J_factors = np.append(J_factors, sources[i][idx])
                break

    J_factors = 10**J_factors

    return names, ra, dec, half_angles, J_factors

J_indices_map = {"0.1": [6], "0.2": [9], "0.5": [12], "10": [15], "min": [6,9,12,15], "max": [15,12,9,6]}

# Source location and information
# Source encoding - 0:name; 1:RA; 2:Dec; 3,4,5:r_h; 6,7,8:J0.1; 9,10,11:J0.2; 12,13,14:J0.5; 15,16,17:J10
sources = np.genfromtxt(args.sources, delimiter=",", dtype=("U20",float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float))
names, ra, dec, half_angles, J_factors = extract_source_info(sources, args.J_type, J_indices_map)
J_factor_sum = np.sum(J_factors)

max_name_length = np.max([len(name) for name in names])

print("Sources used in analysis (%i): "%len(names))
for source in list(zip(names, ra, dec, half_angles, J_factors)):
    print("  - {0:<{1}} (RA = {2:>6.2f} deg, dec = {3:>6.2f} deg), J_factor({4:.1f} deg) = {5:.4f}".format(
          tuple(source)[0], max_name_length, tuple(source)[1]*180/np.pi, tuple(source)[2]*180/np.pi, tuple(source)[3], tuple(source)[4]))

# Create temporary arrays of energy/flux for channels/masses
DM_spectra = np.load(args.spectra, allow_pickle=True)
DM_spectra = DM_spectra.item()

data = [None, None]
data[0] = np.array(DM_spectra[args.channel][args.mass]["Energy"])
data[1] = np.array(DM_spectra[args.channel][args.mass]["EdNdE"] / DM_spectra[args.channel][args.mass]["Energy"])

# Create splined spectrum/spectra
scale_factor = 1/(8*np.pi) / (args.mass**2)
mask = np.where(data[1] > 0)
#spectra = SplinedSpectrum(energies=data[0][mask], fluxes=data[1][mask], E0=data[0][mask][0], order=1)
spectra = SplinedSpectrum(energies=data[0][mask], fluxes=data[1][mask], E0=0.7*args.mass, order=1)

# Load data
llh = []
multillh = MultiPointSourceLLH()
baseline = 0.

dragon = Datasets[sample]
dragon.set_repository_path(args.repo_path)
for season in seasons:
    exp, mc, livetime = dragon.season(season)
    sinDec_bins = dragon.sinDec_bins(season)
    energy_bins = dragon.energy_bins(season)
    
    msg = "   - % 15s (" % season
    msg += "livetime %7.2f days, %6d events" % (livetime, exp.size)
    msg += ", mjd0 %.2f" % min(exp['time'])
    msg += ", mjd1 %.2f)" % max(exp['time'])
    print(msg)
    
    # Set up likelihood
    llh_model = EnergyLLH(twodim_bins=[energy_bins, sinDec_bins], allow_empty=True, spectrum=spectra)
    llh.append(PointSourceLLH(exp, mc, livetime, llh_model, scramble=True, mode="all", nsource_bounds=(0,5000), nsource=len(names))) # multiple sources
    multillh.add_sample(season, llh[-1])
    baseline += llh_model.spectrum(spectra.E0)

baseline /= len(seasons)
print("Average baseline:", baseline)

erange = tuple([min(data[0][mask]), args.mass])

# Injector
inj = PointSourceInjector(spectrum=spectra, E0=spectra.E0, e_range=erange)
inj.fill(dec, multillh.exp, multillh.mc, multillh.livetime, src_w=J_factors)

bkg_file = args.background_trials + "bkg_trials_" + args.channel + "_" + str(args.mass) + ".txt"
bkg_trials = np.genfromtxt(bkg_file, delimiter="\t")
bkg_ns = []
bkg_TS = []
for i in range(len(bkg_trials)):
    bkg_ns.append(bkg_trials[i][-2])
    bkg_TS.append(bkg_trials[i][-1])
bkg_ns = np.array(bkg_ns, dtype=int)
bkg_TS = np.array(bkg_TS)
med_bkg_TS = np.median(bkg_TS)
print("Median background TS:", med_bkg_TS)
#bkg_TS_90 = np.percentile(bkg_TS, 90)
#print("90th percentile background TS:", bkg_TS_90)

sig_file = args.signal_trials + "sig_trials_" + args.channel + "_" + str(args.mass) + ".txt"
sig_trials = np.genfromtxt(sig_file, delimiter="\t")
sig_ni = []
sig_ns = []
sig_TS = []
for i in range(len(sig_trials)):
    sig_ni.append(sig_trials[i][0])
    sig_ns.append(sig_trials[i][-2])
    sig_TS.append(sig_trials[i][-1])
sig_ni = np.array(sig_ni, dtype=int)
sig_ns = np.array(sig_ns, dtype=int)
sig_TS = np.array(sig_TS)

for ni in sorted(np.unique(sig_ni)):
    ni_TS = np.array(sig_TS[sig_ni == ni])
    TS_frac = len(np.where(ni_TS > med_bkg_TS)[0])/len(ni_TS)
    #TS_frac = len(np.where(ni_TS > bkg_TS_90)[0])/len(ni_TS)
    print("Fraction of trials above threshold at %i injected events: %s"%(ni, TS_frac))
    if TS_frac > 0.9:
    #if TS_frac > 0.5:
        flux = inj.mu2flux(ni)
        sigma_v = flux/(baseline * J_factor_sum * scale_factor)
        print("Annihilation cross section for (%s,%i): %s cm^3 s^-1"%(args.channel, args.mass, str(sigma_v)))
        break

try:
    cross_sections = np.load(args.outfile, allow_pickle=True)
    cross_sections = cross_sections.item()
except:
    cross_sections = dict()

if args.channel not in cross_sections.keys():
    cross_sections[args.channel] = dict()
cross_sections[args.channel][args.mass] = sigma_v

np.save(args.outfile, cross_sections)
