import numpy as np
import glob
import argparse
import scipy
from scipy import special, stats
import time
import os,sys

def generate_rand_from_pdf(pdf, x_grid, n=1000):
    cdf = np.cumsum(pdf)
    cdf = cdf / cdf[-1]
    values = np.random.rand(n)
    value_bins = np.searchsorted(cdf, values)
    random_from_cdf = x_grid[value_bins]
    return random_from_cdf

def calculate_angle(dec_source, ra_source, dec_nu, ra_nu):
    return np.arccos(np.cos(dec_source)*np.cos(dec_nu)*np.cos(ra_source - ra_nu) + np.sin(dec_source)*np.sin(dec_nu))

def p2sigma(p):
  """Convert 1-sided p-value (probability for observing >= x) to
  significance defined by equivalent integral on a normal distribution

  Parameters
  ----------
  p : float, np.ndarray
    1-sided p-value for observing >= a given measurement

  Returns
  -------
  1-sided significance from a Gaussian curve
  """

  if isinstance(p, np.ndarray):
    s = np.empty(p.size, dtype=float)

    mask = (p >= 1)
    s[mask] = -1000
    s[~mask] = scipy.special.erfcinv(2 * p[~mask]) * np.sqrt(2)
    return s
  # END if (array)

  if p >= 1:
    return -1000
  return scipy.special.erfcinv(2 * p) * np.sqrt(2)

def ts2sigma(ts, ndf=1.0, eta=0.5, scale=1.0):
  """Convert TS to 1-sided significance from a Gaussian curve.

  Parameters
  ----------
  ts : float
    Test statistic (TS) value
  ndf : float
    Effective degrees of freedom, should match number of free parameters
    in the likelihood in the limit of large statistics
  eta : float
    Fraction of TS > 0
  scale : float
    Scaling applied to TS values

  Returns
  -------
  1-sided significance from a Gaussian curve
  """
  p = ts2p(ts, ndf, eta, scale)
  return p2sigma(p)

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

parser = argparse.ArgumentParser()
parser.add_argument("-d", "--data", type=str, default="/mnt/research/IceCube/datasets/ps_DRAGON/version-001-p00/IC86_201?_exp.npy",
                    dest="data", help="files containing neutrino data")
parser.add_argument("-b", "--bkg_pdfs", type=str, default="/mnt/home/priesbr1/DM_Search/data/pdfs/DRAGON_background_pdfs_360deg.npy",
                    dest="bkg_pdfs", help="file containing background PDFs for each mass, channel, and source")
parser.add_argument("-p", "--sig_pdfs", type=str, default="/mnt/home/priesbr1/DM_Search/data/pdfs/DRAGON_signal_pdfs.npy",
                    dest="sig_pdfs", help="file containing signal PDFs for each mass and channel")
parser.add_argument("-s", "--sources", type=str, default="/mnt/home/priesbr1/DM_Search/data/analysis_sources_ra_dec_jfactors.txt",
                    dest="sources", help="file containing RA, dec, and J-factor information for each source")
parser.add_argument("-j", "--J_type", type=str, choices=["0.1","0.2","0.5","10","min","max"], default="max",
                    dest="J_type", help="J-factor values to use based on opening half-angle")
parser.add_argument("-c", "--channel", type=str, default="b",
                    dest="channel", help="annhiliation channel to run background trials for")
parser.add_argument("-m", "--mass", type=int, default=10,
                    dest="mass", help="WIMP mass to rub background trials for")
parser.add_argument("-t", "--background_trials", type=str, default="/mnt/home/priesbr1/DM_Search/data/trials_results/trials_bkg_360deg_Jmax_100kbe/",
                    dest="background_trials", help="folder containing background trials results")
parser.add_argument("-o", "--outfolder", type=str, default="/mnt/home/priesbr1/DM_Search/data/trials_results/trials_unblinded_360deg_Jmax/",
                    dest="outfolder", help="output folder for trials results")
args = parser.parse_args()

energy_list, energy_ranges, nbins = get_energy_ranges(args.channel, args.mass)
energy_cut = energy_ranges[energy_list.index(args.mass)]

datafiles = sorted(glob.glob(args.data))

nu_decs = []
nu_ras = []

for filename in datafiles:
    season = np.load(filename)
    nu_decs.extend(season["dec"][(season["logE"] >= np.log10(energy_cut[0])) & (season["logE"] <= np.log10(energy_cut[1]))])
    nu_ras.extend(season["ra"][(season["logE"] >= np.log10(energy_cut[0])) & (season["logE"] <= np.log10(energy_cut[1]))])

bkg_pdfs = np.load(args.bkg_pdfs, allow_pickle=True)
bkg_pdfs = bkg_pdfs.item()
sig_pdfs = np.load(args.sig_pdfs, allow_pickle=True)
sig_pdfs = sig_pdfs.item()

def extract_source_info(sources, J_type, J_indices_map):
    names = []
    ra = np.array([])
    dec = np.array([])
    J_factors = np.array([])

    J_indices = J_indices_map[J_type]

    for i in range(len(sources)):
        for idx in J_indices:
            if not np.isnan(sources[i][idx]):
                names.append(sources[i][0])
                ra = np.append(ra, sources[i][1]*np.pi/180.)
                dec = np.append(dec, sources[i][2]*np.pi/180.)
                J_factors = np.append(J_factors, sources[i][idx])
                break

    J_factors = 10**J_factors

    return names, ra, dec, J_factors

J_indices_map = {"0.1": [6], "0.2": [9], "0.5": [12], "10": [15], "min": [6,9,12,15], "max": [15,12,9,6]}

# Source encoding - 0:name; 1:RA; 2:Dec; 3,4,5:r_h; 6,7,8:J0.1; 9,10,11:J0.2; 12,13,14:J0.5; 15,16,17:J10
sources = np.genfromtxt(args.sources, delimiter=",", dtype=("U20",float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float))
names, ra, dec, J_factors = extract_source_info(sources, args.J_type, J_indices_map)

background_file = args.background_trials + "bkg_trials_" + args.channel + "_" + str(args.mass) + ".txt"
bkg_trials = np.genfromtxt(background_file, delimiter="\t")

bkg_TS = np.zeros((len(bkg_trials), len(sources)+1), dtype=float)
for i in range(len(bkg_trials)):
    for j in range(len(sources)+1):
        bkg_TS[i,j] = bkg_trials[i][2*j+1]

outfile = "unblinded_trials_" + args.channel + "_" + str(args.mass) + ".txt"

t1 = time.time()

angles = np.zeros((len(nu_decs), len(sources)), dtype=float)
for j in range(len(sources)):
    angles[:,j] = calculate_angle(dec[j], ra[j], nu_decs, nu_ras)

N = len(nu_decs)
llhs = []

for j, source in enumerate(names):
    N = len(angles[:,j])

    source_llh = []
    
    bkg_bins = np.searchsorted(bkg_pdfs[source][args.channel][args.mass]["Bins"], angles[:,j], "right")-1 # Find left edges of bins
    bkg_probs = np.array(bkg_pdfs[source][args.channel][args.mass]["Counts"][bkg_bins], dtype=float) # Get probabilities for bins
    bkg_probs /= np.sum(bkg_pdfs[source][args.channel][args.mass]["Counts"]) # Normalize probabilities
    sig_bins = np.searchsorted(sig_pdfs[args.channel][args.mass]["Bins"], angles[:,j], "right")-1
    sig_probs = np.array(sig_pdfs[args.channel][args.mass]["Counts"][sig_bins], dtype=float)
    sig_probs /= np.sum(sig_pdfs[args.channel][args.mass]["Counts"])
    
    for ns in range(int(N/10)):
        ns = float(ns)
        combs = (ns/N)*sig_probs + (1-(ns/N))*bkg_probs
        source_llh.append(-1*np.sum(np.log(combs)))
    
    llhs.append(source_llh)

bfs = np.argmin(llhs, axis=1)
tss = []
for L,bf in zip(llhs,bfs):
    tss.append(2*(L[0] - L[bf]))
p_values = []
for j, ts in enumerate(tss):
    perc = scipy.stats.percentileofscore(bkg_TS[:,j], ts)
    p_values.append(1-(perc/100))
sigmas = [np.maximum(0.0,p2sigma(p)) for p in p_values]

bf = int(np.round(np.sum(bfs * J_factors/np.sum(J_factors))))
ts = np.sum(np.array(tss) * J_factors/np.sum(J_factors))
perc = scipy.stats.percentileofscore(bkg_TS[:,-1], ts)
p_value = 1-(perc/100)
sigma = np.maximum(0.0,p2sigma(p_value))

f = open(args.outfolder+outfile, "w")
f.write("# ")
for j, source in enumerate(names):
    f.write("%s_ns \t %s_TS \t %s_sigma \t "%(source, source, source))
f.write("ns \t TS \t sigma \n")

for j, source in enumerate(names):
    f.write("%s \t %s \t %s \t "%(bfs[j], tss[j], sigmas[j]))
f.write("%s \t %s \t %s \n"%(bf, ts, sigma))
f.close()

t2 = time.time()
print("Time: %.2f sec"%(t2-t1))
