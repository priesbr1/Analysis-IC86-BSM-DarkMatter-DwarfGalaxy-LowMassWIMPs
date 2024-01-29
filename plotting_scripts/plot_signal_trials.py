import numpy as np
import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
import argparse
import os, sys
import subprocess
import glob

parser = argparse.ArgumentParser()
parser.add_argument("-t", "--trials", type=str, default="/mnt/home/priesbr1/DM_Search/data/trials_results/trials_sig_360deg_Jmax_distributed/sig_trials_*_*.txt",
                    dest="trials", help="file(s) containing signal trials to plot")
parser.add_argument("-s", "--sources", type=str, default="/mnt/home/priesbr1/DM_Search/data/analysis_sources_ra_dec_jfactors.txt",
                    dest="sources", help="file containing RA, dec, and J-factor information for each source")
parser.add_argument("-j", "--J_type", type=str, choices=["0.1","0.2","0.5","10","min","max"], default="max",
                    dest="J_type", help="J-factor values to use based on opening half-angle")
parser.add_argument("-o", "--outfolder", type=str, default="/mnt/scratch/priesbr1/DM_Search/Plots/signal_trials/TS_360deg_Jmax_distributed/",
                    dest="outfolder", help="folder to save output plots to")
args = parser.parse_args()

sig_trials = sorted(glob.glob(args.trials))
if isinstance(sig_trials,str):
    sig_trials = [sig_trials]

if not os.path.isdir(args.outfolder):
    os.mkdir(args.outfolder)

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
J_total = np.sum(J_factors)

for trial_set in sig_trials:
    assert "J"+args.J_type in trial_set
    
    get_header = subprocess.run(["head", "-1", trial_set], check=True, capture_output=True)
    get_N = subprocess.run(["wc", "-w"], input=get_header.stdout, check=True, capture_output=True)
    N = int(get_N.stdout.decode("utf-8").strip())
    assert int((N-4)/2) == len(names)
    
    print("Plotting TS and ns from %s"%trial_set)
    channel_mass = trial_set[trial_set.find("sig_trials")+11:-4]
    channel, mass = channel_mass.split("_")
    mass = int(mass)
    
    trials = np.genfromtxt(trial_set, delimiter="\t", dtype=None)
    ni = []
    ns_sources = []
    TS_sources = []
    ns_total = []
    TS_total = []
    for i in range(len(trials)):
        ni.append(trials[i][0])
        ns_sources.append([])
        TS_sources.append([])
        for j in range(len(names)):
            ns_sources[i].append(trials[i][2*j - 1])
            TS_sources[i].append(trials[i][2*j])
        ns_sources[i] = np.array(ns_sources[i])
        TS_sources[i] = np.array(TS_sources[i])
        ns_total.append(trials[i][-2])
        TS_total.append(trials[i][-1])
    
    ni = np.array(ni)
    ns_sources = np.array(ns_sources)
    TS_sources = np.array(TS_sources)
    ns_total = np.array(ns_total)
    TS_total = np.array(TS_total)
    
    ni_vals = sorted(np.unique(ni))
    
    for k in range(len(names)):
        ni_source = [int(np.round(n_inj*J_factors[k]/J_total)) for n_inj in ni_vals]
        ns_med, ns_low, ns_high = [], [], []
        TS_med, TS_low, TS_high = [], [], []
        
        for n_inj in ni_vals:
            ns_med.append(np.percentile(ns_sources[:,k][ni == n_inj], 50))
            ns_low.append(np.percentile(ns_sources[:,k][ni == n_inj], 84))
            ns_high.append(np.percentile(ns_sources[:,k][ni == n_inj], 16))
            TS_med.append(np.percentile(TS_sources[:,k][ni == n_inj], 50))
            TS_low.append(np.percentile(TS_sources[:,k][ni == n_inj], 84))
            TS_high.append(np.percentile(TS_sources[:,k][ni == n_inj], 16))
    
        plt.figure()
        plt.plot(ni_source, ns_med, label="Data")
        plt.fill_between(ni_source, ns_low, ns_high, color="blue", alpha=0.2)
        plt.plot(ni_source, ni_source, color="black", linestyle="dashed", label=r"$n_{\mathrm{signal}} = n_{\mathrm{injected}}$")
        if (np.max(ni_vals)/np.min(ni_vals) > 1000):
            plt.xscale("log")
            plt.yscale("log")
        plt.xlabel(r"$n_{\mathrm{injected}}$")
        plt.ylabel(r"$n_{\mathrm{signal}}$")
        plt.title("%s Median Signal Recovery for (%s,%i) WIMPs"%(names[k],channel,mass))
        plt.legend(loc="best")
        plt.tight_layout()
        plt.savefig(args.outfolder + "ns_" + names[k] + "_" + channel + "_" + str(mass) + ".png")
        plt.close()
    
        plt.figure()
        plt.plot(ni_source, TS_med)
        plt.fill_between(ni_source, TS_low, TS_high, color="blue", alpha=0.2)
        if (np.max(ni_vals)/np.min(ni_vals) > 100):
            plt.xscale("log")
        plt.xlabel(r"$n_{\mathrm{injected}}$")
        plt.ylabel("Test Statistic")
        plt.title("%s Median Test Statistics for (%s,%i) WIMPs"%(names[k],channel,mass))
        plt.tight_layout()
        plt.savefig(args.outfolder + "TS_" + names[k] + "_" + channel + "_" + str(mass) + ".png")
        plt.close()
    
    ns_med, ns_low, ns_high = [], [], []
    TS_med, TS_low, TS_high = [], [], []
    
    for n_inj in ni_vals:
        ns_med.append(np.percentile(ns_total[ni == n_inj], 50))
        ns_low.append(np.percentile(ns_total[ni == n_inj], 84))
        ns_high.append(np.percentile(ns_total[ni == n_inj], 16))
        TS_med.append(np.percentile(TS_total[ni == n_inj], 50))
        TS_low.append(np.percentile(TS_total[ni == n_inj], 84))
        TS_high.append(np.percentile(TS_total[ni == n_inj], 16))
        
    plt.figure()
    plt.plot(ni_vals, ns_med, label="Data")
    plt.fill_between(ni_vals, ns_low, ns_high, color="blue", alpha=0.2)
    plt.plot(ni_vals, ni_vals, color="black", linestyle="dashed", label=r"$n_{\mathrm{signal}} = n_{\mathrm{injected}}$")
    if (np.max(ni_vals)/np.min(ni_vals) > 100):
        plt.xscale("log")
        plt.yscale("log")
    plt.xlabel(r"$n_{\mathrm{injected}}$")
    plt.ylabel(r"$n_{\mathrm{signal}}$")
    plt.title("Median Signal Recovery for (%s,%i) WIMPs"%(channel,mass))
    plt.legend(loc="best")
    plt.tight_layout()
    plt.savefig(args.outfolder + "ns_" + channel + "_" + str(mass) + ".png")
    plt.close()
    
    plt.figure()
    plt.plot(ni_vals, TS_med)
    plt.fill_between(ni_vals, TS_low, TS_high, color="blue", alpha=0.2)
    if (np.max(ni_vals)/np.min(ni_vals) > 100):
        plt.xscale("log")
    plt.xlabel(r"$n_{\mathrm{injected}}$")
    plt.ylabel("Test Statistic")
    plt.title("Median Test Statistics for (%s,%i) WIMPs"%(channel,mass))
    plt.tight_layout()
    plt.savefig(args.outfolder + "TS_" + channel + "_" + str(mass) + ".png")
    plt.close()
    
print("Finished plotting")
