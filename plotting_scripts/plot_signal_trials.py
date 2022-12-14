import numpy as np
import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
import argparse
import os, sys
import glob

parser = argparse.ArgumentParser()
parser.add_argument("-t", "--trials", type=str, default="/mnt/home/priesbr1/DM_Search/data/trials_results/trials_sig_360deg_Jmax_distributed/sig_trials_*_*.txt",
                    dest="trials", help="file(s) containing signal trials to plot")
parser.add_argument("-o", "--outfolder", type=str, default="/mnt/scratch/priesbr1/DM_Search/Plots/signal_trials/TS_360deg_Jmax_distributed/",
                    dest="outfolder", help="folder to save output plots to")
args = parser.parse_args()

sig_trials = sorted(glob.glob(args.trials))
if isinstance(sig_trials,str):
    sig_trials = [sig_trials]

if not os.path.isdir(args.outfolder):
    os.mkdir(args.outfolder)

for trial_set in sig_trials:
    print("Plotting TS and ns from %s"%trial_set)
    channel_mass = trial_set[trial_set.find("sig_trials")+11:-4]
    channel, mass = channel_mass.split("_")
    mass = int(mass)
    
    trials = np.genfromtxt(trial_set, delimiter="\t", dtype=None)
    ni = []
    ns = []
    TS = []
    for i in range(len(trials)):
        ni.append(trials[i][0])
        ns.append(trials[i][-2])
        TS.append(trials[i][-1])
    
    ni = np.array(ni)
    ns = np.array(ns)
    TS = np.array(TS)
    
    ni_vals = sorted(np.unique(ni))
    ns_med, ns_low, ns_high = [], [], []
    TS_med, TS_low, TS_high = [], [], []
    
    for n_inj in ni_vals:
        ns_med.append(np.percentile(ns[ni == n_inj], 50))
        ns_low.append(np.percentile(ns[ni == n_inj], 84))
        ns_high.append(np.percentile(ns[ni == n_inj], 16))
        TS_med.append(np.percentile(TS[ni == n_inj], 50))
        TS_low.append(np.percentile(TS[ni == n_inj], 84))
        TS_high.append(np.percentile(TS[ni == n_inj], 16))
        
    plt.figure()
    plt.plot(ni_vals, ns_med, label="Data")
    plt.fill_between(ni_vals, ns_low, ns_high, color="blue", alpha=0.2)
    plt.plot(ni_vals, ni_vals, color="black", linestyle="dashed", label=r"$n_{\mathrm{signal}} = n_{\mathrm{injected}}$")
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
    plt.xlabel(r"$n_{\mathrm{injected}}$")
    plt.ylabel("Test Statistic")
    plt.title("Median Test Statistics for (%s,%i) WIMPs"%(channel,mass))
    plt.tight_layout()
    plt.savefig(args.outfolder + "TS_" + channel + "_" + str(mass) + ".png")
    plt.close()
    
print("Finished plotting")
