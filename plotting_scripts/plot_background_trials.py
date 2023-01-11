import numpy as np
import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
import argparse
import os, sys
import glob

parser = argparse.ArgumentParser()
parser.add_argument("-t", "--trials", type=str, default="/mnt/home/priesbr1/DM_Search/data/trials_results/trials_bkg_360deg_Jmax/bkg_trials_*_*.txt",
                    dest="trials", help="file(s) containing background trials to plot")
parser.add_argument("-o", "--outfolder", type=str, default="/mnt/scratch/priesbr1/DM_Search/Plots/background_trials/TS_360deg_Jmax/",
                    dest="outfolder", help="folder to save output plots to")
args = parser.parse_args()

bkg_trials = sorted(glob.glob(args.trials))
if isinstance(bkg_trials,str):
    bkg_trials = [bkg_trials]

if not os.path.isdir(args.outfolder):
    os.mkdir(args.outfolder)

for trial_set in bkg_trials:
    print("Plotting TS distribution from %s"%trial_set)
    channel_mass = trial_set[trial_set.find("bkg_trials")+11:-4]
    channel, mass = channel_mass.split("_")
    mass = int(mass)
    
    trials = np.genfromtxt(trial_set, delimiter="\t", dtype=None)
    ns = []
    TS = []
    for i in range(len(trials)):
        ns.append(trials[i][-2])
        TS.append(trials[i][-1])
    
    ns = np.array(ns)
    TS = np.array(TS)
    
    print("Median ns for %i (%s,%i) background trials: %i"%(len(ns), channel, mass, np.median(ns)))
    print("Median TS for %i (%s,%i) background trials: %.8f"%(len(TS), channel, mass, np.median(TS)))
    
    plt.figure()
    plt.hist(TS, bins=30)
    #plt.xscale("symlog", linthresh=10**(-4))
    plt.yscale("log")
    plt.xlabel("Test Statistic")
    plt.ylabel("Counts")
    plt.title("Background Test Statistic for (%s,%i) WIMPs"%(channel,mass))
    plt.tight_layout()
    plt.savefig(args.outfolder + "TS_" + channel + "_" + str(mass) + ".png")
    plt.close()

print("Finished plotting")
