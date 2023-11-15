import numpy as np
import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
import argparse
import os, sys
import glob

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", type=str, default="/mnt/home/priesbr1/DM_Search/data/cross_section_results/cross_section_results_360deg_Jmax_distributed_erange.npy",
                    dest="filename", help="filename with cross section results")
parser.add_argument("-c", "--comparison", type=str, default="/mnt/home/priesbr1/DM_Search/data/comparison_data/",
                    dest="comparison", help="folder containing cross seciton results to compare to")
parser.add_argument("-o", "--output", type=str, default="/mnt/scratch/priesbr1/DM_Search/Plots/cross_sections/",
                    dest="output", help="output folder for plots")
args = parser.parse_args()

params = args.filename[args.filename.rfind("/")+len("cross_section_results")+2:args.filename.rfind(".")]

cross_sections = np.load(args.filename, allow_pickle=True)
cross_sections = cross_sections.item()

legends = {"b":r"$b\bar{b}$", "Mu":r"$\mu^+\mu^-$", "Tau":r"$\tau^+\tau^-$", "W":r"$W^+W^-$", "Nu":r"$\nu\bar{\nu}$", "Nue":r"$\nu_e\bar{\nu}_e$", "NuMu":r"$\nu_{\mu}\bar{\nu}_{\mu}$",
           "NuTau":r"$\nu_{\tau}\bar{\nu}_{\tau}$"}
colors = {"b":"mediumturquoise", "Mu":"forestgreen", "Tau":"dodgerblue", "W":"darkmagenta", "Nu":"mediumvioletred", "Nue":"darkorchid", "NuMu":"orchid", "NuTau":"lightpink"}

fig = plt.figure()
ax = fig.add_subplot(111)
fig.subplots_adjust(top=0.85)
fig.patch.set_facecolor("white")

for channel in cross_sections.keys():
    m_WIMP = []
    sigma_v = []
    for mass in cross_sections[channel].keys():
        m_WIMP.append(mass)
        sigma_v.append(cross_sections[channel][mass])
    
    ax.plot(m_WIMP, sigma_v, color=colors[channel], linestyle="-", label=legends[channel], linewidth=2)
ax.semilogx()
ax.semilogy()
ax.set_xlabel(r"$m_{WIMP}$ [GeV/c$^{2}$]", fontsize=12)
if ("unblind" in args.filename) or ("upper_limit" in args.filename):
    ax.set_ylabel(r"$\langle \sigma v \rangle$ Limits [cm$^{3}$ s$^{-1}$]", fontsize=12)
else:
    ax.set_ylabel(r"$\langle \sigma v \rangle$ Sensitivity [cm$^{3}$ s$^{-1}$]", fontsize=12)
ax.legend(loc="best", ncol=2, prop={"size":12}, labelspacing=0.5)
ax.grid(color="grey", linestyle="-", linewidth=0.5, alpha=0.5)
plt.tight_layout()
plt.savefig(args.output + "cross_section_results_" + params + ".png")
print("Finished solo plot")

comps = dict()
comp_files = sorted(glob.glob(args.comparison+"*.csv"))
for cf in comp_files:
    print("Using comparison", cf)
    channel_result = cf[cf.rfind("/")+1:cf.rfind(".")]
    channel, detector, year, location, confidence = channel_result.split("_")
    channel = channel[len(channel)//2:]
    confidence = confidence[2:]
    
    comp_results = np.genfromtxt(cf, delimiter=",", dtype=None)
    if detector+"_"+year+"_"+location+"_"+confidence not in comps.keys():
        comps[detector+"_"+year+"_"+location+"_"+confidence] = dict()
    comps[detector+"_"+year+"_"+location+"_"+confidence][channel] = dict()
    comps[detector+"_"+year+"_"+location+"_"+confidence][channel]["m_WIMP"] = np.array([comp_results[i][0] for i in range(len(comp_results))])
    comps[detector+"_"+year+"_"+location+"_"+confidence][channel]["sigma_v"] = np.array([comp_results[i][1] for i in range(len(comp_results))])

comp_linestyles = ["--", ":", "-."]
if len(comps.keys()) > len(comp_linestyles):
    raise RuntimeError("Too many comparisons -- cannot plot with separate linestyles")

fig = plt.figure(figsize=(10,6))
ax = fig.add_subplot(111)
fig.subplots_adjust(top=0.85)
fig.patch.set_facecolor("white")

all_channels = []
all_channels.extend(list(cross_sections.keys()))
for result in comps.keys():
    all_channels.extend(list(comps[result].keys()))
all_channels = list(dict.fromkeys(all_channels).keys()) # Preserves order of insertion

for channel in all_channels:
    ax.plot([], [], color=colors[channel], linestyle="-", label=legends[channel], linewidth=2)

for channel in cross_sections.keys():
    m_WIMP = []
    sigma_v = []
    for mass in cross_sections[channel].keys():
        m_WIMP.append(mass)
        sigma_v.append(cross_sections[channel][mass])
    
    ax.plot(m_WIMP, sigma_v, color=colors[channel], linestyle="-", linewidth=2)
if ("unblind" in args.filename) or ("upper_limit" in args.filename):
    ax.plot([], [], color="black", linestyle="-", label="Current Limits (29DG, 90% CL)", linewidth=2)
else:
    ax.plot([], [], color="black", linestyle="-", label="Current Sensitivities (29DG, 90% CL)", linewidth=2)

for i, result in enumerate(comps.keys()):
    detector, year, location, confidence = result.split("_")
    ax.plot([], [], color="black", linestyle=comp_linestyles[i], label=detector+" "+year+" ("+location+", "+confidence+"% CL)", linewidth=2)
    for channel in comps[result].keys():
        ax.plot(comps[result][channel]["m_WIMP"], comps[result][channel]["sigma_v"], color=colors[channel], linestyle=comp_linestyles[i], linewidth=2)

ax.semilogx()
ax.semilogy()
ax.set_xlabel(r"$m_{WIMP}$ [GeV/c$^{2}$]", fontsize=12)
ax.set_ylabel(r"$\langle \sigma v \rangle$ [cm$^{3}$ s$^{-1}$]", fontsize=12)
ax.legend(bbox_to_anchor=(1.04,0.5), loc="center left", prop={"size":12}, labelspacing=0.5)
ax.grid(color="grey", linestyle="-", linewidth=0.5, alpha=0.5)
plt.tight_layout()
plt.savefig(args.output + "cross_section_results_comparison_" + params + ".png")
print("Finished comparison plot")
