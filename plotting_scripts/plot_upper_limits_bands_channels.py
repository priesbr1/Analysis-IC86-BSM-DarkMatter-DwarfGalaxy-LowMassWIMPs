import numpy as np
import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
import argparse
import os, sys
import glob

parser = argparse.ArgumentParser()
parser.add_argument("-u", "--upper_limits", type=str, default="/mnt/home/priesbr1/DM_Search/data/cross_section_results/cross_section_results_upper_limits_360deg_Jmax_distributed_erange.npy",
                    dest="limits", help="filename with upper limits")
parser.add_argument("-b", "--bands", type=str, default="/mnt/home/priesbr1/DM_Search/data/cross_section_results/cross_section_results_bands_360deg_Jmax_distributed_erange.npy",
                    dest="bands", help="filename with Brazil bands")
parser.add_argument("-c", "--comparison", type=str, default="/mnt/home/priesbr1/DM_Search/data/comparison_data_IceCube/",
                    dest="comparison", help="folder containing cross section results to compare to")
parser.add_argument("-n", "--channel", type=str, default="b",
                    dest="channel", help="channel to plot")
parser.add_argument("-o", "--output", type=str, default="/mnt/scratch/priesbr1/DM_Search/Plots/cross_sections/",
                    dest="output", help="output folder for plots")
args = parser.parse_args()

supported_channels = ["b", "Mu", "Tau", "W", "Nu", "Nue", "NuMu", "NuTau"]
if args.channel not in supported_channels:
    raise ValueError("Channel (%s) not supported"%args.channel)

if "bands" not in args.bands:
    raise RuntimeError("Cross sections file does not have bands")

params = args.limits[args.limits.rfind("/")+len("cross_section_results_upper_limits")+2:args.limits.rfind(".")]

upper_limits = np.load(args.limits, allow_pickle=True)
upper_limits = upper_limits.item()

bands = np.load(args.bands, allow_pickle=True)
bands = bands.item()

legends = {"b":r"$b\bar{b}$", "Mu":r"$\mu^+\mu^-$", "Tau":r"$\tau^+\tau^-$", "W":r"$W^+W^-$", "Nu":r"$\nu\bar{\nu}$", "Nue":r"$\nu_e\bar{\nu}_e$", "NuMu":r"$\nu_{\mu}\bar{\nu}_{\mu}$",
           "NuTau":r"$\nu_{\tau}\bar{\nu}_{\tau}$"}
colors = ["tab:blue", "tab:orange", "tab:green", "tab:red", "tab:purple", "tab:cyan", "tab:olive", "tab:pink", "tab:brown"]
linestyles = {"b":"-", "Mu":"-", "Tau":"-", "W":"-", "Nu":"-", "Nue":"--", "NuMu":":", "NuTau":"-."}
nu_channels = ["Nu", "Nue", "NuMu", "NuTau"]

fig = plt.figure(figsize=(10,6))
ax = fig.add_subplot(111)
fig.subplots_adjust(top=0.85)
fig.patch.set_facecolor("white")

m_WIMP_lims = []
m_WIMP_bands = []
sigma_v_lims = []
sigma_v_bands = []
for mass in upper_limits[args.channel].keys():
    m_WIMP_lims.append(mass)
    sigma_v_lims.append(upper_limits[args.channel][mass])
    m_WIMP_bands.append(mass)
    sigma_v_bands.append([])
    for sigma in bands[args.channel][mass]["Cross Section"].keys():
        sigma_v_bands[-1].append(bands[args.channel][mass]["Cross Section"][sigma])
    sigma_v_bands[-1] = sorted(sigma_v_bands[-1])

ax.plot(m_WIMP_bands, [sigma_v_bands[i][2] for i in range(len(sigma_v_bands))], color="black", linestyle="-", linewidth=2)
ax.fill_between(m_WIMP_bands, [sigma_v_bands[i][1] for i in range(len(sigma_v_bands))], [sigma_v_bands[i][3] for i in range(len(sigma_v_bands))], color="tab:gray", alpha=0.5, linestyle="-", linewidth=2)
ax.fill_between(m_WIMP_bands, [sigma_v_bands[i][0] for i in range(len(sigma_v_bands))], [sigma_v_bands[i][4] for i in range(len(sigma_v_bands))], color="tab:gray", alpha=0.5, linestyle="-", linewidth=2)
ax.plot(m_WIMP_lims, sigma_v_lims, color="black", linestyle="--", linewidth=2)
ax.plot([], [], color="black", linestyle="--", label="Current %s Limits (29DG, 90%% CL)"%legends[args.channel], linewidth=2)
ax.plot([], [], color="black", linestyle="-", label="Current %s Sensitivities (29DG, 90%% CL)"%legends[args.channel], linewidth=2)
ax.semilogx()
ax.semilogy()
ax.set_xlabel(r"$m_{WIMP}$ [GeV/c$^{2}$]", fontsize=12)
ax.set_ylabel(r"$\langle \sigma v \rangle$ [cm$^{3}$ s$^{-1}$]", fontsize=12)
ax.legend(loc="best", prop={"size":12}, labelspacing=0.5)
ax.grid(color="grey", linestyle="-", linewidth=0.5, alpha=0.5)
plt.tight_layout()
plt.savefig(args.output + "cross_section_results_upper_limits_bands_%s_"%args.channel + params + ".png")
print("Finished solo plot")

comps = dict()
comp_files = sorted(glob.glob(args.comparison+"*.csv"))
for cf in comp_files:
    channel_result = cf[cf.rfind("/")+1:cf.rfind(".")]
    channel, detector, year, location, confidence = channel_result.split("_")
    channel = channel[len(channel)//2:]
    confidence = confidence[2:]
    
    if (channel == args.channel) or ((args.channel in nu_channels) and (args.channel in channel)):
        print("Using comparison", cf)
        comp_results = np.genfromtxt(cf, delimiter=",", dtype=None)
        if detector+"_"+year+"_"+location+"_"+confidence not in comps.keys():
            comps[detector+"_"+year+"_"+location+"_"+confidence] = dict()
        comps[detector+"_"+year+"_"+location+"_"+confidence][channel] = dict()
        comps[detector+"_"+year+"_"+location+"_"+confidence][channel]["m_WIMP"] = np.array([comp_results[i][0] for i in range(len(comp_results))])
        comps[detector+"_"+year+"_"+location+"_"+confidence][channel]["sigma_v"] = np.array([comp_results[i][1] for i in range(len(comp_results))])

if (len(comps.keys()) > len(colors)):
    raise RuntimeError("Too many comparisons -- cannot plot with separate colors")

fig = plt.figure(figsize=(10,6))
ax = fig.add_subplot(111)
fig.subplots_adjust(top=0.85)
fig.patch.set_facecolor("white")

if (args.channel in nu_channels):
    for channel in nu_channels:
        ax.plot([], [], color="tab:gray", linestyle=linestyles[channel], label=legends[channel], linewidth=2)
else:
    ax.plot([], [], color="tab:gray", linestyle=linestyles[args.channel], label=legends[args.channel], linewidth=2)

m_WIMP_lims = []
m_WIMP_bands = []
sigma_v_lims = []
sigma_v_bands = []
for mass in upper_limits[args.channel].keys():
    m_WIMP_lims.append(mass)
    sigma_v_lims.append(upper_limits[args.channel][mass])
    m_WIMP_bands.append(mass)
    sigma_v_bands.append([])
    for sigma in bands[args.channel][mass]["Cross Section"]:
        sigma_v_bands[-1].append(bands[args.channel][mass]["Cross Section"][sigma])
    sigma_v_bands[-1] = sorted(sigma_v_bands[-1])
    
ax.plot(m_WIMP_bands, [sigma_v_bands[i][2] for i in range(len(sigma_v_bands))], color="black", linestyle="-", linewidth=2)
ax.fill_between(m_WIMP_bands, [sigma_v_bands[i][1] for i in range(len(sigma_v_bands))], [sigma_v_bands[i][3] for i in range(len(sigma_v_bands))], color="gray", alpha=0.5, linestyle="-", linewidth=2)
ax.fill_between(m_WIMP_bands, [sigma_v_bands[i][0] for i in range(len(sigma_v_bands))], [sigma_v_bands[i][4] for i in range(len(sigma_v_bands))], color="gray", alpha=0.5, linestyle="-", linewidth=2)
ax.plot(m_WIMP_lims, sigma_v_lims, color="black", linestyle="--", linewidth=2)
ax.plot([], [], color="black", linestyle="--", label="Current Limits (29DG, 90% CL)", linewidth=2)
ax.plot([], [], color="black", linestyle="-", label="Current Sensitivities (29DG, 90% CL)", linewidth=2)

for i, result in enumerate(comps.keys()):
    detector, year, location, confidence = result.split("_")
    ax.plot([], [], color=colors[i], linestyle="-", label=detector+" "+year+" ("+location+", "+confidence+"% CL)", linewidth=2)
    for channel in comps[result].keys():
        ax.plot(comps[result][channel]["m_WIMP"], comps[result][channel]["sigma_v"], color=colors[i], linestyle=linestyles[channel], linewidth=2)

ax.semilogx()
ax.semilogy()
ax.set_xlabel(r"$m_{WIMP}$ [GeV/c$^{2}$]", fontsize=12)
ax.set_ylabel(r"$\langle \sigma v \rangle$ [cm$^{3}$ s$^{-1}$]", fontsize=12)
ax.legend(bbox_to_anchor=(1.04,0.5), loc="center left", prop={"size":12}, labelspacing=0.5)
ax.grid(color="grey", linestyle="-", linewidth=0.5, alpha=0.5)
plt.tight_layout()
plt.savefig(args.output + "cross_section_results_upper_limits_bands_comparison_%s_"%args.channel + params + ".png")
print("Finished comparison plot")
