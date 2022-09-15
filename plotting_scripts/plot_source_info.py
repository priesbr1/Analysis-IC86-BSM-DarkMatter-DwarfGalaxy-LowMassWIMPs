import numpy as np
import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
import argparse
import healpy as hp

parser = argparse.ArgumentParser()
parser.add_argument("-s", "--sources", type=str, default="/mnt/home/priesbr1/DM_Search/data/analysis_sources_ra_dec_jfactors.txt",
                    dest="sources", help="file containing source information")
parser.add_argument("-o", "--outfolder", type=str, default="/mnt/scratch/priesbr1/DM_Search/Plots/source_info/",
                    dest="outfolder", help="output folder to save plots to")
args = parser.parse_args()

# Source encoding - 0:name; 1:RA; 2:Dec; 3,4,5:r_h; 6,7,8:J0.1; 9,10,11:J0.2; 12,13,14:J0.5; 15,16,17:J10
sources = np.genfromtxt(args.sources, delimiter=",", dtype=("U20",float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float))
names, ras, decs, rhs, J0_1, J0_2, J0_5, J10 = [], [], [], [], [], [], [], []

for i in range(len(sources)):
    names.append(sources[i][0])
    ras.append(np.radians(sources[i][1]))
    decs.append(np.radians(sources[i][2]))
    
    if not np.isnan(sources[i][4]):
        rhs.append([sources[i][3],sources[i][4],sources[i][5]])
    else:
        rhs.append([sources[i][3],None,None])
    
    if not np.isnan(sources[i][6]):
        J0_1.append([sources[i][6],sources[i][7],sources[i][8]])
    else:
        J0_1.append([None,None,None])
    
    if not np.isnan(sources[i][9]):
        J0_2.append([sources[i][9],sources[i][10],sources[i][11]])
    else:
        J0_2.append([None,None,None])
    
    if not np.isnan(sources[i][12]):
        J0_5.append([sources[i][12],sources[i][13],sources[i][14]])
    else:
        J0_5.append([None,None,None])
    
    if not np.isnan(sources[i][15]):
        J10.append([sources[i][15],sources[i][16],sources[i][17]])
    else:
        J10.append([None,None,None])
    
ras = np.array(ras)
decs = np.array(decs)
rhs = np.array(rhs)
J0_1 = np.array(J0_1)
J0_2 = np.array(J0_2)
J0_5 = np.array(J0_5)
J10 = np.array(J10)

for i in range(len(names)):
    print("%s: (%.2f, %.2f)"%(names[i], ras[i]*180./np.pi, decs[i]*180./np.pi))

nside = 32
npix = hp.nside2npix(nside)
hpx_map = np.zeros(npix, dtype=int)
J_factor_sets = [J0_1, J0_2, J0_5, J10]
J_factor_colors = ["blue", "orange", "teal", "purple"]
J_factor_extensions = [0.1, 0.2, 0.5, 10]
max_extensions = []

fig = plt.figure(num=1, figsize=(11,6))
hp.mollview(hpx_map, fig=1, coord="C", rot=180, title="Dwarf Galaxy Sources", norm="hist", cbar=False)
ax1 = fig.get_axes()[0]
ax1.text(2.05,0., r"$0^\circ$", ha="left", va="center")
ax1.text(-2.05,0., r"$360^\circ$", ha="right", va="center")
ax1.text(1.9,0.40, r"$30^\circ$", ha="left", va="center")
ax1.text(1.4,0.775, r"$60^\circ$", ha="left", va="center")
ax1.text(1.9,-0.40, r"$-30^\circ$", ha="left", va="center")
ax1.text(1.4,-0.775, r"$-60^\circ$", ha="left", va="center")
hp.graticule()
hp.projscatter(np.pi/2-decs, ras, color="red", marker="+", s=75, label="Dwarf Galaxies")
plt.legend(loc="best")
plt.savefig(args.outfolder + "source_skymap.png", bbox_inches="tight")
plt.close(1)

fig = plt.figure(figsize=(12,8))
for i in range(len(names)):
    if rhs[i][1] is not None:
        plt.errorbar(i+0.5, rhs[i][0], yerr=np.array([rhs[i][1],rhs[i][2]]).reshape(2,1), color="black", marker="o",
                     markersize=5, capsize=2)
        max_extensions.append(rhs[i][0]+rhs[i][1])
    else:
        plt.errorbar(i+0.5, rhs[i][0], color="black", marker="o", markersize=5)
        max_extensions.append(rhs[i][0])
plt.errorbar([], [], color="black", marker="o", markersize=5, capsize=2, label="Angular Radius")
for j, theta in enumerate(np.array(J_factor_extensions)[J_factor_extensions < 1.5/60*np.amax([max_extensions])]):
    plt.axhline(theta*60, color=J_factor_colors[j], linestyle="dashed", label=r"J(%s$^{\circ}$)"%theta)
plt.xticks(np.arange(len(names))+0.5, names, rotation=90)
plt.ylabel("Angular Size [arcmin]")
plt.legend(loc="best")
plt.tight_layout()
plt.savefig(args.outfolder + "source_sizes.png")
plt.close()

fig = plt.figure(figsize=(12,8))
for i in range(len(names)):
    for j, J_factor in enumerate(J_factor_sets):
        if J_factor[i][0] is not None:
            plt.errorbar(i+0.2*(j+1), J_factor[i][0], yerr=np.array([J_factor[i][1],J_factor[i][2]]).reshape(2,1),
                         color=J_factor_colors[j], marker="o", markersize=5, capsize=2)
        if (i % 2 == 1):
            plt.axvspan(i,i+1,color="gray",alpha=0.1)
for j, color in enumerate(J_factor_colors):
    plt.errorbar([], [], color=color, marker="o", markersize=5, capsize=2, label=r"J(%s$^{\circ}$)"%J_factor_extensions[j])
    
plt.xticks(np.arange(len(names))+0.5, names, rotation=90)
plt.ylabel(r"$\log_{10}(J) \left[\mathrm{GeV^{2}/cm^{5}}\right]$")
plt.legend(loc="best")
plt.tight_layout()
plt.savefig(args.outfolder + "source_Jfactors.png")
plt.close()
