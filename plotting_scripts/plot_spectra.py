import numpy as np
import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
import argparse
import os, sys

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", type=str, default=None,
                    dest="filename", help="filename for spectra information")
parser.add_argument("-c", "--channels", type=str, nargs='+', default=["b", "Tau", "Mu", "W", "Nu"],
                    dest="channels", help="annihilation channels to use in plots")
parser.add_argument("-m", "--masses", type=int, nargs='+', default=[5, 10, 20, 30, 50, 80, 90, 100, 150, 200, 300, 500],
                    dest="masses", help="WIMP masses to use in plots")
parser.add_argument("-o", "--output", type=str, default=None,
                    dest="output", help="output folder for plots")
args = parser.parse_args()

if "neutrinos_e" in args.filename:
    observed = "NuE"
elif "neutrinos_mu" in args.filename:
    observed = "NuMu"
elif "neutrinos_tau" in args.filename:
    observed = "NuTau"
else:
    raise ValueError("Cannot find neutrino type from filename %s (looking for \"neutrinos_flavor\")"%filename)

out_folder = args.output + observed + "_spectra/"
if not os.path.isdir(out_folder):
    os.mkdir(out_folder)

data = np.load(args.filename, allow_pickle=True)
data = data.item()

annihilation_channels = data.keys()
WIMP_masses = data[list(annihilation_channels)[0]].keys()

for mass in args.masses:
    if mass not in WIMP_masses:
        print("Skipping mass %i -- mass not in WIMP masses: %s"%(mass, WIMP_masses))
        continue
    
    plt.figure(1) # LogX vs. EdNdE
    plt.figure(2) # E vs. dNdE
    
    for channel in sorted(args.channels):
        if channel not in annihilation_channels:
            print("Skipping channel %s -- channel not in keys: %s"%(channel, annihilation_channels))
            continue
        spectra = data[channel][mass]
        
        plt.figure(1)
        plt.plot(spectra["LogX"], spectra["EdNdE"], label=channel)
        
        plt.figure(2)
        plt.plot(spectra["Energy"], spectra["EdNdE"]/spectra["Energy"], label=channel)

    plt.figure(1)
    plt.xlabel("Log(E_particle/M_WIMP)")
    plt.ylabel("EdNdE")
    plt.title("%s DM Spectra, M_WIMP=%i GeV"%(observed, mass))
    plt.legend(loc="best")
    plt.tight_layout()
    plt.savefig(out_folder + "spectra_logX_WIMP_%i.png"%mass)

    plt.figure(2)
    plt.xlabel("Energy [GeV]")
    plt.xscale("log")
    plt.ylabel("dNdE [GeV^-1]")
    plt.title("%s DM Spectra, M_WIMP=%i GeV"%(observed, mass))
    plt.legend(loc="best")
    plt.tight_layout()
    plt.savefig(out_folder + "spectra_energy_WIMP_%i.png"%mass)

    plt.close(1)
    plt.close(2)
    print("Finished plots for WIMP mass = %i GeV"%mass)
