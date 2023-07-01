import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-d", "--datafolder", type=str, default="/data/ana/BSM/IC86_LE_WIMP_dwarfgalaxy/energy_optimization_MC/",
                    dest="datafolder", help="folder containing i3 MC datafiles")
parser.add_argument("-s", "--spectra", type=str, default="/data/ana/BSM/IC86_LE_WIMP_dwarfgalaxy/annihilation_spectra/DM_CirelliSpectrum_dictNoEW_neutrinos_oscillated.npy",
                    dest="spectra", help="file containing WIMP annihilation spectra")
parser.add_argument("-c", "--channel", type=str, default="b",
                    dest="channel", help="annihilation channel to optimize energy range for")
parser.add_argument("-m", "--mass", type=int, default=10,
                    dest="mass", help="WIMP mass to optimize energy range for")
parser.add_argument("-o", "--outfolder", type=str, default="/data/ana/BSM/IC86_LE_WIMP_dwarfgalaxy/energy_optimization_results/",
                    dest="outfolder", help="output directory for energy range optimization results")
args = parser.parse_args()


for part in range(10):
    os.system("""./submit_npx3_10GB_1CPU.sh python energy_range_optimization.py --datafolder %s --spectra %s --channel %s --mass %s --part %s --outfolder %s"""%(args.datafolder, args.spectra, args.channel, args.mass, part, args.outfolder))

print("Submitted")
