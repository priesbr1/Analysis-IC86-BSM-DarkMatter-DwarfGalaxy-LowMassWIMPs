import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-d", "--data", type=str, default="/data/ana/BSM/IC86_LE_WIMP_dwarfgalaxy/ps_DRAGON/version-001-p00/IC86_201?_MC.npy",
                    dest="data", help="MC file(s) containing data to generate signal pdf from")
parser.add_argument("-p", "--spectra", type=str, default="/data/ana/BSM/IC86_LE_WIMP_dwarfgalaxy/annihilation_spectra/DM_CirelliSpectrum_dict_neutrinos_oscillated.npy",
                    dest="spectra", help="file containing WIMP anniliation spectra")
parser.add_argument("-c", "--channel", type=str, default="b",
                    dest="channel", help="annihilation channel to generate signal PDF for")
parser.add_argument("-o", "--output", type=str, default="/data/ana/BSM/IC86_LE_WIMP_dwarfgalaxy/pdfs/DRAGON_signal_pdfs.npy",
                    dest="outfile", help="file to save signal PDFs to (will also try to load previously-calculated PDFs)")
args = parser.parse_args()

os.system("""./submit_npx3_10GB_1CPU.sh python generate_signal_pdfs.py --data %s --spectra %s --channel %s --output %s"""%(args.data, args.spectra, args.channel, args.output))

print("Submitted")
