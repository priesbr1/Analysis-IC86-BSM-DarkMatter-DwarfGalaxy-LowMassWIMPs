import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-d", "--data", type=str, default="/data/ana/BSM/IC86_LE_WIMP_dwarfgalaxy/ps_DRAGON/version-001-p00/IC86_201?_exp.npy",
                    dest="data", help="file(s) containing data to generate background pdf from")
parser.add_argument("-s", "--sources", type=str, default="/data/ana/BSM_IC86_LE_WIMP_dwarfgalaxy/analysis_sources_ra_dec_jfactors.txt",
                    dest="sources", help="file containing RA, dec, and J-factor information for each source")
parser.add_argument("-c", "--channel", type=str, default="b",
                    dest="channel", help="annihilation channel to generate background PDFs for")
parser.add_argument("-m", "--mass", type=int, default=10,
                    dest="mass", help="WIMP mass to generate background PDFs for")
parser.add_argument("-b", "--band_size", type=int, default=360,
                    dest="band_size", help="size of declination band in degrees")
parser.add_argument("-e", "--seed", type=int, default=None,
                    dest="seed", help="RNG seed for reproducibility")
parser.add_argument("-o", "--output", type=str, default="/data/ana/BSM/IC86_LE_WIMP_dwarfgalaxy/pdfs/DRAGON_background_pdfs_360deg.npy",
                    dest="outfile", help="file to save background PDFs to (will also try to load previously-calculated PDFs)")
args = parser.parse_args()

if (args.seed == None):
    os.system("""./submit_npx_10GB_1CPU.sh python generate_background_pdfs.py --data %s --sources %s --channel %s --mass %s --band_size %s --output %s"""%(args.data, args.sources, args.channel, args.mass, args.band_size, args.output))
else:
    os.system("""./submit_npx_10GB_1CPU.sh python generate_background_pdfs.py --data %s --sources %s --channel %s --mass %s --band_size %s --seed %s --output %s"""%(args.data, args.sources, args.channel, args.mass, args.band_size, args.seed, args.output))

print("Submitted")
