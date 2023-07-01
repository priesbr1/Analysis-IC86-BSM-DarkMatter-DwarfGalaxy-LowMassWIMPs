import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-j", "--num_jobs", type=int, default-300,
                    dest="num_jobs", help="number of jobs to submit"
parser.add_argument("-b", "--bkg_pdfs", type=str, default="/data/ana/BSM/IC86_LE_WIMP_dwarfgalaxy/pdfs/DRAGON_background_pdfs_360deg.npy",
                    dest="bkg_pdfs", help="file containing background PDFs for each mass, channel, and source")
parser.add_argument("-p", "--sig_pdfs", type=str, default="/data/ana/BSM/IC86_LE_WIMP_dwarfgalaxy/pdfs/DRAGON_signal_pdfs.npy",
                    dest="sig_pdfs", help="file containing signal PDFs for each mass and channel")
parser.add_argument("-s", "--sources", type=str, default="/data/ana/BSM/IC86_LE_WIMP_dwarfgalaxy/analysis_sources_ra_dec_jfactors.txt",
                    dest="sources", help="file containing RA, dec, and J-factor information for each source")
parser.add_argument("-j", "--J_type", type=str, choices=["0.1","0.2","0.5","10","min","max"], default="max",
                    dest="J_type", help="J-factor values to use based on opening half-angles")
parser.add_argument("-c", "--channel", type=str, default="b",
                    dest="channel", help="annhiliation channel to run background trials for")
parser.add_argument("-m", "--mass", type=int, default=10,
                    dest="mass", help="WIMP mass to rub background trials for")
parser.add_argument("-n", "--num_trials", type=int, default=1,
                    dest="num_trials", help="number of signal trials to run per core")
parser.add_argument("-d", "--seed", type=int, default=None,
                    dest="seed", help="RNG seed for reproducibility")
parser.add_argument("-o", "--outfolder", type=str, default="/data/ana/BSM/IC86_LE_WIMP_dwarfgalaxy/trials_results/trials_sig_360deg_Jmax_distributed/",
                    dest="outfolder", help="folder to save signal trials to")
parser.add_argument("-f", "--bkg_folder", type=str, default="/data/ana/BSM/IC86_LE_WIMP_dwarfgalaxy/trials_results/trials_bkg_360deg_Jmax/",
                    dest="bkg_folder", help="folder containing results from corresponding background trials")
args = parser.parse_args()

if (args.seed == None):
    for i in range(args.num_jobs):
        os.system("""./submit_npx3_10GB_1CPU.sh python signal_pdf_trials.py --bkg_pdfs %s --sig_pdfs %s --sources %s --J_type %s --channel %s --mass %s --num_trials %s --outfolder %s --bkg_folder %s"""%(args.bkg_pdfs, args.sig_pdfs, args.sources, args.J_type, args.channel, args.mass, args.num_trials, args.outfolder, args.bkg_folder))
else:
    for i in range(args.num_jobs):
        os.system("""./submit_npx3_10GB_1CPU.sh python signal_pdf_trials.py --bkg_pdfs %s --sig_pdfs %s --sources %s --J_type %s --channel %s --mass %s --num_trials %s --seed %s --outfolder %s --bkg_folder %s"""%(args.bkg_pdfs, args.sig_pdfs, args.sources, args.J_type, args.channel, args.mass, args.num_trials, args.seed, args.outfolder, args.bkg_folder))

print("Submitted")
