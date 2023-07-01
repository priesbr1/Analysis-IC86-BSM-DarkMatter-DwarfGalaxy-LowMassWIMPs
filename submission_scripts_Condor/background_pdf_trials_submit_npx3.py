import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-j", "--num_jobs", type=int, default=500,
                    dest="num_jobs", help="number of jobs to submit")
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
parser.add_argument("-e", "--num_events", type=int, default=1000,
                    dest="num_events", help="number of background events to sample per source")
parser.add_argument("-n", "--num_trials", type=int, default=10,
                    dest="num_trials", help="number of background trials to run per core")
parser.add_argument("-x", "--max_trials", type=int, default=10000,
                    dest="max_trials", help="maximum number of trials per set")
parser.add_argument("-d", "--seed", type=int, default=None,
                    dest="seed", help="RNG seed for reproducibility")
parser.add_argument("-o", "--outfolder", type=str, default="/data/ana/BSM/IC86_LE_WIMP_dwarfgalaxy/trials_results/trials_bkg_360deg_Jmax_1kbe/",
                    dest="outfolder", help="folder to save background trials to")
args = parser.parse_args()

if (args.seed == None):
    for i in range(args.num_jobs):
        os.system("""./submit_npx3_10GB_2CPU.sh python background_pdf_trials.py --bkg_pdfs %s --sig_pdfs %s --sources %s --J_type %s --channel %s --mass %s --num_events %s --num_trials %s --max_trials %s --outfolder %s"""%(args.bkg_pdfs, args.sig_pdfs, args.sources, args.J_type, args.channel, args.mass, args.num_events, args.num_trials, args.max_trials, args.outfolder))
else:
    for i in range(args.num_jobs):
        os.system("""./submit_npx3_10GB_2CPU.sh python background_pdf_trials.py --bkg_pdfs %s --sig_pdfs %s --sources %s --J_type %s --channel %s --mass %s --num_events %s --num_trials %s --max_trials %s --seed %s --outfolder %s"""%(args.bkg_pdfs, args.sig_pdfs, args.sources, args.J_type, args.channel, args.mass, args.num_events, args.num_trials, args.max_trials, args.seed, args.outfolder))

print("Submitted")
