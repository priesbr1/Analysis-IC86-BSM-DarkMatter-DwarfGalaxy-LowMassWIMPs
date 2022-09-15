import numpy as np
import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-p", "--pdfs", type=str, default="/mnt/home/priesbr1/DM_Search/data/pdfs/DRAGON_background_pdfs_360deg.npy",
                    dest="pdfs", help="file containing PDFs to plot")
parser.add_argument("-b", "--background", type=int, default=1,
                    dest="background", help="whether or not PDFs are background (1) or signal (0)")
parser.add_argument("-o", "--output", type=str, default="/mnt/scratch/priesbr1/DM_Search/Plots/PDFs/",
                    dest="output", help="output folder for plots")
args = parser.parse_args()

pdfs = np.load(args.pdfs, allow_pickle=True)
pdfs = pdfs.item()

if (args.background == True):
    filename_start = args.pdfs.find("background")
    band_size_idx = [char.isdigit() for char in args.pdfs[filename_start:]].index(True)
    band_size = args.pdfs[filename_start:][band_size_idx:-4]
    print("Declination band size:", band_size)
    for source in pdfs.keys():
        for channel in pdfs[source].keys():
            for mass in pdfs[source][channel].keys():
                
                counts, bins = pdfs[source][channel][mass]["Counts"], pdfs[source][channel][mass]["Bins"]
                
                plt.figure()
                plt.hist(bins[:-1], bins=bins, weights=counts/np.sum(counts))
                plt.xlabel(r"Angle from source [$\mathrm{rad}$]")
                plt.ylabel("Probability Density")
                plt.title("%s Background PDF for (%s,%i) WIMPs"%(source, channel, mass))
                plt.savefig(args.output + "background_pdf_" + source + "_" + channel + "_" + str(mass) + "_" + band_size  + ".png")
                plt.close()
        
        print("Finished plotting background PDFs for %s"%source)

else:
    for channel in pdfs.keys():
        for mass in pdfs[channel].keys():
            
            counts, bins = pdfs[channel][mass]["Counts"], pdfs[channel][mass]["Bins"]
            
            plt.figure()
            plt.hist(bins[:-1], bins=bins, weights=counts/np.sum(counts))
            plt.xlabel(r"Energy-weighted angular separation [$\mathrm{rad}$]")
            plt.ylabel("Probability Density")
            plt.title("Signal PDF for (%s,%i) WIMPs"%(channel, mass))
            plt.savefig(args.output + "signal_pdf_" + channel + "_" + str(mass) + ".png")
            plt.close()
        
        print("Finished plotting signal PDFs for %s annihilation channel"%channel)
