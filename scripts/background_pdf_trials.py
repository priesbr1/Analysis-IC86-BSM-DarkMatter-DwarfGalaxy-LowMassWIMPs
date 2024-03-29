import numpy as np
import argparse
import os,time
from multiprocessing import Process

def generate_rand_from_pdf(pdf, x_grid, n=1000):
    cdf = np.cumsum(pdf)
    cdf = cdf / cdf[-1]
    values = np.random.rand(n)
    value_bins = np.searchsorted(cdf, values)
    random_from_cdf = x_grid[value_bins]
    return random_from_cdf

def calculate_angle(dec_source, ra_source, dec_nu, ra_nu):
    return np.arccos(np.cos(dec_source)*np.cos(dec_nu)*np.cos(ra_source - ra_nu) + np.sin(dec_source)*np.sin(dec_nu))

parser = argparse.ArgumentParser()
parser.add_argument("-b", "--bkg_pdfs", type=str, default="/mnt/home/priesbr1/DM_Search/data/pdfs/DRAGON_background_pdfs_360deg.npy",
                    dest="bkg_pdfs", help="file containing background PDFs for each mass, channel, and source")
parser.add_argument("-p", "--sig_pdfs", type=str, default="/mnt/home/priesbr1/DM_Search/data/pdfs/DRAGON_signal_pdfs.npy",
                    dest="sig_pdfs", help="file containing signal PDFs for each mass and channel")
parser.add_argument("-s", "--sources", type=str, default="/mnt/home/priesbr1/DM_Search/data/analysis_sources_ra_dec_jfactors.txt",
                    dest="sources", help="file containing RA, dec, and J-factor information for each source")
parser.add_argument("-j", "--J_type", type=str, choices=["0.1","0.2","0.5","10","min","max"], default="max",
                    dest="J_type", help="J-factor values to use based on opening half-angle")
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
parser.add_argument("-o", "--outfolder", type=str, default="/mnt/home/priesbr1/DM_Search/data/trials_results/DRAGON_background_trials/",
                    dest="outfolder", help="folder to save background trials to")
args = parser.parse_args()

if (args.seed != None):
    np.random.seed(args.seed)

bkg_pdfs = np.load(args.bkg_pdfs, allow_pickle=True)
bkg_pdfs = bkg_pdfs.item()
sig_pdfs = np.load(args.sig_pdfs, allow_pickle=True)
sig_pdfs = sig_pdfs.item()

def extract_source_info(sources, J_type, J_indices_map):
    names = []
    ra = np.array([])
    dec = np.array([])
    J_factors = np.array([])
    
    J_indices = J_indices_map[J_type]
    
    for i in range(len(sources)):
        for idx in J_indices:
            if not np.isnan(sources[i][idx]):
                names.append(sources[i][0])
                ra = np.append(ra, sources[i][1]*np.pi/180.)
                dec = np.append(dec, sources[i][2]*np.pi/180.)
                J_factors = np.append(J_factors, sources[i][idx])
                break
    
    J_factors = 10**J_factors
    
    return names, ra, dec, J_factors

J_indices_map = {"0.1": [6], "0.2": [9], "0.5": [12], "10": [15], "min": [6,9,12,15], "max": [15,12,9,6]}

# Source encoding - 0:name; 1:RA; 2:Dec; 3,4,5:r_h; 6,7,8:J0.1; 9,10,11:J0.2; 12,13,14:J0.5; 15,16,17:J10
sources = np.genfromtxt(args.sources, delimiter=",", dtype=("U20",float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float))
names, ra, dec, J_factors = extract_source_info(sources, args.J_type, J_indices_map)
#print(J_factors/np.sum(J_factors))

outfile = "bkg_trials_" + args.channel + "_" + str(args.mass) + ".txt"

def run_bkg_trials(bkg_pdfs, sig_pdfs, names, J_factors, args, outfile, num_bkg_trials, max_trials):
    t1 = time.time()
    
    for t in range(num_bkg_trials):
        t1_trial = time.time()
        data = []
        for source in names:
            counts, bins = bkg_pdfs[source][args.channel][args.mass]["Counts"], bkg_pdfs[source][args.channel][args.mass]["Bins"]
            data.append(generate_rand_from_pdf(counts,bins[:-1],args.num_events))
        
        # Fit
        N = int(np.sum([len(samp) for samp in data])/len(names))
        llhs = []
        
        for j, source in enumerate(names):
            source_llh = []
            
            bkg_bins = np.searchsorted(bkg_pdfs[source][args.channel][args.mass]["Bins"], data[j], "right")-1 # Find left edges of bins
            bkg_probs = np.array(bkg_pdfs[source][args.channel][args.mass]["Counts"][bkg_bins], dtype=float) # Get probabilities for bins
            bkg_probs /= np.sum(bkg_pdfs[source][args.channel][args.mass]["Counts"]) # Normalize probabilities
            sig_bins = np.searchsorted(sig_pdfs[args.channel][args.mass]["Bins"], data[j], "right")-1
            sig_probs = np.array(sig_pdfs[args.channel][args.mass]["Counts"][sig_bins], dtype=float)
            sig_probs /= np.sum(sig_pdfs[args.channel][args.mass]["Counts"])
            
            for ns in range(int(N/100)):
                ns = float(ns)
                combs = (ns/N)*sig_probs + (1-(ns/N))*bkg_probs
                source_llh.append(-1*np.sum(np.log(combs)))
            
            llhs.append(source_llh)
        
        bfs = np.argmin(llhs, axis=1)
        tss = []
        for L,bf in zip(llhs,bfs):
            tss.append(2*(L[0] - L[bf]))
        
        bf = int(np.round(np.sum(bfs * J_factors/np.sum(J_factors))))
        #print("Best-fit signal events:", bf)
        ts = np.sum(np.array(tss) * J_factors/np.sum(J_factors))
        #print("TS:", ts)
        
        t2_trial = time.time()
        #print("Trial time: %.2f s"%(t2_trial-t1_trial))
        
        written = False
        while (written == False):
            if not os.path.isfile(args.outfolder+outfile):
                f = open(args.outfolder+outfile,"w")
                f.write("# ")
                for j, source in enumerate(names):
                    f.write("%s_ns \t %s_TS \t "%(source, source))
                f.write("ns \t TS \n")
                f.close()
            else:
                try:
                    os.rename(args.outfolder+outfile,args.outfolder+outfile)
                    num_complete = sum(1 for line in open(args.outfolder+outfile)) - 1
                    if (num_complete < max_trials):
                        f = open(args.outfolder+outfile,"a")
                        for j, source in enumerate(names):
                            f.write("%s \t %s \t "%(bfs[j], tss[j]))
                        f.write("%s \t %s \n"%(bf,ts))
                        f.close()
                    else:
                        print("Desired number of trials reached")
                        return
                    written = True
                except:
                    print("Waiting for file access...")
                    time.sleep(1)
    
    t2 = time.time()
    print("Time: %.2f sec"%(t2-t1))

p1 = Process(target=run_bkg_trials, args=(bkg_pdfs, sig_pdfs, names, J_factors, args, outfile, args.num_trials, args.max_trials,))
p2 = Process(target=run_bkg_trials, args=(bkg_pdfs, sig_pdfs, names, J_factors, args, outfile, args.num_trials, args.max_trials,))
p1.start()
p2.start()
p1.join()
p2.join()

print("Done")
