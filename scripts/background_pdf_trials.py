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
parser.add_argument("-b", "--bkg_pdfs", type=str, default="/mnt/home/priesbr1/DM_Search/data/pdfs/DRAGON_background_pdfs_5deg.npy",
                    dest="bkg_pdfs", help="file containing background PDFs for each mass, channel, and source")
parser.add_argument("-p", "--sig_pdfs", type=str, default="/mnt/home/priesbr1/DM_Search/data/pdfs/DRAGON_signal_pdfs.npy",
                    dest="sig_pdfs", help="file containing signal PDFs for each mass and channel")
parser.add_argument("-s", "--sources", type=str, default="/mnt/home/priesbr1/DM_Search/data/analysis_sources_ra_dec_jfactors.txt",
                    dest="sources", help="file containing RA, dec, and J-factor information for each source")
parser.add_argument("-c", "--channel", type=str, default="b",
                    dest="channel", help="annhiliation channel to run background trials for")
parser.add_argument("-m", "--mass", type=int, default=10,
                    dest="mass", help="WIMP mass to rub background trials for")
parser.add_argument("-n", "--num_trials", type=int, default=10,
                    dest="num_trials", help="number of background trials to run per core")
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

# Source encoding - 0:name; 1:RA; 2:Dec; 3,4,5:r_h; 6,7,8:J0.1; 9,10,11:J0.2; 12,13,14:J0.5; 15,16,17:J10
sources = np.genfromtxt(args.sources, delimiter=",", dtype=("U20",float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float))
names = []
ra  = np.array([])
dec = np.array([])
J_factors = np.array([])
for i in range(len(sources)):
    for idx in [15,12,9,6]:
        if not np.isnan(sources[i][idx]):
            names.append(sources[i][0])
            ra = np.append(ra, sources[i][1]*np.pi/180.)
            dec = np.append(dec, sources[i][2]*np.pi/180.)
            J_factors = np.append(J_factors, sources[i][idx])
            break
J_factors = 10**J_factors
#print(J_factors/np.sum(J_factors))

outfile = "bkg_trials_" + args.channel + "_" + str(args.mass) + ".txt"

def run_bkg_trials(bkg_pdfs, sig_pdfs, names, J_factors, args, outfile, num_bkg_trials):
    t1 = time.time()
    
    for t in range(num_bkg_trials):
        t1_trial = time.time()
        data = []
        for source in names:
            counts, bins = bkg_pdfs[source][args.channel][args.mass]["Counts"], bkg_pdfs[source][args.channel][args.mass]["Bins"]
            data.append(generate_rand_from_pdf(counts,bins[:-1],1000))
        
        # Fit
        N = int(np.sum([len(samp) for samp in data])/len(names))
        llhs = []
        
        for j, source in enumerate(names):
            source_llh = []
            for i in range(N):
                ll = 0
                i = float(i)
                for n in data[j]:
                    bkg_bin = np.max(np.where((bkg_pdfs[source][args.channel][args.mass]["Bins"] <= n) == True)) # Find left edge of bin
                    bkg_prob = bkg_pdfs[source][args.channel][args.mass]["Counts"][bkg_bin] # Get probability of being in that bin
                    bkg_prob /= np.sum(bkg_pdfs[source][args.channel][args.mass]["Counts"]) # Normalize
                    sig_bin = np.max(np.where((sig_pdfs[args.channel][args.mass]["Bins"] <= n) == True))
                    sig_prob = sig_pdfs[args.channel][args.mass]["Counts"][sig_bin]
                    sig_prob /= np.sum(sig_pdfs[args.channel][args.mass]["Counts"])
                    comb = (i/N)*sig_prob + (1-(i/N))*bkg_prob
                    ll += np.log(comb)
                
                source_llh.append(-ll)
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
                    f = open(args.outfolder+outfile,"a")
                    for j, source in enumerate(names):
                        f.write("%s \t %s \t "%(bfs[j], tss[j]))
                    f.write("%s \t %s \n"%(bf,ts))
                    f.close()
                    written = True
                except:
                    print("Waiting for file access...")
                    time.sleep(1)
    
    t2 = time.time()
    print("Time: %.2f sec"%(t2-t1))

p1 = Process(target=run_bkg_trials, args=(bkg_pdfs, sig_pdfs, names, J_factors, args, outfile, args.num_trials,))
p2 = Process(target=run_bkg_trials, args=(bkg_pdfs, sig_pdfs, names, J_factors, args, outfile, args.num_trials,))
p1.start()
p2.start()
p1.join()
p2.join()

print("Done")
