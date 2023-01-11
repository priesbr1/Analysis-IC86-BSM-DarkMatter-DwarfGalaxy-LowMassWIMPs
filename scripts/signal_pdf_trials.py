import numpy as np
import argparse
import os,time
np.seterr(all="raise")

def generate_rand_from_pdf(pdf, x_grid, n):
    cdf = np.cumsum(pdf)
    cdf = cdf / cdf[-1]
    values = np.random.rand(n)
    value_bins = np.searchsorted(cdf, values)
    random_from_cdf = x_grid[value_bins]
    return random_from_cdf

parser = argparse.ArgumentParser()
parser.add_argument("-b", "--bkg_pdfs", type=str, default="/mnt/home/priesbr1/DM_Search/data/pdfs/DRAGON_background_pdfs_360deg.npy",
                    dest="bkg_pdfs", help="file containing background PDFs for each mass, channel, and source")
parser.add_argument("-p", "--sig_pdfs", type=str, default="/mnt/home/priesbr1/DM_Search/data/pdfs/DRAGON_signal_pdfs.npy",
                    dest="sig_pdfs", help="file containing signal PDFs for each mass and channel")
parser.add_argument("-s", "--sources", type=str, default="/mnt/home/priesbr1/DM_Search/data/analysis_sources_ra_dec_jfactors.txt",
                    dest="sources", help="file containing RA, dec, and J-factor information for each source")
parser.add_argument("-c", "--channel", type=str, default="b",
                    dest="channel", help="annhiliation channel to run background trials for")
parser.add_argument("-m", "--mass", type=int, default=10,
                    dest="mass", help="WIMP mass to rub background trials for")
parser.add_argument("-n", "--num_trials", type=int, default=1,
                    dest="num_trials", help="number of signal trials to run per core")
parser.add_argument("-d", "--seed", type=int, default=None,
                    dest="seed", help="RNG seed for reproducibility")
parser.add_argument("-o", "--outfolder", type=str, default="/mnt/home/priesbr1/DM_Search/data/trials_results/trials_sig_360deg_Jmax_distributed/",
                    dest="outfolder", help="folder to save signal trials to")
parser.add_argument("-f", "--bkg_folder", type=str, default="/mnt/home/priesbr1/DM_Search/data/trials_results/trials_bkg_360deg_Jmax/",
                    dest="bkg_folder", help="folder containing results from corresponding background trials")
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

outfile = "sig_trials_" + args.channel + "_" + str(args.mass) + ".txt"

def run_one_trial(ni, bkg_pdfs, sig_pdfs, source, channel, mass):
    counts, bins = sig_pdfs[channel][mass]["Counts"], sig_pdfs[channel][mass]["Bins"]
    data = generate_rand_from_pdf(counts,bins[:-1],ni)
    
    # Recover signal
    N = len(data)
    llh = []
    ll = 0
    
    for i in range(int(N*1.05)):
        ll = 0
        i = float(i)
        for n in data:
            bkg_bin = np.max(np.where((bkg_pdfs[source][channel][mass]["Bins"] <= n) == True)) # Find left edge of bin
            bkg_prob = bkg_pdfs[source][channel][mass]["Counts"][bkg_bin] # Get probability of being in that bin
            bkg_prob /= np.sum(bkg_pdfs[source][channel][mass]["Counts"])
            sig_bin = np.max(np.where((sig_pdfs[channel][mass]["Bins"] <= n) == True))
            sig_prob = sig_pdfs[channel][mass]["Counts"][sig_bin]
            sig_prob /= np.sum(sig_pdfs[channel][mass]["Counts"])
            comb = (i/N)*sig_prob + (1-(i/N))*bkg_prob
            try:
                ll += np.log(comb)
            except:
                ll = 0 # Datapoint is guaranteed to be background --> comb = 1 --> log(comb) = 0
        
        llh.append(-ll)
    
    bf = np.argmin(llh)
    ts = 2*(llh[0] - llh[bf])
    
    return bf, ts

# Load in background trial information
bkg_ns = []
bkg_ts = []
for i, source in enumerate(names):
    source_ns = []
    source_ts = []
    bkg_trials = np.genfromtxt(args.bkg_folder+outfile.replace("sig","bkg"), delimiter="\t")
    for t in range(len(bkg_trials)):
        source_ns.append(bkg_trials[t][2*i])
        source_ts.append(bkg_trials[t][2*i+1])
    bkg_ns.append(np.median(source_ns))
    bkg_ts.append(np.median(source_ts))

#in_events = np.concatenate((np.arange(1,10,1,dtype=int), np.linspace(10,50,20,dtype=int)))
in_events = np.concatenate((np.arange(1,10,1,dtype=int), np.linspace(10,30,20,dtype=int,endpoint=False), np.arange(30,300+30/2,30,dtype=int)))

t1 = time.time()

for ni in in_events:
    rec = []
    ts_rec = []
    
    for x in range(args.num_trials):
        t1_trial = time.time()
        bfs = []
        tss = []
        
        for source,jf in zip(names,J_factors):
            ni_source = int(round(ni*jf/np.sum(J_factors)))
            
            if (ni_source == 0): # Background trial
                bf_source = bkg_ns[names.index(source)]
                ts_source = bkg_ts[names.index(source)]
            else: # Signal trial
                bf_source, ts_source = run_one_trial(ni_source, bkg_pdfs, sig_pdfs, source, args.channel, args.mass)
            
            bfs.append(bf_source)
            tss.append(ts_source)
        
        bf_total = int(round(np.sum(bfs)))
        ts_total = np.sum(np.array(tss))
        rec.append(bf_total)
        ts_rec.append(ts_total)
        
        t2_trial = time.time()
        print("Trial time for n_inj=%i: %.2f sec"%(ni, t2_trial-t1_trial))
        
        written = False
        while (written == False):
            if not os.path.isfile(args.outfolder+outfile):
                f = open(args.outfolder+outfile,"w")
                f.write("# ni \t ")
                for j, source in enumerate(names):
                    f.write("%s_ns \t %s_TS \t "%(source, source))
                f.write("ns \t TS \n")
                f.close()
            else:
                try:
                    os.rename(args.outfolder+outfile,args.outfolder+outfile)
                    f = open(args.outfolder+outfile,"a")
                    f.write("%s \t "%ni)
                    for j, source in enumerate(names):
                        f.write("%s \t %s \t "%(bfs[j], tss[j]))
                    f.write("%s \t %s \n"%(bf_total,ts_total))
                    f.close()
                    written = True
                except:
                    print("Waiting for file access...")
                    time.sleep(1)
    
t2 = time.time()
print("Time: %.2f sec"%(t2-t1))

print("Done")
