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
parser.add_argument("-d", "--seed", type=int, default=None,
                    dest="seed", help="RNG seed for reproducibility")
parser.add_argument("-o", "--outfolder", type=str, default="/mnt/home/priesbr1/DM_Search/data/trials_results/trials_sig_360deg_Jmax_modified/",
                    dest="outfolder", help="folder to save signal trials to")
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

#in_events = np.concatenate((np.arange(1,10,1,dtype=int), np.linspace(10,50,20,dtype=int)))
in_events = np.arange(30,300+30/2,30,dtype=int)
num_sig_trials = 1

t1 = time.time()

for ni in in_events:
    rec = []
    ts_rec = []
    ni_sources = np.random.choice(len(sources), ni, replace=True, p=J_factors/np.sum(J_factors)) # Choose sources based on source weight
    data = []
    
    for x in range(num_sig_trials):
        t1_trial = time.time()
        counts, bins = sig_pdfs[args.channel][args.mass]["Counts"], sig_pdfs[args.channel][args.mass]["Bins"]
        for j, source in enumerate(names):
            data.append(generate_rand_from_pdf(counts,bins[:-1],ni_sources[j]))
        
        # Recover signal 
        N = len(data)
        llhs = []
        
        for j, source in enumerate(names):
            source_llh = []
            for i in range(int(N*1.05)):
                ll = 0
                i = float(i)
                for n in data[j]:
                    bkg_bin = np.max(np.where((bkg_pdfs[source][args.channel][args.mass]["Bins"] <= n) == True)) # Find left edge of bin
                    bkg_prob = bkg_pdfs[source][args.channel][args.mass]["Counts"][bkg_bin] # Get probability of being in that bin
                    bkg_prob /= np.sum(bkg_pdfs[source][args.channel][args.mass]["Counts"])
                    sig_bin = np.max(np.where((sig_pdfs[args.channel][args.mass]["Bins"] <= n) == True))
                    sig_prob = sig_pdfs[args.channel][args.mass]["Counts"][sig_bin]
                    sig_prob /= np.sum(sig_pdfs[args.channel][args.mass]["Counts"])
                    comb = (i/N)*sig_prob + (1-(i/N))*bkg_prob
                    try:
                        ll += np.log(comb)
                    except:
                        pass # Datapoint is guaranteed to be background --> comb = 1 --> log(comb) = 0
            
                source_llh.append(-ll)
            llhs.append(source_llh)
          
        bfs = np.argmin(llhs, axis=1)
        tss = []
        for L,bf in zip(llhs,bfs):
            tss.append(2*(L[0] - L[bf]))
        
        bf = int(np.round(np.sum(bfs * J_factors/np.sum(J_factors))))
        ts = np.sum(np.array(tss))
        rec.append(bf)
        ts_rec.append(ts)
        
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
                    f.write("%s \t %s \n"%(bf,ts))
                    f.close()
                    written = True
                except:
                    print("Waiting for file access...")
                    time.sleep(1)
    
t2 = time.time()
print("Time: %.2f sec"%(t2-t1))

print("Done")
