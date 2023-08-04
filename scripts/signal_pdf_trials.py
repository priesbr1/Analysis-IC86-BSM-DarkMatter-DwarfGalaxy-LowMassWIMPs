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
parser.add_argument("-j", "--J_type", type=str, choices=["0.1","0.2","0.5","10","min","max"], default="max",
                    dest="J_type", help="J-factor values to use based on opening half-angle")
parser.add_argument("-c", "--channel", type=str, default="b",
                    dest="channel", help="annhiliation channel to run background trials for")
parser.add_argument("-m", "--mass", type=int, default=10,
                    dest="mass", help="WIMP mass to rub background trials for")
parser.add_argument("-e", "--num_bkg_events", type=int, default=1000,
                    dest="num_bkg_events", help="number of background events to sample")
parser.add_argument("-n", "--num_trials", type=int, default=1,
                    dest="num_trials", help="number of signal trials to run per core")
parser.add_argument("-d", "--seed", type=int, default=None,
                    dest="seed", help="RNG seed for reproducibility")
parser.add_argument("-o", "--outfolder", type=str, default="/mnt/home/priesbr1/DM_Search/data/trials_results/trials_sig_360deg_Jmax_distributed/",
                    dest="outfolder", help="folder to save signal trials to")
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

outfile = "sig_trials_" + args.channel + "_" + str(args.mass) + ".txt"

def run_one_trial(ni_source, data, bkg_pdfs, sig_pdfs, source, channel, mass, J_factors):
    
    # Recover signal
    N = len(data)
    llh = []
    ll = 0
    
    bkg_bins = np.searchsorted(bkg_pdfs[source][channel][mass]["Bins"], data, "right")-1 # Find left edges of bins
    bkg_probs = np.array(bkg_pdfs[source][channel][mass]["Counts"][bkg_bins], dtype=float) # Get probabilities for bins
    bkg_probs /= np.sum(bkg_pdfs[source][channel][mass]["Counts"]) # Normalize probabilities
    sig_bins = np.searchsorted(sig_pdfs[channel][mass]["Bins"], data, "right")-1
    sig_probs = np.array(sig_pdfs[channel][mass]["Counts"][sig_bins], dtype=float)
    sig_probs /= np.sum(sig_pdfs[channel][mass]["Counts"])
    
    if (ni_source != 0):
        for ns in range(int(ni_source*1.05)):
            ns = float(ns)
            combs = (ns/N)*sig_probs + (1-(ns/N))*bkg_probs
            llh.append(-1*np.sum(np.log(combs)))
        
        bf = np.argmin(llh)
        ts = 2*(llh[0] - llh[bf])
    
    else:
        for ns in range(int(len(data)/100)):
            ns = float(ns)
            combs = (ns/N)*sig_probs + (1-(ns/N))*bkg_probs
            llh.append(-1*np.sum(np.log(combs)))
    
        bf = np.argmin(llh)
        bf = int(np.round(bf*J_factors[names.index(source)]/np.sum(J_factors)))
        ts = 2*(llh[0] - llh[bf])
        ts *= J_factors[names.index(source)]/np.sum(J_factors)
    
    return bf, ts

#in_events = np.concatenate((np.arange(1,10,1,dtype=int), np.linspace(10,50,20,dtype=int)))
#in_events = np.concatenate((np.arange(1,10,1,dtype=int), np.linspace(10,30,20,dtype=int,endpoint=False), np.arange(30,300+30/2,30,dtype=int)))
in_events = np.concatenate((np.arange(10,100,10,dtype=int),np.arange(100,1000+30/2,30,dtype=int)))
#in_events = np.arange(100,1000+30/2,30,dtype=int)
#in_events = np.concatenate((np.arange(10,100,10,dtype=int),np.arange(100,1000,30,dtype=int),np.arange(1000,10000,1000,dtype=int),np.arange(10000,100000+10000/2,10000,dtype=int)))
#in_events = np.concatenate((np.arange(100,1000,30,dtype=int),np.arange(1000,10000,1000,dtype=int),np.arange(10000,100000+10000/2,10000,dtype=int)))

t1 = time.time()

data_bkg = []
for source in names:
    counts, bins = bkg_pdfs[source][args.channel][args.mass]["Counts"], bkg_pdfs[source][args.channel][args.mass]["Bins"]
    data_bkg.append(generate_rand_from_pdf(counts,bins[:-1],args.num_bkg_events))

for ni in in_events:
    rec = []
    ts_rec = []
    
    for x in range(args.num_trials):
        t1_trial = time.time()
        bfs = []
        tss = []
        
        for j, source, in enumerate(names):
            ni_source = int(round(ni*J_factors[j]/np.sum(J_factors)))
            counts, bins = sig_pdfs[args.channel][args.mass]["Counts"], sig_pdfs[args.channel][args.mass]["Bins"]
            data_sig = generate_rand_from_pdf(counts,bins[:-1],ni_source)
            data = np.concatenate((data_bkg[j], data_sig))

            bf_source, ts_source = run_one_trial(ni_source, data, bkg_pdfs, sig_pdfs, source, args.channel, args.mass, J_factors)
            
            bfs.append(bf_source)
            tss.append(ts_source)
        
        bf_total = int(np.round(np.sum(bfs)))
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
