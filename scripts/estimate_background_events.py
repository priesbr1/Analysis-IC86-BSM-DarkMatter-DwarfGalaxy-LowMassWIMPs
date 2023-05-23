import numpy as np
import glob
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-d", "--data", type=str, default=None,
                    dest="data", help="file(s) containing data to generate background pdf from")
parser.add_argument("-s", "--sources", type=str, default="/mnt/home/priesbr1/DM_Search/data/analysis_sources_ra_dec_jfactors.txt",
                    dest="sources", help="file containing RA, dec, and J-factor information for each source")
parser.add_argument("-o", "--sigma", type=int, default=2,
                    dest="sigma", help="sigma value for angular error tolerance per event")
parser.add_argument("-e", "--seed", type=int, default=None,
                    dest="seed", help="RNG seed for reproducibility")
args = parser.parse_args()

if (args.seed != None):
    np.random.seed(args.seed)

def calculate_angle(dec_source, ra_source, dec_nu, ra_nu):
    return np.arccos(np.cos(dec_source)*np.cos(dec_nu)*np.cos(ra_source - ra_nu) + np.sin(dec_source)*np.sin(dec_nu))

# Source encoding - 0:name; 1:RA; 2:Dec; 3,4,5:r_h; 6,7,8:J0.1; 9,10,11:J0.2; 12,13,14:J0.5; 15,16,17:J10
sources = np.genfromtxt(args.sources, delimiter=",", dtype=("U20",float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float))
names = []
ra  = np.array([])
dec = np.array([])
for i in range(len(sources)):
    names.append(sources[i][0])
    ra = np.append(ra, sources[i][1]*np.pi/180.)
    dec = np.append(dec, sources[i][2]*np.pi/180.)
max_name_length = np.max([len(name) for name in names])

background_counters = np.zeros(len(names), dtype=int)
event_counter = 0

if "*" in args.data or "?" in args.data:
    filenames = sorted(glob.glob(args.data))
else:
    filenames = [args.data]
print("Using data from:")
for filename in filenames:
    print(" - "+str(filename))
    data = np.load(filename)
    
    for i in range(len(data["run"])):
        nu_dec = data["dec"][i]
        nu_ra = np.random.uniform(0, 2*np.pi)
        nu_err = data["angErr"][i]
        
        for j in range(len(names)):
            if (calculate_angle(dec[j], ra[j], nu_dec, nu_ra) <= args.sigma*nu_err):
                background_counters[j] += 1
        event_counter += 1

print("-"*10)

print("Estimated background events:")
for j in range(len(names)):
    print("  - {0:<{1}}: {2}/{3}   ({4:.4f})".format(names[j], max_name_length, background_counters[j], event_counter, float(background_counters[j]/event_counter)))
