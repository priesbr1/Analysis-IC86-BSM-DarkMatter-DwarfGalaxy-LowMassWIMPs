import numpy as np
from math import *
import argparse
import glob
import time
import os, sys
from icecube import icetray, dataclasses, dataio
from I3Tray import *
from icecube.dataclasses import *

parser = argparse.ArgumentParser()
parser.add_argument("-d", "--datafolder", type=str, default="/mnt/research/IceCube/jpandre/Matt/level5p/numu/14674/",
                    dest="datafolder", help="folder containing MC datafiles")
parser.add_argument("-f", "--files", type=str, default="/Level5p_IC86.2013_genie_numu.014674.??????.i3.bz2",
                    dest="files", help="names of MC files")
args = parser.parse_args()

Infile_List = sorted(glob.glob(args.datafolder + args.files))
print("Number of files: %i"%len(Infile_List))

weights = []

def get_weights(frame):
    mctree = frame["I3MCTree"]
    mcparticle = get_most_energetic_primary(mctree)
    mcweightdict = frame["I3MCWeightDict"]
    
    if (abs(mcparticle.type == 14)):
        infile_list_length = len(Infile_List)
        ftype = 0
    else:
        return
    
    if (mcparticle.type > 0):
        typeweight = 0.7
        qtype = 0
    if (mcparticle.type < 0):
        typeweight = 0.3
        qtype = 1
    
    w = frame["NeutrinoWeights"]["OscillatedRate"]/(infile_list_length*typeweight/0.5)
    weights.append(w)

tray = I3Tray()
tray.AddModule("I3Reader", "reader", FilenameList=Infile_List)
tray.AddModule(get_weights, "get_weights")
tray.AddModule("TrashCan", "done")

t1 = time.time()
tray.Execute()
t2 = time.time()
tray.Finish()

print("Execution time: %.2f minutes"%((t2-t1)/60))

weights = np.array(weights)
print("Number of events: %i"%len(weights))
weights = weights[~np.isnan(weights)]
print("Number of events with non-NaN weights: %i"%len(weights))

livetime = np.sum(weights)/np.sum(weights**2)
print("MC effective livetime: %.2f seconds = %.2f days = %.2f years"%(livetime, livetime/86400., livetime/86400./365.))
