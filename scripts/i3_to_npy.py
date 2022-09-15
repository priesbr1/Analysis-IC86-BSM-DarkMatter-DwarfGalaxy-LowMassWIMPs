#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py2-v3.0.1/icetray-start
#METAPROJECT: simulation/V06-01-00-RC4

import os, sys
import glob
import numpy as np
import argparse
from icecube import icetray, dataio, dataclasses, astro
from I3Tray import I3Units

np.seterr(divide = "raise") # Turn DivideByZero warnings into errors for error handling

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--files", type=str, default=None,
                    dest="input_files", help="names of input files with wildcard fillers")
parser.add_argument("-i", "--inpath", type=str, default=None,
                    dest="input_path", help="path to input files")
parser.add_argument("-o", "--outpath", type=str, default=None,
                    dest="output_path", help="destination path for output file")
parser.add_argument("-m", "--MC", type=int, default=0,
                    dest="MC", help="whether or not to process as MC reco")
args = parser.parse_args()

input_files = str(args.input_path + args.input_files)
out_path = args.output_path
MC = bool(args.MC)

def MC_angular_error(zen_reco, azi_reco, zen_mc, azi_mc):
    return np.arccos(np.sin(zen_reco)*np.sin(zen_mc)*np.cos(azi_reco-azi_mc) + np.cos(zen_reco)*np.cos(zen_mc))

def read_files(filename_list, path_to_output_file):
    
    # Set up data array
    if MC:
        data = np.array([], dtype=[('run','i8'), ('event','i8'), ('subevent','i8'), ('ra','f8'), ('dec','f8'), ('azi','f8'), ('zen','f8'), ('time','f8'), ('trueRa','f8'), ('trueDec','f8'), ('trueE','f8'), ('logE','f8'), ('angErr','f8'), ('ow','f8')])
    else:
        data = np.array([], dtype=[('run','i8'), ('event','i8'), ('subevent','i8'), ('ra','f8'), ('dec','f8'), ('azi','f8'), ('zen','f8'), ('time','f8'), ('logE','f8'), ('angErr','f8')])
        grl = np.array([], dtype=[('run','i8'), ('start','f8'), ('stop','f8'), ('livetime','f8'), ('events','i8')])
        run_groups = dict()
    outfile = None    

    for i3_filename in filename_list:
        if not outfile:
            index = i3_filename.find("Level5p")
            year = i3_filename[i3_filename.find("20"):i3_filename.find("20")+4]
            buff = len("Level5p_IC86.")
            if MC:
                MC_set = i3_filename[i3_filename.find("14"):i3_filename.find("14")+5]
                outfile = path_to_output_file + i3_filename[index:index+buff] + MC_set + '.' + year + "_data_MC.npy"
            else:
                outfile = path_to_output_file + i3_filename[index:index+buff] + year + "_data_exp.npy"
                if not os.path.isdir(path_to_output_file + "GRL/"):
                    os.mkdir(path_to_output_file + "GRL/")
                grl_outfile = path_to_output_file + "GRL/" + i3_filename[index:index+buff] + year + "_data_exp.npy"
            print("Saving output to {}".format(outfile))
            if not MC:
                print("Saving GRL output to {}".format(grl_outfile))

        print("Reading file {}".format(i3_filename[i3_filename.find("Level5p"):]))
        i3_file = dataio.I3File(i3_filename)
        frame_counter = 0

        while i3_file.more():

            # Get next P frame if it exists
            try:
                pframe = i3_file.pop_physics()
            except:
                continue

            event = pframe["I3EventHeader"]
            run_id = event.run_id
            event_id = event.event_id
            subevent_id = event.sub_event_id
            mjd = event.start_time.mod_julian_day_double

            if not MC:
                if run_id not in run_groups.keys():
                    run_groups[run_id] = []
                run_groups[run_id].append(event)
            
            if MC:
                try:
                    nu = pframe["IC86_Dunkman_L6_PegLeg_MultiNest8D_NumuCC"]
                except KeyError:
                    print('Skipping frame with no "IC86_Dunkman_L6_PegLeg_MultiNest8D_NumuCC"')
                    continue
                try:
                    logE = np.log10(nu.energy / I3Units.GeV)
                except FloatingPointError:
                    print("Skipping event with energy = 0 GeV")
                    continue
                azi = nu.dir.azimuth / I3Units.rad
                zen = nu.dir.zenith / I3Units.rad
                ra, dec = astro.dir_to_equa(zen, azi, mjd)

                MC_nu = dataclasses.get_most_energetic_primary(pframe["I3MCTree"])
                trueE = MC_nu.energy / I3Units.GeV
                trueAzi = MC_nu.dir.azimuth / I3Units.rad
                trueZen = MC_nu.dir.zenith / I3Units.rad
                trueRA, trueDec = astro.dir_to_equa(trueZen, trueAzi, mjd)
                angErr = MC_angular_error(azi, zen, trueAzi, trueZen)
                
                try:
                    MC_weightdict = pframe["I3MCWeightDict"]
                except KeyError:
                    print("Skipping event without OneWeight")
                    continue
                    
                if MC_nu.pdg_encoding < 0:
                    nu_nubar = 0.3 # Anti-neutrino
                else:
                    nu_nubar = 0.7 # Neutrino
                ow = MC_weightdict["OneWeight"]/MC_weightdict["NEvents"]/len(filename_list)/nu_nubar * 2

                data = np.append(data, np.array([(run_id, event_id, subevent_id, ra, dec, azi, zen, mjd, trueRA, trueDec, trueE, logE, angErr, ow)], dtype=data.dtype))

            else:
                try:
                    nu = pframe["IC86_Dunkman_L6_PegLeg_MultiNest8D_"]
                except KeyError:
                    print('Skipping frame with no "IC86_Dunkman_L6_PegLeg_MultiNest8D_"')
                    continue
                try:
                    logE = np.log10(nu.energy / I3Units.GeV)
                except FloatingPointError:
                    print("Skipping event with energy = 0 GeV")
                    continue
                azi = nu.dir.azimuth / I3Units.rad
                zen = nu.dir.zenith / I3Units.rad
                ra, dec = astro.dir_to_equa(zen, azi, mjd)
                angErr = 0.89318577 * np.exp(-0.03439216 * nu.energy / I3Units.GeV) + 0.08273342 # Curvefit from https://wiki.icecube.wisc.edu/index.php/File:Ang_res_tot.png

                data = np.append(data, np.array([(run_id, event_id, subevent_id, ra, dec, azi, zen, mjd, logE, angErr)], dtype=data.dtype))

            frame_counter += 1

        if not MC:
            runs = sorted(list(run_groups.keys()))
            for run_id in runs:
                run_events = run_groups[run_id]
                start_time = np.min([event.start_time.mod_julian_day_double for event in run_events])
                stop_time = np.max([event.end_time.mod_julian_day_double for event in run_events])
                grl = np.append(grl, np.array([(run_id, start_time, stop_time, stop_time-start_time, len(run_events))], dtype=grl.dtype))

        print("Finished reading file {} ({} P frames)".format(i3_filename[i3_filename.find("Level5p"):], frame_counter))

    print("Saving output as {} ({} files)".format(outfile[outfile.find("Level5p"):], len(filename_list)))
    np.save(outfile, data)
    if not MC:
        print("Saving GRL output as {}".format(grl_outfile[grl_outfile.find("GRL"):]))
        np.save(grl_outfile, grl)

if '*' in input_files or '?' in input_files:
    input_files = sorted(glob.glob(input_files))

if isinstance(input_files, str):
    input_files = [input_files]
if not isinstance(input_files, list):
    raise TypeError("Input files are not a string or list")

read_files(input_files, out_path)
