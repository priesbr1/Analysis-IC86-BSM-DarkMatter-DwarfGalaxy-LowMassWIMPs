import numpy as np
import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
import healpy as hp
import argparse
import os, sys

def cat2hpx(ra, dec, nside):
    """
    Modified from https://stackoverflow.com/questions/50483279/make-a-2d-histogram-with-healpix-pixellization-using-healpy
    
    Convert datapoints to a HEALPix map as a 2D histogram.

    Inputs
    ----------
    ra, dec : (ndarray, ndarray)
        Coordinates of the sources in degrees. Assumes ICRS RA/declination coordinates.

    nside : int
        HEALPix nside of the target map

    Returns
    ------
    hpx_map : ndarray
        HEALPix map of the bin counts in ICRS coordinates
    """

    npix = hp.nside2npix(nside)

    # convert to theta, phi
    theta = np.radians(90. - dec)
    phi = np.radians(ra)

    # convert to HEALPix indices
    indices = hp.ang2pix(nside, theta, phi)

    idx, counts = np.unique(indices, return_counts=True)

    # fill the fullsky map
    hpx_map = np.zeros(npix, dtype=int)
    hpx_map[idx] = counts

    return hpx_map

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input_file", type=str, default=None,
                    dest="input_file", help="full path and name of input .npy file")
parser.add_argument("-o", "--output_folder", type=str, default=None,
                    dest="output_folder", help="parent directory for output plots")
parser.add_argument("-d", "--degrees", type=int, default=1,
                    dest="degrees", help="whether or not to convert raidans to degrees")
parser.add_argument("-r", "--ra_hms", type=int, default=0,
                    dest="ra_hms", help="whether or not to convert RA to hms format")
parser.add_argument("-p", "--ps", type=int, default=0,
                    dest="ps", help="whether or not sample is from point_source_tracks")
parser.add_argument("-u", "--use_sources", type=int, default=0,
                    dest="use_sources", help="whether or not to use sources information in plotting")
parser.add_argument("-s", "--sources_file", type=str, default="None",
                    dest="sources_file", help="full path and name of file with known sources")
parser.add_argument("-e", "--seed", type=int, default=None,
                    dest="seed", help="RNG seed for reproducibility")
args = parser.parse_args()

if (args.seed != None):
    np.random.seed(args.seed)

filename = args.input_file
save_folder = args.output_folder
deg = bool(args.degrees)
hms = bool(args.ra_hms)
ps = bool(args.ps)
use_sources = bool(args.use_sources)
sources_file = args.sources_file

if ps:
    if save_folder.find("ps_tracks") != -1:
        print("Overwriting output_folder to save to appropriate point source directory")
        save_folder.replace("i3", "ps", 1)

data = np.load(filename)
MC = bool(filename.find("MC") >= 0)
year = filename[filename.find("20"):filename.find("20")+4]
if ps:
    version = filename[filename.find("version"):filename.rfind('/')]
    version = version.replace("ersion", '').replace('-', '')
    save_folder += version + '_' + year + "_data/"
elif MC:
    MC_set = filename[filename.find("14"):filename.find("14")+5]
    save_folder += year + '_' + MC_set + "_data/"
else:
    save_folder += year + "_data/"
if os.path.isdir(save_folder) != True:
    os.mkdir(save_folder)

if ps and MC:
    print("Plotting point source MC distributions from %s version %s"%(year, version))
elif ps:
    print("Plotting point source data distributions from %s version %s"%(year, version))
elif MC:
    print("Plotting MC distributions from %s set %s"%(year, MC_set))
else:
    print("Plotting data distributions from %s"%year)
print("Saving plots to %s"%save_folder)

if use_sources and sources_file != "None":
    # Source encoding - 0:name; 1:RA; 2:Dec; 3,4,5:r_h; 6,7,8:J0.1; 9,10,11:J0.2; 12,13,14:J0.5; 15,16,17:J10
    sources = np.genfromtxt(args.sources, delimiter=",", dtype=("U20",float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float))
    source_names = []
    source_ra  = np.array([])
    source_dec = np.array([])
    for i in range(len(sources)):
        source_names.append(sources[i][0])
        source_ra = np.append(source_ra, np.radians(sources[i][1]))
        source_dec = np.append(source_dec, np.radians(sources[i][2]))
    print("Using sources from %s"%sources_file)

    dwarf_galaxies_file = open(save_folder + "source_names.txt", 'w')
    dwarf_galaxies_file.write("Known dwarf galaxies used in latest plots:\n\n")
    for name in source_names:
        dwarf_galaxies_file.write("%s\n"%name)
    dwarf_galaxies_file.close()

if deg:
    angle_units = "deg"
    data["ra"] = np.degrees(data["ra"])
    ra_range = (0,360)
    data["dec"] = np.degrees(data["dec"])
    dec_range = (-90,90)
    data["azi"] = np.degrees(data["azi"])
    azi_range = (0,360)
    data["zen"] = np.degrees(data["zen"])
    zen_range = (0,180)
    data["angErr"] = np.degrees(data["angErr"])
    if MC:
        data["trueRa"] = np.degrees(data["trueRa"])
        data["trueDec"] = np.degrees(data["trueDec"])
    if use_sources and sources_file != "None":
        source_ra = np.degrees(source_ra)
        source_dec = np.degrees(source_dec)
    if hms:
        ra_units = "hrs"
        data["ra"] *= 24.0/360.0
        ra_range = (0,24)
        if MC:
            data["trueRa"] *= 24.0/360.0
        if sources_file != "None":
            source_ra *= 24.0/360.0
    else:
        ra_units = angle_units
else:
    angle_units = "rad"
    ra_range = (0,2*np.pi)
    dec_range = (-np.pi/2,np.pi/2)
    azi_range = (0,2*np.pi)
    zen_range = (0,np.pi)
    if hms:
        ra_units = "hrs"
        data["ra"] *= 12.0*np.pi # 180/np.pi*24/360
        ra_range = (0,24)
        if MC:
            data["trueRa"] *= 12.0*np.pi
        if use_sources and sources_file != "None":
            source_ra *= 12.0*np.pi
    else:
        ra_units = angle_units
    
plt.figure()
plt.hist(data["run"], bins=100)
plt.xlabel("Run Number")
plt.ylabel("Counts")
plt.title("Run Number Distribution, n=%i"%len(data["run"]))
plt.tight_layout()
plt.savefig(save_folder + "run.png")
plt.close()

plt.figure()
plt.hist(data["event"], bins=100)
plt.xlabel("Event Number")
plt.ylabel("Counts")
plt.title("Event Number Distribution, n=%i"%len(data["event"]))
plt.tight_layout()
plt.savefig(save_folder + "event.png")
plt.close()

plt.figure()
plt.hist(data["subevent"], bins=100)
plt.xlabel("Subevent Number")
plt.ylabel("Counts")
plt.title("Subevent Number Distribution, n=%i"%len(data["subevent"]))
plt.tight_layout()
plt.savefig(save_folder + "subevent.png")
plt.close()

plt.figure()
if MC:
    plt.hist(data["ra"], bins=100, range=ra_range, alpha=0.5, label="Reco")
    plt.hist(data["trueRa"], bins=100, range=ra_range, alpha=0.5, label="MC Truth")
    plt.legend(loc="best")
else:
    plt.hist(data["ra"], bins=100, range=ra_range, label="Reco")
plt.xlabel("Right Ascension [" + ra_units + ']')
plt.ylabel("Counts")
plt.title("Right Ascension Distribution, n=%i"%len(data["ra"]))
plt.tight_layout()
plt.savefig(save_folder + "ra_" + ra_units + ".png")
plt.close()

plt.figure()
if MC:
    plt.hist(data["dec"], bins=100, range=dec_range, alpha=0.5, label="Reco")
    plt.hist(data["trueDec"], bins=100, range=dec_range, alpha=0.5, label="MC Truth")
    plt.legend(loc="best")
else:
    plt.hist(data["dec"], bins=100, range=dec_range, label="Reco")
plt.xlabel("Declination [" + angle_units + ']')
plt.ylabel("Counts")
plt.title("Declination Distribution, n=%i"%len(data["dec"]))
plt.tight_layout()
plt.savefig(save_folder + "dec_" + angle_units + ".png")
plt.close()

plt.figure()
plt.hist(data["azi"], bins=100, range=azi_range)
plt.xlabel("Azimuth [" + angle_units + ']')
plt.ylabel("Counts")
plt.title("Azimuth Distribution, n=%i"%len(data["azi"]))
plt.tight_layout()
plt.savefig(save_folder + "azimuth_" + angle_units + ".png")
plt.close()

plt.figure()
plt.hist(data["zen"], bins=100, range=zen_range)
plt.xlabel("Zenith [" + angle_units + ']')
plt.ylabel("Counts")
plt.title("Zenith Distribution, n=%i"%len(data["zen"]))
plt.tight_layout()
plt.savefig(save_folder + "zenith_" + angle_units + ".png")
plt.close()

plt.figure()
plt.hist(data["time"], bins=100)
plt.xlabel("Event Start Time [MJD]")
plt.ylabel("Counts")
plt.title("MJD Event Time Distribution, n=%i"%len(data["time"]))
plt.tight_layout()
plt.savefig(save_folder + "mjd.png")
plt.close()

plt.figure()
if MC:
    plt.hist(data["logE"], bins=100, alpha=0.5, label="Reco")
    plt.hist(np.log10(data["trueE"]), bins=100, alpha=0.5, label="MC Truth")
    plt.legend(loc="best")
else:
    plt.hist(data["logE"], bins=100, label="Reco")
plt.xlabel("log(Energy/GeV)")
plt.ylabel("Counts")
plt.title("log(Energy) Distribution, n=%i"%len(data["logE"]))
plt.tight_layout()
plt.savefig(save_folder + "logEnergy.png")
plt.close()

plt.figure()
plt.hist(data["angErr"], bins=100, range=(0,np.max(data["angErr"])))
plt.xlabel("Angular Error [" + angle_units + ']')
plt.ylabel("Counts")
plt.title("Angular Error Distribution, n=%i"%len(data["angErr"]))
plt.tight_layout()
plt.savefig(save_folder + "angError_" + angle_units + ".png")
plt.close()

if deg:
    sindec = np.sin(data["dec"]*np.pi/180)
else:
    sindec = np.sin(data["dec"])

plt.figure()
plt.scatter(sindec, data["angErr"], marker=".", s=1/10)
plt.xlabel("Sin(Declination)")
plt.ylabel("Angular Error [" + angle_units + ']')
plt.title("Angular Error vs. Sin(Declination)")
plt.tight_layout()
plt.savefig(save_folder + "angError_" + angle_units + "_sindec.png")
plt.close()

plt.figure()
scrambled = np.random.uniform(ra_range[0], ra_range[1], len(data["ra"]))
plt.hist2d(scrambled, data["dec"], bins=100, norm=matplotlib.colors.LogNorm())
plt.xlabel("Right Ascension [" + ra_units + ']')
plt.ylabel("Declination [" + angle_units + ']')
bar = plt.colorbar()
bar.set_label("Counts")
if use_sources and sources_file != "None":
    plt.scatter(source_ra[0], source_dec[0], color="red", marker="+", label="Dwarf Galaxies")
    for i in range(1, len(source_ra)):
        plt.scatter(source_ra[i], source_dec[i], color="red", marker="+", s=75)
    plt.legend(loc="best")
    plt.title("Scrambled RA and Declination with Dwarf Galaxies")
    plt.tight_layout()
    plt.savefig(save_folder + "ra_" + ra_units + "_dec_" + angle_units + "_sources.png")
else:
    plt.title("Scrambled RA and Declination")
    plt.tight_layout()
    plt.savefig(save_folder + "ra_" + ra_units + "_dec_" + angle_units + ".png")
plt.close()

fig = plt.figure(num=1, figsize=(11,6))
if angle_units == "deg":
    hpx_map = cat2hpx(scrambled, data["dec"], nside=32)
else:
    hpx_map = cat2hpx(np.degrees(scrambled), np.degrees(data["dec"]), nside=32)
if use_sources and sources_file != "None":
    hp.mollview(hpx_map, fig=1, coord='C', rot=180, title="Scrambled Data with Sources", norm="hist")
else:
    hp.mollview(hpx_map, fig=1, coord='C', rot=180, title="Scrambled Data", norm="hist")
plt.xlabel("Scrambled Right Ascension [hrs]")
plt.ylabel("Declination [degrees]")
ax1 = fig.get_axes()[0]
ax1.text(2.05,0., r"$0^\circ$", ha="left", va="center")
ax1.text(-2.05,0., r"$360^\circ$", ha="right", va="center")
ax1.text(1.9,0.45, r"$30^\circ$", ha="left", va="center")
ax1.text(1.4,0.8, r"$60^\circ$", ha="left", va="center")
ax1.text(1.9,-0.45, r"$-30^\circ$", ha="left", va="center")
ax1.text(1.4,-0.8, r"$-60^\circ$", ha="left", va="center")
hp.graticule()
if use_sources and sources_file != "None":
    if angle_units == "rad":
        hp.projscatter(np.pi/2-source_dec, source_ra, color="red", marker='+', s=75, label="Dwarf Galaxies")
    else:
        hp.projscatter(np.pi/2-np.radians(source_dec), np.radians(source_ra), color="red", marker='+', s=75, label="Dwarf Galaxies")
    plt.legend(loc="best")
    plt.tight_layout()
    plt.savefig(save_folder + "skymap_sources.png")
else:
    plt.tight_layout()
    plt.savefig(save_folder + "skymap.png")
plt.close(1)

if MC:

    plt.figure()
    plt.hist(data["ow"], bins=100)
    plt.yscale("log")
    plt.ylabel("Counts")
    plt.xlabel("OneWeight/N_events [GeV cm^2 sr]")
    plt.title("MC OneWeight/N_events Distribution, n=%i"%len(data["ow"]))
    plt.tight_layout()
    plt.savefig(save_folder + "ow.png")
    plt.close()

    if deg:
        sindec_MC = np.sin(data["trueDec"]*np.pi/180)
    else:
        sindec_MC = np.sin(data["trueDec"])

    plt.figure()
    plt.scatter(sindec_MC, data["angErr"], marker=".", s=1/10)
    plt.xlabel("Sin(Declination)")
    plt.ylabel("Angular Error [" + angle_units + ']')
    plt.title("Angular Error vs. MC Sin(Declination)")
    plt.tight_layout()
    plt.savefig(save_folder + "angError_" + angle_units + "_sindec_mc.png")
    plt.close()

    plt.figure()
    scrambled = np.random.uniform(ra_range[0], ra_range[1], len(data["trueRa"]))
    plt.hist2d(scrambled, data["trueDec"], bins=100, norm=matplotlib.colors.LogNorm())
    plt.xlabel("Right Ascension [" + ra_units + ']')
    plt.ylabel("Declination [" + angle_units + ']')
    bar = plt.colorbar()
    bar.set_label("Counts")
    if use_sources and sources_file != "None":
        plt.scatter(source_ra[0], source_dec[0], color="red", marker="+", label="Dwarf Galaxies")
        for i in range(1, len(source_ra)):
            plt.scatter(source_ra[i], source_dec[i], color="red", marker="+", s=75)
        plt.legend(loc="best")
        plt.title("Scrambled MC RA and Declination with Dwarf Galaxies")
        plt.tight_layout()
        plt.savefig(save_folder + "ra_" + ra_units + "_dec_" + angle_units + "_mc_sources.png")
    else:
        plt.title("Scrambled MC RA and Declination")
        plt.tight_layout()
        plt.savefig(save_folder + "ra_" + ra_units + "_dec_" + angle_units + "_mc.png")
    plt.close()

    fig = plt.figure(num=2, figsize=(11,6))
    hpx_map = cat2hpx(scrambled, data["trueDec"], nside=32)
    if use_sources and sources_file != "None":
        hp.mollview(hpx_map, fig=2, coord='C', rot=180, title="Scrambled MC Data with Sources", norm="hist")
    else:
        hp.mollview(hpx_map, fig=2, coord='C', rot=180, title="Scrambled MC Data", norm="hist")
    ax1 = fig.get_axes()[0]
    ax1.text(2.05,0., r"$0^\circ$", ha="left", va="center")
    ax1.text(-2.05,0., r"$360^\circ$", ha="right", va="center")
    ax1.text(1.9,0.45, r"$30^\circ$", ha="left", va="center")
    ax1.text(1.4,0.8, r"$60^\circ$", ha="left", va="center")
    ax1.text(1.9,-0.45, r"$-30^\circ$", ha="left", va="center")
    ax1.text(1.4,-0.8, r"$-60^\circ$", ha="left", va="center")
    hp.graticule()
    if use_sources and sources_file != "None":
        hp.projscatter(np.pi/2-source_dec, source_ra, color="red", marker='+', s=75, label="Dwarf Galaxies")
        plt.legend(loc="best")
        plt.tight_layout()
        plt.savefig(save_folder + "skymap_mc_sources.png")
    else:
        plt.tight_layout()
        plt.savefig(save_folder + "skymap_mc.png")
    plt.close(2)
