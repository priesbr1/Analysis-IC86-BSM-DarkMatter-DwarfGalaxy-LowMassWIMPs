import numpy as np
import matplotlib as mpl
mpl.use("Agg")
mpl.rc("font", family="serif", size="15")
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import argparse
import pickle

parser = argparse.ArgumentParser()
parser.add_argument("-d", "--datafolder", type=str, default="/mnt/home/priesbr1/DM_Search/data/energy_optimization/",
                    dest="datafolder", help="folder containing results from energy optimization calculations")
parser.add_argument("-c", "--channel", type=str, default="b",
                    dest="channel", help="annihilation channel to create plots for")
parser.add_argument("-m", "--masses", type=int, nargs="+", default=[10,20,30,50,80,100,150,200,300],
                    dest="masses", help="WIMP masses to create plots for")
parser.add_argument("-o", "--outfolder", type=str, default="/mnt/scratch/priesbr1/DM_Search/Plots/energy_optimization/",
                    dest="outfolder", help="output folder for plots")
args = parser.parse_args()

result_list = []

for mass in args.masses:
    print(args.channel, mass)
    signal_list = []
    background_list = []
    for i in range(0,10):
        filename = args.datafolder + "weight_vs_cutoff_" + args.channel + "_" + str(mass) + "_part" + str(i) + ".txt"
        try:
            f = open(filename,"rb")
        except:
            print("Cannot find %s -- skipping..."%filename)
            continue
        garbage = pickle.load(f)
        garbage = pickle.load(f)
        signal_list.append(pickle.load(f))
        background_list.append(pickle.load(f))
        f.close()
    
    for i in range(0,10):
        counter = 0
        for lower_limit in range(0,500):
            for upper_limit in range(lower_limit+1,500):
                signal_list[i][counter] *= 0.1
                background_list[i][counter] *= 0.1
                counter += 1
	
    for i in range(1,10):
        counter = 0
        for lower_limit in range(0,500):
            for upper_limit in range(lower_limit+1,500):
                signal_list[0][counter] += signal_list[i][counter]
                background_list[0][counter] += background_list[i][counter]
                counter += 1
    
    ratio = []
    ll = []
    ul = []
    r_full = []
    
    counter = 0
    for lower_limit in range(0,500):
        for upper_limit in range(lower_limit+1,500):
            if (background_list[0][counter] > 0):
                ratio.append(signal_list[0][counter]/np.sqrt(background_list[0][counter]))
            else:
                ratio.append(1e-38)
            counter += 1
    
    counter = 0
    for lower_limit in range(0,500):
        for upper_limit in range(0,500):
            ll.append(lower_limit)
            ul.append(upper_limit)
            if (upper_limit > lower_limit) and (background_list[0][counter] > 0):
                r_full.append(signal_list[0][counter]/np.sqrt(background_list[0][counter]))
                counter += 1
            else:
                r_full.append(0)
     
    max_ratio = max(ratio)
    max_index = ratio.index(max_ratio)
    
    counter = 0
    for lower_limit in range(0,500):
        for upper_limit in range(lower_limit+1,500):
            if (counter == max_index):
                le = lower_limit
                he = upper_limit
                print(lower_limit,upper_limit)
            counter += 1
        
    data = np.empty(len(r_full), dtype=[("lo",float),("hi",float),("ratio",float)])
    data["lo"] = ll
    data["hi"] = ul
    data["ratio"] = r_full

    print(len(data[data["ratio"] != 0]))

    s = plt.scatter(data["lo"], data["hi"], c=data["ratio"], norm=colors.Normalize(vmin=min(ratio)+1e-30, vmax=max(ratio)), cmap="magma")
    plt.plot(le, he, marker="x", linestyle="None", markersize=12, label="Optimum")

    #plt.xlim([0,2*mass])
    #plt.ylim([0,30])
    cbar = plt.colorbar(s)
    cbar.set_label("Weights")
    plt.xlabel(r"Lower Cut [$\mathrm{GeV}$]")
    plt.ylabel(r"High Cut [$\mathrm{GeV}$]")
    plt.legend()
    plt.tight_layout()
    plt.savefig(args.outfolder + "optimal_energy_range_" + args.channel + "_" + str(mass) + ".png")
    plt.close()
    #plt.show()
    #exit()
    #print("optimal", np.sqrt(background_list[0][max_index]/background_list[0][mass-1])/(signal_list[0][max_index]/signal_list[0][mass-1]))
    #print("optimal index", max_index)
    #print("optimal signal", signal_list[0][max_index])
    #print("current signal", signal_list[0][mass-1])
    #print("result", signal_list[0][mass-1]/(np.sqrt(background_list[0][max_index]/background_list[0][mass-1])/(signal_list[0][max_index]/signal_list[0][mass-1])))
    #result_list.append(signal_list[0][max_index])
