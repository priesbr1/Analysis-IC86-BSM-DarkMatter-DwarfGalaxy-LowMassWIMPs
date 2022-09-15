import os, sys
import math
import numpy as np
import argparse
from DM_FluxComputation import *

import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
from matplotlib import colors
from mpl_toolkits.axes_grid1 import make_axes_locatable

parser = argparse.ArgumentParser()
parser.add_argument("-o", "--oscillated", type=int, default=True,
                    dest="oscillated", help="whether or not neutrino spectra have already been oscillated")
parser.add_argument("-p", "--output", type=str, default=None,
                    dest="output", help="output folder for plots")
parser.add_argument("-u", "--flux_units", type=int, default=0,
                    dest="units", help="whether or not to use flux units of GeV^-1 (0) or GeV^-1 cm^-2 s^-1 (1)")
args = parser.parse_args()

def oscillate_PPPC4(spectra, channels, masses, nu_types):
    
    oscillated = dict()
    
    for channel in channels:
        print("Doing channel %s"%channel)
        oscillated[channel] = dict()

        for mass in masses:
            oscillated[channel][mass] = dict()

            # Oscillate neutrinos
            osc_nu = oscillate_spectra(spectra[channel][mass], ["nu_mu", "nu_e", "nu_tau"])

            # Spectra after oscillation
            for var in osc_nu.keys():
                oscillated[channel][mass][var] = np.array(osc_nu[var])
            
    return oscillated

if (args.oscillated == False):
    nu_e = np.load("/mnt/home/priesbr1/DM_Search/DM_CirelliSpectrum_dict_neutrinos_e.npy", allow_pickle=True)
    nu_mu = np.load("/mnt/home/priesbr1/DM_Search/DM_CirelliSpectrum_dict_neutrinos_mu.npy", allow_pickle=True)
    nu_tau = np.load("/mnt/home/priesbr1/DM_Search/DM_CirelliSpectrum_dict_neutrinos_tau.npy", allow_pickle=True)
    
    nu_e = nu_e.item()
    nu_mu = nu_mu.item()
    nu_tau = nu_tau.item()
    
    channels = list(nu_e.keys())
    masses = list(nu_e[channels[0]].keys())
    
    osc_nu_e = oscillate_PPPC4(nu_e, channels=channels, masses=masses, nu_types=["nu_mu", "nu_e", "nu_tau"])
    osc_nu_mu = oscillate_PPPC4(nu_mu, channels=channels, masses=masses, nu_types=["nu_mu", "nu_e", "nu_tau"])
    osc_nu_tau = oscillate_PPPC4(nu_tau, channels=channels, masses=masses, nu_types=["nu_mu", "nu_e", "nu_tau"])
    
    nu_tot = dict()
    osc_nu_tot = dict()
    
    for channel in nu_e.keys():
        nu_tot[channel] = dict()
        osc_nu_tot[channel] = dict()
        
        for mass in nu_e[channel].keys():
            nu_tot[channel][mass] = dict()
            osc_nu_tot[channel][mass] = dict()
            
            for var in nu_e[channel][mass].keys():
                nu_tot[channel][mass][var] = (nu_e[channel][mass][var] + nu_mu[channel][mass][var] + nu_tau[channel][mass][var])/3
                if "dNdE" in var:
                    osc_nue_var = (osc_nu_e[channel][mass][var+"_nu_e_osc"] + osc_nu_mu[channel][mass][var+"_nu_e_osc"] + osc_nu_tau[channel][mass][var+"_nu_e_osc"])/3
                    osc_numu_var = (osc_nu_e[channel][mass][var+"_nu_mu_osc"] + osc_nu_mu[channel][mass][var+"_nu_mu_osc"] + osc_nu_tau[channel][mass][var+"_nu_mu_osc"])/3
                    osc_nutau_var = (osc_nu_e[channel][mass][var+"_nu_tau_osc"] + osc_nu_mu[channel][mass][var+"_nu_tau_osc"] + osc_nu_tau[channel][mass][var+"_nu_tau_osc"])/3
                    osc_nu_tot[channel][mass][var] = (osc_nue_var + osc_numu_var + osc_nutau_var)/3
                else:
                    osc_nu_tot[channel][mass][var] = (osc_nu_e[channel][mass][var] + osc_nu_mu[channel][mass][var] + osc_nu_tau[channel][mass][var])/3
    
    np.save("/mnt/home/priesbr1/DM_Search/DM_CirelliSpectrum_dict_neutrinos_all.npy", nu_tot)
    np.save("/mnt/home/priesbr1/DM_Search/DM_CirelliSpectrum_dict_neutrinos_oscillated.npy", osc_nu_tot)

else:
    nu_tot = np.load("/mnt/home/priesbr1/DM_Search/DM_CirelliSpectrum_dict_neutrinos_all.npy", allow_pickle=True)
    osc_nu_tot = np.load("/mnt/home/priesbr1/DM_Search/DM_CirelliSpectrum_dict_neutrinos_oscillated.npy", allow_pickle=True)
    
    nu_tot = nu_tot.item()
    osc_nu_tot = osc_nu_tot.item()

plot_masses = [5, 10, 20, 30, 50, 80, 90, 100, 150, 200, 300, 500, 600, 1000, 9000]
plot_channels = ["b", "Mu", "Tau", "W", "Nu"] #plot_channels = ["b", "Mu", "Tau", "W", "Nu", "Nue", "NuMu", "NuTau"]
legends = {"b":r"$b\bar{b}$", "Mu":r"$\mu^+\mu^-$", "Tau":r"$\tau^+\tau^-$", "W":r"$W^+W^-$", "Nu":r"$\nu\bar{\nu}$", "Nue":r"$\nu_e\bar{\nu_e}$", "NuMu":r"$\nu_{\mu}\bar{\nu_{\mu}}$", 
           "NuTau":r"$\nu_{\tau}\bar{\nu_{\tau}}$"}
colors = {"b":"mediumturquoise", "Mu":"forestgreen", "Tau":"dodgerblue", "W":"darkmagenta", "Nu":"mediumvioletred", "Nue":"darkorchid", "NuMu":"orchid", "NuTau":"lightpink"}

nu_type = "nu_mu"
type_legend = {"nu_e":r"$\nu_\mathrm{e}$", "nu_mu":r"$\nu_\mathrm{\mu}$", "nu_tau":r"$\nu_\mathrm{\tau$}"}

for mass in plot_masses:
    fig = plt.figure()
    ax = fig.add_subplot(111)
    fig.subplots_adjust(top=0.85)
    fig.patch.set_facecolor("white")
    for channel in plot_channels:
        if (channel == "W" and mass <100) or (channel == "b" and mass <10):
            continue
        
        ax.plot(osc_nu_tot[channel][mass]["Energy"], osc_nu_tot[channel][mass]["EdNdE"]/osc_nu_tot[channel][mass]["Energy"], color=colors[channel], linestyle="-", label=legends[channel], 
                linewidth=2)
        ax.plot(nu_tot[channel][mass]["Energy"], nu_tot[channel][mass]["EdNdE"]/nu_tot[channel][mass]["Energy"], color=colors[channel], linestyle="--", linewidth=2)
        ax.semilogy()
        ax.set_xlabel(r"$E_\mathrm{\nu}$ [GeV]", fontsize=12)
        ax.set_ylabel(r"$\mathrm{d}N/\mathrm{d}E_\mathrm{\nu}$ [GeV$^{-1}$]", fontsize=12)
        ax.legend(loc="best", ncol=2, prop={"size":12}, labelspacing=0.5)
        ax.grid(color="grey", linestyle="-", linewidth=0.5, alpha=0.5)
    plt.savefig(args.output + "spectra_energy_WIMP_%i.png"%mass)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    fig.subplots_adjust(top=0.85)
    fig.patch.set_facecolor("white")
    for channel in plot_channels:
        if (channel == "W" and mass <100) or (channel == "b" and mass <10):
            continue

        ax.plot(np.log10(osc_nu_tot[channel][mass]["Energy"]/mass), osc_nu_tot[channel][mass]["EdNdE"], color=colors[channel], linestyle="-", label=legends[channel], linewidth=2)
        ax.plot(np.log10(nu_tot[channel][mass]["Energy"]/mass), nu_tot[channel][mass]["EdNdE"], color=colors[channel], linestyle="--", linewidth=2)
        ax.semilogy()
        ax.set_xlabel(r"$\log_\mathrm{10}\left(E_\mathrm{\nu}/M_\mathrm{WIMP}\right)$", fontsize=12)
        ax.set_ylabel(r"$E \mathrm{d}N/\mathrm{d}E_\mathrm{\nu}$", fontsize=12)
        ax.legend(loc="best", ncol=2, prop={"size":12}, labelspacing=0.5)
        ax.grid(color="grey", linestyle="-", linewidth=0.5, alpha=0.5)
    plt.savefig(args.output + "spectra_logX_WIMP_%i.png"%mass)

if (args.units == True):
    mc = np.load("/mnt/research/IceCube/datasets/ps_DRAGON/IC86_2013_MC.npy")
    livetime = 370.9925128750838 # days
    
    # Calculate effective area
    decs = np.radians(np.asarray([[-90., 90.]]))
    sin_decs = np.sin(decs)
    E_bins = np.logspace(0., 9., 61)
    logE_bins = np.log10(E_bins)
    dlog_E = np.diff(logE_bins)
    
    for ii, (low_dec, high_dec) in enumerate(sin_decs):
        d_omega = 2.*np.pi*np.abs(high_dec - low_dec)
        dec_msk = np.sin(mc['trueDec']) > low_dec
        dec_msk *= np.sin(mc['trueDec']) < high_dec
        mc_cut = mc[dec_msk]
        weights = mc_cut['ow'] / (1e4 * mc_cut['trueE'] * dlog_E[np.digitize(np.log10(mc_cut['trueE']), bins = logE_bins) -1] * d_omega * np.log(10.))
        A_eff, bins = np.histogram(mc_cut['trueE'], bins=E_bins, weights=weights)
        print(bins)
        A_eff /= 100**2 # Convert m^2 to cm^2

    # Create temporary arrays of energy/flux for channels/masses
    data_tot = dict()
    data_osc = dict()
    
    for channel in nu_tot.keys():
        data_tot[channel] = dict()
        for mass in nu_tot[channel].keys():
            data_tot[channel][mass] = dict()
            data_tot[channel][mass]["Energy"] = np.array(nu_tot[channel][mass]["Energy"])
            data_tot[channel][mass]["dNdE"] = np.array(nu_tot[channel][mass]["EdNdE"] / nu_tot[channel][mass]["Energy"])
            
            for i, energy in enumerate(data_tot[channel][mass]["Energy"]):
                for j, A in enumerate(A_eff):
                    if (bins[j] < energy) and (energy < bins[j+1]):
                        data_tot[channel][mass]["dNdE"][i] /= A
            data_tot[channel][mass]["dNdE"] /= (livetime * 86400) # Convert days to s
    
    for channel in osc_nu_tot.keys():
        data_osc[channel] = dict()
        for mass in osc_nu_tot[channel].keys():
            data_osc[channel][mass] = dict()
            data_osc[channel][mass]["Energy"] = np.array(osc_nu_tot[channel][mass]["Energy"])
            data_osc[channel][mass]["dNdE"] = np.array(osc_nu_tot[channel][mass]["EdNdE"] / osc_nu_tot[channel][mass]["Energy"])
    
            for i, energy in enumerate(data_osc[channel][mass]["Energy"]):
                for j, A in enumerate(A_eff):
                    if (bins[j] < energy) and (energy < bins[j+1]):
                        data_osc[channel][mass]["dNdE"][i] /= A
            data_osc[channel][mass]["dNdE"] /= (livetime * 86400)
    
    for mass in plot_masses:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        fig.subplots_adjust(top=0.85)
        fig.patch.set_facecolor("white")
        for channel in plot_channels:
            if (channel == "W" and mass <100) or (channel == "b" and mass <10):
                continue
            
            ax.plot(data_osc[channel][mass]["Energy"], data_osc[channel][mass]["dNdE"], color=colors[channel], linestyle="-", label=legends[channel], linewidth=2)
            ax.plot(data_tot[channel][mass]["Energy"], data_tot[channel][mass]["dNdE"], color=colors[channel], linestyle="--", linewidth=2)
            ax.semilogy()
            ax.set_xlabel(r"$E_\mathrm{\nu}$ [GeV]", fontsize=12)
            ax.set_ylabel(r"$\mathrm{d}N/\mathrm{d}E_\mathrm{\nu}$ [GeV$^{-1}$ cm$^{-2}$ s$^{-1}$]", fontsize=12)
            ax.legend(loc="best", ncol=2, prop={"size":12}, labelspacing=0.5)
            ax.grid(color="grey", linestyle="-", linewidth=0.5, alpha=0.5)
        plt.savefig(args.output + "spectra_energy_flux_WIMP_%i.png"%mass)
        
        fig = plt.figure()
        ax = fig.add_subplot(111)
        fig.subplots_adjust(top=0.85)
        fig.patch.set_facecolor("white")
        for channel in plot_channels:
            if (channel == "W" and mass <100) or (channel == "b" and mass <10):
                continue
            
            ax.plot(np.log10(data_osc[channel][mass]["Energy"]/mass), data_osc[channel][mass]["dNdE"]*data_osc[channel][mass]["Energy"], color=colors[channel], linestyle="-", 
                    label=legends[channel], linewidth=2)
            ax.plot(np.log10(data_tot[channel][mass]["Energy"]/mass), data_tot[channel][mass]["dNdE"]*data_tot[channel][mass]["Energy"], color=colors[channel], linestyle="--", linewidth=2)
            ax.semilogy()
            ax.set_xlabel(r"$\log_\mathrm{10}\left(E_\mathrm{\nu}/M_\mathrm{WIMP}\right)$", fontsize=12)
            ax.set_ylabel(r"$E \mathrm{d}N/\mathrm{d}E_\mathrm{\nu}$ [cm$^{-2}$ s$^{-1}$]", fontsize=12)
            ax.legend(loc="best", ncol=2, prop={"size":12}, labelspacing=0.5)
            ax.grid(color="grey", linestyle="-", linewidth=0.5, alpha=0.5)
        plt.savefig(args.output + "spectra_logX_flux_WIMP_%i.png"%mass)
