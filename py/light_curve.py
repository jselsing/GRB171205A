#!/anaconda3/envs/py36

# Plotting
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as pl
import seaborn as sns
sns.set_style('ticks')


font = {'family' : 'serif',
        'size'   : 18}
pl.rc('font', **font)
pl.rc('text', usetex=True)

import numpy as np
import pandas as pd
from astropy.io import fits
from scipy import interpolate
from scipy.signal import medfilt
import pyphot
lib = pyphot.get_library()
import astropy.units as u
import datetime


def main():
    synphot = pd.read_csv("../data/synphot.csv")

    fig, ax1 = pl.subplots(figsize=(5.5, 4.5))
    colors = ["#4C72B0", "#C44E52", "#bd7f7f"]
    ax1.errorbar(synphot["dt"]+0.2, synphot["J"]-2, yerr=[synphot["J_el"], synphot["J_eh"]], fmt=".", label="J-2", color=colors[2])
    ax1.plot(synphot["dt"]+0.2, synphot["J"].values-2, color=colors[2])
    ax1.errorbar(synphot["dt"]+0.1, synphot["R"], yerr=[synphot["R_el"], synphot["R_eh"]], fmt=".", label="R", color=colors[1])
    ax1.plot(synphot["dt"]+0.1, synphot["R"].values, color=colors[1])
    ax1.errorbar(synphot["dt"], synphot["B"]+2, yerr=[synphot["B_el"], synphot["B_eh"]], fmt=".", label="B+2", color=colors[0])
    ax1.plot(synphot["dt"], synphot["B"].values+2, color=colors[0])


    # grond_phot = pd.read_csv("../data/grond_phot.csv")
    # ax1.errorbar(grond_phot["dt"], grond_phot["J mag"].values -36.0983 - 2, yerr=grond_phot["J magerr"], fmt=".", label="J-2", color=colors[2])
    # ax1.plot(grond_phot["dt"], grond_phot["J mag"].values -36.0983 - 2, color=colors[2])
    # ax1.errorbar(grond_phot["dt"], grond_phot["r' mag"].values - 36.0983 , yerr=grond_phot["r' magerr"], fmt=".", label="r'", color=colors[1])
    # ax1.plot(grond_phot["dt"], grond_phot["r' mag"].values - 36.0983, color=colors[1])

    # ax1.errorbar(grond_phot["dt"], grond_phot["g' mag"].values-36.0983 + 2, yerr=grond_phot["g' magerr"], fmt=".", label="g'+2", color=colors[0])
    # ax1.plot(grond_phot["dt"], grond_phot["g' mag"].values-36.0983 + 2, color=colors[0])

    # print(grond_phot)

    # ax1.errorbar(grond_phot["dt"], grond_phot["g' mag"]-36.0983 + 2, yerr=grond_phot["g' magerr"], fmt=".", label="g'", color=colors[0])
    # ax1.errorbar(grond_phot["dt"], grond_phot["r' mag"]-36.0983 , yerr=grond_phot["r' magerr"], fmt=".", label="r'", color=colors[1])
    # ax1.errorbar(grond_phot["dt"], grond_phot["i' mag"]-36.0983, yerr=grond_phot["i' magerr"], fmt=".", label="i'")
    # ax1.errorbar(grond_phot["dt"], grond_phot["z' mag"]-36.0983, yerr=grond_phot["r' magerr"], fmt=".", label="z'")
    # ax1.errorbar(grond_phot["dt"], grond_phot["J mag"] -36.0983 - 2, yerr=grond_phot["J magerr"], fmt=".", label="J", color=colors[2])
    # ax1.errorbar(grond_phot["dt"], grond_phot["H mag"] -36.0983, yerr=grond_phot["H magerr"], fmt=".", label="H")



    ax1.set_xlabel("Time since GRB trigger (days)")
    ax1.set_ylabel("Absolute magnitude (mag)")
    ax1.set_xlim((0, 69))
    ax1.set_ylim((-21.5, -8.5,))
    ax1.set_yticks([-10, -12, -14, -16, -18, -20])
    ax1.set_yticklabels([-10, -12, -14, -16, -18, -20])

    pl.legend()
    ax1.invert_yaxis()
    pl.tight_layout()
    pl.savefig("../figures/synthetic_lightcurve.pdf")
    pl.show()


if __name__ == '__main__':
    main()
