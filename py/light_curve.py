#!/anaconda3/envs/py36

# Plotting
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as pl
import seaborn as sns
sns.set_style('ticks')


font = {'family' : 'serif',
        'size'   : 14}
pl.rc('font', **font)

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

    fig, ax1 = pl.subplots(figsize=(4.5, 4.5))


    ax1.errorbar(synphot["dt"], synphot["B"]+2, yerr=[synphot["B_el"], synphot["B_eh"]], fmt=".", label="B+2")
    ax1.errorbar(synphot["dt"]+0.1, synphot["R"], yerr=[synphot["R_el"], synphot["R_eh"]], fmt=".", label="R")
    ax1.errorbar(synphot["dt"]+0.2, synphot["J"]-2, yerr=[synphot["J_el"], synphot["J_eh"]], fmt=".", label="J-2")


    grond_phot = pd.read_csv("../data/grond_phot.csv")

    print(grond_phot)
    ax1.errorbar(grond_phot["dt"], grond_phot["J mag"] - 36.0983 - 2, yerr=grond_phot["J magerr"], fmt=".", label="J-2")
    ax1.errorbar(grond_phot["dt"], grond_phot["H mag"] - 36.0983 - 2, yerr=grond_phot["H magerr"], fmt=".", label="H-2")

    ax1.errorbar(grond_phot["dt"], grond_phot["g' mag"] - 36.0983 + 2, yerr=grond_phot["g' magerr"], fmt=".", label="g'+2'")
    ax1.errorbar(grond_phot["dt"], grond_phot["r' mag"] - 36.0983, yerr=grond_phot["r' magerr"], fmt=".", label="r'")

    # ax1.errorbar(grond_phot["dt"], grond_phot["i' mag"] - 36.0983, yerr=grond_phot["i' magerr"], fmt=".", label="i'")

    # ax1.errorbar(grond_phot["dt"], grond_phot["z' mag"] - 36.0983, yerr=grond_phot["r' magerr"], fmt=".", label="z'")



    ax1.set_xlabel("Time since GRB trigger (days)")
    ax1.set_ylabel("Absolute magnitude (mag)")
    ax1.set_xlim((0, 69))
    ax1.set_ylim((-21, -9,))

    pl.legend()
    ax1.invert_yaxis()
    pl.savefig("../figures/synthetic_lightcurve.pdf")
    pl.show()


if __name__ == '__main__':
    main()
