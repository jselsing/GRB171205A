#!/anaconda3/envs/py36

# Plotting
import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as pl
import seaborn as sns
sns.set_style('ticks')

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

    fig, ax1 = pl.subplots()


    ax1.errorbar(synphot["dt"], synphot["B"]+2, yerr=[synphot["B_el"], synphot["B_eh"]], fmt=".", label="B+2")
    ax1.errorbar(synphot["dt"]+0.1, synphot["R"], yerr=[synphot["R_el"], synphot["R_eh"]], fmt=".", label="R")
    ax1.errorbar(synphot["dt"]+0.2, synphot["J"]-2, yerr=[synphot["J_el"], synphot["J_eh"]], fmt=".", label="J-2")
    ax1.set_xlabel("Time since GRB trigger (days)")
    ax1.set_ylabel("Absolute magnitude (mag)")
    ax1.set_xlim((0, 279))
    ax1.set_ylim((-21, -9,))

    pl.legend()
    ax1.invert_yaxis()
    pl.savefig("../figures/synthetic_lightcurve.pdf")
    pl.show()


if __name__ == '__main__':
    main()
