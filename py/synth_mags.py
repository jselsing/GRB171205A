#!/anaconda3/envs/py36


# Plotting
import matplotlib; matplotlib.use('Qt5Agg')
import matplotlib.pyplot as pl
import seaborn as sns; sns.set_style('ticks')

import numpy as np
import pandas as pd
from astropy.io import fits
from scipy import interpolate
from scipy.signal import medfilt
import pyphot
lib = pyphot.get_library()
import astropy.units as u
import datetime


def timediff(t1, t2):

    time1 = datetime.datetime(int(t1[0:4]), int(t1[5:7]), int(t1[8:10]), hour=int(round(float(t1[11:13]))), minute=int(round(float(t1[14:16]))), second=int(round(float(t1[17:19]))))
    time2 = datetime.datetime(int(t2[0:4]), int(t2[5:7]), int(t2[8:10]), hour=int(round(float(t2[11:13]))), minute=int(round(float(t2[14:16]))), second=int(round(float(t2[17:19]))))

    dt = time2 - time1
    return dt.total_seconds()/(3600)/(24)


def get_data(path, OB, filter_bad = False):

    fac = 10
    # print(int(OB[2:]))
    # if int(OB[2:]) > 5:
    #     fac = 1


    f = fits.open("%s/UVB%s.fits"%(path, OB))
    # print(f[1].columns)
    wl = fac*f[1].data.field("WAVE").flatten()
    q = f[1].data.field("QUAL").flatten()
    mask_wl = (wl > 3200) & (wl < 5600)
    mask_qual = ~q.astype("bool")
    flux = interpolate.interp1d(wl[mask_qual], f[1].data.field("FLUX").flatten()[mask_qual], bounds_error=False, fill_value=0)
    error = interpolate.interp1d(wl[mask_qual], f[1].data.field("ERR").flatten()[mask_qual], bounds_error=False, fill_value=0)
    wl_plot = wl[mask_wl]
    flux = flux(wl_plot)
    error = error(wl_plot)
    wl_UVB = wl_plot
    flux_UVB = flux
    error_UVB = error

    f = fits.open("%s/VIS%s.fits"%(path, OB))
    wl = fac*f[1].data.field("WAVE").flatten()
    q = f[1].data.field("QUAL").flatten()
    try:
        t = f[1].data.field("TRANS").flatten()
    except:
        t = np.ones_like(f[1].data.field("QUAL").flatten())
    mask_wl = (wl > 5600) & (wl < 10200)
    mask_qual = ~q.astype("bool")

    flux = interpolate.interp1d(wl[mask_qual], f[1].data.field("FLUX").flatten()[mask_qual], bounds_error=False, fill_value=0)
    error = interpolate.interp1d(wl[mask_qual], f[1].data.field("ERR").flatten()[mask_qual], bounds_error=False, fill_value=0)
    wl_plot = wl[mask_wl]
    flux = flux(wl_plot)/t[mask_wl]
    error = error(wl_plot)/t[mask_wl]
    wl_VIS = wl_plot
    flux_VIS = flux
    error_VIS = error

    f = fits.open("%s/NIR%s.fits"%(path, OB))
    try:
        f_nod = fits.open("%s/NIR%s_NOD.fits"%(path, OB))
    except:
        f_nod = fits.open("%s/NIR%s.fits"%(path, OB))
    wl = fac*f[1].data.field("WAVE").flatten()
    f[1].data["FLUX"] = np.concatenate((f[1].data["FLUX"].flatten()[(wl <= 22500)], f_nod[1].data["FLUX"].flatten()[(wl > 22500)]))
    f[1].data["ERR"] = np.concatenate((f[1].data["ERR"].flatten()[(wl <= 22500)], f_nod[1].data["ERR"].flatten()[(wl > 22500)]))

    q = f[1].data.field("QUAL").flatten()
    try:
        t = f[1].data.field("TRANS").flatten()
    except:
        t = np.ones_like(f[1].data.field("QUAL").flatten())
    mask_wl = (wl > 10200)
    mask_qual = ~q.astype("bool")
    flux = interpolate.interp1d(wl[mask_qual], f[1].data.field("FLUX").flatten()[mask_qual], bounds_error=False, fill_value=0)
    error = interpolate.interp1d(wl[mask_qual], f[1].data.field("ERR").flatten()[mask_qual], bounds_error=False, fill_value=0)
    wl_plot = wl[mask_wl]
    flux = flux(wl_plot)/t[mask_wl]
    error = error(wl_plot)/t[mask_wl]
    wl_NIR = wl_plot
    flux_NIR = flux
    error_NIR = error



    wl = np.concatenate((wl_UVB, wl_VIS, wl_NIR))
    flux = np.concatenate((flux_UVB, flux_VIS, flux_NIR))
    error = np.concatenate((error_UVB, error_VIS, error_NIR))

    if filter_bad:

        wl, flux, error = filter_bad_values(wl, flux, error)

    return wl, flux, error, f[0].header["DATE"]


def filter_bad_values(wl, flux, error):
    medfilter = medfilt(flux, 501)
    mask = np.logical_and(abs(flux - medfilter) < 3*error, ~np.isnan(flux))
    f = interpolate.interp1d(wl[mask], flux[mask], bounds_error=False, fill_value = np.nan)
    g = interpolate.interp1d(wl[mask], error[mask], bounds_error=False, fill_value = np.nan)
    return wl, f(wl), g(wl)

def synpassflux(wl, flux, band):
    # Calculate synthetic magnitudes
    n_bands = len(band)

    filt = np.genfromtxt("/Users/jonatanselsing/github/iPTF16geu/data/passbands/%s"%band)
    lamb_T, T = filt[:,0], filt[:,1]
    f = pyphot.Filter(lamb_T, T, name=band, dtype='photon', unit='Angstrom')
    fluxes = f.get_flux(wl, flux, axis=0)
    synmags = -2.5 * np.log10(fluxes) - f.AB_zero_mag
    cwav = np.mean(lamb_T)
    cwave = (float(max(lamb_T[T > np.percentile(T, 10)] - cwav)), float(cwav - min(lamb_T[T > np.percentile(T, 10)])))
    synmag_flux = ((synmags*u.ABmag).to((u.erg/(u.s * u.cm**2 * u.AA)), u.spectral_density(cwav * u.AA))).value
    return synmag_flux, cwav, cwave, synmags



# GRB121027A
def main():

    z = 0.0368
    OBs = ["OB1", "OB2", "OB3", "OB4", "OB5", "OB6", "OB7"]

    ra = 167.447
    dec = -12.604
    t_trigger = "2017:12:05T07:20:43"

    dist_mod = 36.0983


    fig, ax1 = pl.subplots()


    for ii, OB in enumerate(OBs):

        wl, flux, error, t_obs = get_data("../data", OB, filter_bad=True)
        dt = np.around(timediff(t_trigger, t_obs), 1)


        passbands = ["Bessell_B.dat", "Bessell_R.dat", "2MASS_J.dat"]
        pass_names = ["B", "R", "J"]

        # passbands = ["SDSS_u.dat", "SDSS_g.dat", "SDSS_r.dat", "SDSS_i.dat", "SDSS_z.dat", "2MASS_J.dat", "2MASS_H.dat", "2MASS_Ks.dat"]
        # pass_names = ["u'", "g'", "r'", "i'", "z'", "J", "H", "Ks"]



        # bands = photometry.loc[photometry['GRB'] == row["GRB"]]
        # print(bands)
        # exit()
        pl.plot(wl[::10], medfilt(flux, 51)[::10], color="black", alpha= 0.7, linestyle="steps-mid")
        for pp, kk in enumerate(passbands):
            # print(bands['INS.FILT1.NAME'])
            # meas_mag = bands.loc[bands['INS.FILT1.NAME'] == "%s_prime"%kk[-5]]







            synmag_flux, cwav, cwave, synmag = synpassflux(wl, flux, kk)
            synmag_error, cwav, cwave, synmag_err = synpassflux(wl, error, kk)
            _, _, _, synmag_err_up = synpassflux(wl, flux + error, kk)
            _, _, _, synmag_err_low = synpassflux(wl, flux - error, kk)
            e_u = (synmag - synmag_err_up)
            e_l = (synmag_err_low - synmag)
            pl.errorbar(cwav, synmag_flux, xerr = [[cwave]], yerr = synmag_error,  fmt = 'o', zorder = 10, ms = 5, elinewidth=1.7, label = "%s = %s$^{+%s}_{-%s}$"%(pass_names[pp], np.around(synmag - dist_mod, 1), np.around(e_u, 1), np.around(e_l, 1)))

            # if meas_mag.shape[0] > 0:
            #     meas_flux = ((meas_mag["MAG"].values*u.ABmag).to((u.erg/(u.s * u.cm**2 * u.AA)), u.spectral_density(cwav * u.AA))).value
            #     meas_flux_up =  (((meas_mag["MAG"].values + meas_mag["ERROR"].values)*u.ABmag).to((u.erg/(u.s * u.cm**2 * u.AA)), u.spectral_density(cwav * u.AA))).value
            #     meas_flux_do =  (((meas_mag["MAG"].values - meas_mag["ERROR"].values)*u.ABmag).to((u.erg/(u.s * u.cm**2 * u.AA)), u.spectral_density(cwav * u.AA))).value

            #     pl.errorbar(cwav, meas_flux, xerr = [[cwave]], yerr = [[(meas_flux_do, meas_flux_up)]],  fmt = 'o', ms = 15)




        pl.title("dt = %s days"%dt)
        pl.legend()
        pl.ylabel(r'$\mathrm{F}_{\lambda}~[10^{-17}~\mathrm{erg}~\mathrm{s}^{-1}~\mathrm{cm}^{-2}~\mathrm{\AA}^{-1}]$')
        pl.xlabel(r"Observed Wavelength ($\mathrm{\AA}$)")
        print(np.percentile(medfilt(flux, 51), 90)*2)
        pl.ylim((-1e-18, np.nanpercentile(medfilt(flux, 51), 90)*2))
        pl.xlim(3000, 25000)
        pl.savefig("../figures/SN2017iuk_%s.pdf"%dt)
        pl.clf()

            # with open("density_prof_%s.dat"%OB, "wb") as f:
            #   f.write(b'1e4 second\n')
            #   np.savetxt(f, list(zip(np.arange(len(v)), v, rho.value)), fmt="%i %10.3f %1.4e")


if __name__ == '__main__':
    main()
