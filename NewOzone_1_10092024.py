import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline, interp1d
from scipy.constants import h, c, sigma
from math import e
from scipy.integrate import simpson
import time
import os
from datetime import datetime


def sun_temp(gyears: float) -> float:
    return 5621 + 107 * gyears + (-51.4) * gyears ** 2 + 14.7 * gyears ** 3 + (
        -2.06) * gyears ** 4 + 0.127 * gyears ** 5 + \
        (-2.59 * 10 ** (-3)) * gyears ** 6


def sun_radius(gyears: float) -> float:
    return (0.87 + 0.0299 * gyears + (-1.91) * 10 ** (-3) * gyears ** 2 + 3.5 * 10 ** (-4) * gyears ** 3) \
        * 695700.6371 * 10 ** 3


def flux_bolometric_calc(gyears):
    return ((sun_radius(gyears) ** 2) / ((1.496 * 10 ** 11) ** 2)) * sigma * sun_temp(gyears) ** 4


def smoothing_spline_interpolation_solar(age, age_lst, flux_lst, s):
    interpolated_flux = np.zeros(flux_lst.shape[1])
    for i in range(flux_lst.shape[1]):
        spline = UnivariateSpline(age_lst, flux_lst[:, i], s=s)
        interpolated_flux[i] = spline(age)
    return interpolated_flux


def solar_to_atmosphere_flux(flux_given, tau_lst):
    flux_post_atm = []
    for m in range(len(flux_given)):
        flux_post_atm_item = flux_given[m] * e ** ((-1) * tau_lst[m])
        if flux_post_atm_item >= 10 ** (-5):
            flux_post_atm.append(flux_post_atm_item)
        else:
            flux_post_atm.append(0.0)
    return flux_post_atm


def solar_to_atmosphere_flux_photons(flux_given, tau_lst):
    flux_post_atm = []
    for m in range(len(flux_given)):
        flux_post_atm_item = flux_given[m] * e ** ((-1) * tau_lst[m])
        if flux_post_atm_item >= 1000:
            flux_post_atm.append(flux_post_atm_item)
        else:
            flux_post_atm.append(0.0)
    return flux_post_atm


def atmosphere_to_water_flux(flux_given, coeff_lst, depth):
    flux_after_wat = flux_given * e ** (coeff_lst * depth * (-1))
    for m in range(len(flux_after_wat)):
        if flux_after_wat[m] <= 10 ** (-5):
            flux_after_wat[m] = 0.0
        if depth < 0:
            print(f"Problem with depth: {depth}")
        if coeff_lst[m] < 0:
            print(f"Problem with coefficient: {coeff_lst[m]}")
    return flux_after_wat


def atmosphere_to_mixed_water_flux(flux_given, coeff_lst, depth):
    flux_after_wat = ((flux_given) / (coeff_lst * depth)) * (1 - e ** (coeff_lst * depth * (-1)))
    for m in range(len(flux_after_wat)):
        if flux_after_wat[m] <= 10 ** (-5):
            flux_after_wat[m] = 0.0
        if depth < 0:
            print(f"Problem with depth: {depth}")
        if coeff_lst[m] < 0:
            print(f"Problem with coefficient: {coeff_lst[m]}")
    return flux_after_wat


def atmosphere_to_water_flux_photons(flux_given_photons, coeff_lst, depth):
    flux_after_wat_photons = flux_given_photons * e ** (coeff_lst * depth * (-1))
    for m in range(len(flux_after_wat_photons)):
        if flux_after_wat_photons[m] <= 1000:
            flux_after_wat_photons[m] = 0.0
        if depth < 0:
            print(f"Problem with depth: {depth}")
        if coeff_lst[m] < 0:
            print(f"Problem with coefficient: {coeff_lst[m]}")
    return flux_after_wat_photons


def atmosphere_to_mixed_water_flux_photons(flux_given_photons, coeff_lst, depth):
    flux_after_wat_photons = ((flux_given_photons) / (coeff_lst * depth)) * (1 - e ** (coeff_lst * depth * (-1)))
    for m in range(len(flux_after_wat_photons)):
        if flux_after_wat_photons[m] <= 1000:
            flux_after_wat_photons[m] = 0.0
        if depth < 0:
            print(f"Problem with depth: {depth}")
        if coeff_lst[m] < 0:
            print(f"Problem with coefficient: {coeff_lst[m]}")
    return flux_after_wat_photons


# region Intro Print, Input Statements, Time Start
print("Note: UVC is 100nm to 280nm; UVB is 280nm to 315nm; UVA is 315nm to 400nm.")
print("Disclaimer: For photons units, numbers less than 1000 will be considered negligible (rounded to 0).")
print("Disclaimer: The bond energy of the C-N bond is for an amide (peptide) bond, not a normal C-N bond.")
print("Note: Taking average mass of 1 amino acid as 110 Da.")
print("Note: Taking an average of 300 amino acids per each chain.")

double_CO_bond_wavelength = 161.0  # in nm
double_CC_bond_wavelength = 199.0  # in nm
# single_CP_bond_wavelength = 233.0  # in nm
single_CH_bond_wavelength = 291.0  # in nm
single_CO_bond_wavelength = 334.0  # in nm
single_CC_bond_wavelength = 346.0  # in nm
single_PO_bond_wavelength = 357.0  # in nm
single_CN_bond_wavelength = 392.0  # in nm

print("Wavelengths needed to break each bond below")
print("C=O: 161 nm \nC=C: 195 nm \nC-H: 290 nm \nC-O: 334 nm \nC-C: 345 nm \nP-O: 357 nm\nC-N: 392 nm ")

target_age = float(input("Age of the Sun (in Gyr): "))
# gyears_ago = 4.603 - target_age
bond_break = input("Do you want to break a bond? \nIf no, skip; if yes, what bond do you want to break (ex. for "
                   "typing: C-C, C=O): ").upper()

wavelength_decision = input("Do you want the flux interval relating to bond (y/n) (skip if skipped previous question):"
                            " ").lower()
wavelength_integration_start = 100.0
wavelength_integration_end = 100.0
if "y" in wavelength_decision:
    if bond_break == "C=O":
        wavelength_integration_end = double_CO_bond_wavelength
    elif bond_break == "C=C":
        wavelength_integration_end = double_CC_bond_wavelength
    # elif bond_break == "C-P":
    #     wavelength_integration_end = single_CP_bond_wavelength
    elif bond_break == "C-H":
        wavelength_integration_end = single_CH_bond_wavelength
    elif bond_break == "C-O":
        wavelength_integration_end = single_CO_bond_wavelength
    elif bond_break == "C-C":
        wavelength_integration_end = single_CC_bond_wavelength
    elif bond_break == "P-O":
        wavelength_integration_end = single_PO_bond_wavelength
    elif bond_break == "C-N":
        wavelength_integration_end = single_CN_bond_wavelength
    else:
        bond_break = ""
else:
    wavelength_integration_start = float(input("Wavelength Integration Start: "))
    wavelength_integration_end = float(input("Wavelength Integration End: "))
# fvw_decision = input("Do you want a plot of Flux vs Wavelength (y/n): ").lower()
# fvt_decision = input("Do you want a plot of Flux vs Time (y/n): ").lower()

water_depth = float(input("Ocean water depth in m: "))
mixed_layer_depth = float(input("Ocean mixing layer depth in m: "))
print(f"Radius of sun at {target_age}: {sun_radius(target_age)}")
start_time = time.time()
# endregion

# region Reading data files for Sun
file_path_fut_56 = r'Claire_Data\Sun5.6Gyr.txt'
data_fut_56 = pd.read_csv(file_path_fut_56, sep='\\s+', comment='#', names=['wavelength', 'flux'],
                          skiprows=103, nrows=5726)

file_path_00 = r'Claire_Data\SunModern.txt'
data_00 = pd.read_csv(file_path_00, sep='\\s+', comment='#', names=['wavelength', 'flux'],
                      skiprows=103, nrows=5726)

file_path_06 = r'Claire_Data\Sun0.6Ga.txt'
data_06 = pd.read_csv(file_path_06, sep='\\s+', comment='#', names=['wavelength', 'flux'],
                      skiprows=103, nrows=5726)

file_path_18 = r'Claire_Data\Sun1.8Ga.txt'
data_18 = pd.read_csv(file_path_18, sep='\\s+', comment='#', names=['wavelength', 'flux'],
                      skiprows=103, nrows=5726)

file_path_24 = r'Claire_Data\Sun2.4Ga.txt'
data_24 = pd.read_csv(file_path_24, sep='\\s+', comment='#', names=['wavelength', 'flux'],
                      skiprows=103, nrows=5726)

file_path_38 = r'Claire_Data\Sun3.8Ga.txt'
data_38 = pd.read_csv(file_path_38, sep='\\s+', comment='#', names=['wavelength', 'flux'],
                      skiprows=103, nrows=5726)

file_path_44 = r'Claire_Data\Sun4.4Ga.txt'
data_44 = pd.read_csv(file_path_44, sep='\\s+', comment='#', names=['wavelength', 'flux'],
                      skiprows=103, nrows=5726)
# endregion

# region Creating wavelengths and fluxes lists from data
wavelengths = np.array(data_00['wavelength'])
wavelengths_conv = 10 ** (-9) * np.array(data_06['wavelength'])
flux_fut_56 = h * c * 10 ** 4 * np.array(data_fut_56['flux']) / wavelengths_conv
flux_00 = h * c * 10 ** 4 * np.array(data_00['flux']) / wavelengths_conv
flux_06 = h * c * 10 ** 4 * np.array(data_06['flux']) / wavelengths_conv
flux_18 = h * c * 10 ** 4 * np.array(data_18['flux']) / wavelengths_conv
flux_24 = h * c * 10 ** 4 * np.array(data_24['flux']) / wavelengths_conv
flux_38 = h * c * 10 ** 4 * np.array(data_38['flux']) / wavelengths_conv
flux_44 = h * c * 10 ** 4 * np.array(data_44['flux']) / wavelengths_conv

flux_fut_56_photons = 10 ** 4 * np.array(data_fut_56['flux'])
flux_00_photons = 10 ** 4 * np.array(data_00['flux'])
flux_06_photons = 10 ** 4 * np.array(data_06['flux'])
flux_18_photons = 10 ** 4 * np.array(data_18['flux'])
flux_24_photons = 10 ** 4 * np.array(data_24['flux'])
flux_38_photons = 10 ** 4 * np.array(data_38['flux'])
flux_44_photons = 10 ** 4 * np.array(data_44['flux'])
# endregion

# region Defining regions for integrations
only_uva = [315, 400]
until_uva = [100, 400]
only_uvb = [280, 315]
until_uvb = [100, 315]
only_uvc = [100, 280]
until_uvc = [100, 280]
# endregion

# region Finding flux at Sun (before atmosphere)
ages_of_Sun = np.array([0.203, 0.803, 2.203, 2.803, 4.003, 4.603, 5.600])
flux_data = np.vstack([flux_44, flux_38, flux_24, flux_18, flux_06, flux_00, flux_fut_56])
flux_data_photons = np.vstack([flux_44_photons, flux_38_photons, flux_24_photons, flux_18_photons,
                               flux_06_photons, flux_00_photons, flux_fut_56_photons])
flux_needed = []
flux_needed_photons = []
smoothing_factor = 0

if target_age in ages_of_Sun:
    if target_age == 0.203:
        flux_needed = flux_44
        flux_needed_photons = flux_44_photons
    if target_age == 0.803:
        flux_needed = flux_38
        flux_needed_photons = flux_38_photons
    if target_age == 2.203:
        flux_needed = flux_24
        flux_needed_photons = flux_24_photons
    if target_age == 2.803:
        flux_needed = flux_18
        flux_needed_photons = flux_18_photons
    if target_age == 4.003:
        flux_needed = flux_06
        flux_needed_photons = flux_06_photons
    if target_age == 4.603:
        flux_needed = flux_00
        flux_needed_photons = flux_00_photons
    if target_age == 5.6:
        flux_needed = flux_fut_56
        flux_needed_photons = flux_fut_56_photons
else:
    flux_needed = smoothing_spline_interpolation_solar(target_age, ages_of_Sun, flux_data, smoothing_factor)
    flux_needed_photons = smoothing_spline_interpolation_solar(target_age, ages_of_Sun, flux_data_photons,
                                                               smoothing_factor)
# endregion

# region Atmosphere Modeling
mass_CO2 = 7.308032317 * 10 ** (-26)  # kg
mass_N2 = 4.652830391 * 10 ** (-26)
mass_O2 = 5.313724929 * 10 ** (-26)
mass_03 = 7.970587394 * 10 ** (-26)
g_Earth = 9.8  # N/kg

# region Reading Atmosphere files
file_path_CO2_acs = r'ACS_Data\CO2_ACS.csv'
data_CO2_acs = pd.read_csv(file_path_CO2_acs, comment='#', names=['wavelength', 'acs'], nrows=322)
data_CO2_acs_wl = data_CO2_acs['wavelength'].astype(float)
data_CO2_acs_acs = data_CO2_acs['acs'].astype(float)

file_path_N2_acs = r'ACS_Data\N2_ACS.csv'
data_N2_acs = pd.read_csv(file_path_N2_acs, comment='#', names=['wavelength', 'acs'], nrows=12)
data_N2_acs_wl = data_N2_acs['wavelength'].astype(float)
data_N2_acs_acs = data_N2_acs['acs'].astype(float)

file_path_O2_acs = r'ACS_Data\O2_ACS.csv'
data_O2_acs = pd.read_csv(file_path_O2_acs, comment='#', names=['wavelength', 'acs'], nrows=322)
data_O2_acs_wl = data_O2_acs['wavelength'].astype(float)
data_O2_acs_acs = data_O2_acs['acs'].astype(float)

file_path_O2_time = r'ACS_Data\O2_vs_Time.csv'
data_O2_time = pd.read_csv(file_path_O2_time, delimiter=",", comment='#', names=['time', 'pressure'], nrows=177)
data_O2_time_time = data_O2_time['time'].astype(float)
print(f"New O2 time data: {data_O2_time_time}")
data_O2_time_pressure = data_O2_time['pressure'].astype(float)
print(f"New O2 pressure data (in %PAL): {data_O2_time_pressure}")

file_path_O3_acs = r'ACS_Data\O3_ACS.csv'
data_O3_acs = pd.read_csv(file_path_O3_acs, comment='#', names=['wavelength', 'acs'], nrows=562)
data_O3_acs_wl = data_O3_acs['wavelength'].astype(float)
data_O3_acs_acs_special = data_O3_acs['acs'].astype(float)

file_path_O3_sigma = r'ACS_Data\O3_Sigma.csv'  # EDIT stuff here to implement ozone here
data_O3_sigma = pd.read_csv(file_path_O3_sigma, comment="#", names=['O2_amount', 'O3_sigma'])
data_O3_sigma_O2 = data_O3_sigma['O2_amount'].astype(float)
data_O3_sigma_sigma = data_O3_sigma['O3_sigma'].astype(float)
data_O3_sigma_sigma = data_O3_sigma_sigma * 2.687 * 10 ** 20 * mass_03

# # Amount of ozone: O3_3 < O3_1 < O3_2
# file_path_O3_density = r'ACS_Data\O3_Density_3.csv'  # Change to needed density
# # data_O3_density = pd.read_csv(file_path_O3_density, comment='#', names=['altitude', 'density'], nrows=147)  # For O3_Density_1
# # data_O3_density = pd.read_csv(file_path_O3_density, comment='#', names=['altitude', 'density'], nrows=163)  # For O3_Density_2
# data_O3_density = pd.read_csv(file_path_O3_density, comment='#', names=['altitude', 'density'], nrows=263)  # For O3_Density_3
# data_O3_density_altitude = data_O3_density['altitude'].astype(float)
# data_O3_density_altitude = data_O3_density_altitude * 1000  # For O3_Density_2 and 3
# data_O3_density_density = data_O3_density['density'].astype(float)
# # data_O3_density_density = data_O3_density_density * 10**(-3)  # For O3_Density_1
# data_O3_density_density = data_O3_density_density * 10 ** 6 * 0.0480 / (6.02 * 10 ** 23)  # For O3_Density_2 and 3
# O3_sigma = simpson(y=data_O3_density_density, x=data_O3_density_altitude)
# # O3_sigma = 0.005  # Trial and error sigma
# print(f"Sigma: {O3_sigma}")
# # endregion

# region Find ACS lists
CO2_interpolator = interp1d(data_CO2_acs_wl, data_CO2_acs_acs, kind='linear', bounds_error=False,
                            fill_value="extrapolate")
CO2_acs_lst = CO2_interpolator(wavelengths)
CO2_threshold_wavelength = 317.1748803
# noinspection PyUnresolvedReferences
CO2_acs_lst[wavelengths > CO2_threshold_wavelength] = 0

N2_interpolator = interp1d(data_N2_acs_wl, data_N2_acs_acs, kind='linear', bounds_error=False,
                           fill_value="extrapolate")
N2_acs_lst = N2_interpolator(wavelengths)
N2_threshold_wavelength = 122.250272
# noinspection PyUnresolvedReferences
N2_acs_lst[wavelengths > N2_threshold_wavelength] = 0

O2_interpolator = interp1d(data_O2_acs_wl, data_O2_acs_acs, kind='linear', bounds_error=False,
                           fill_value="extrapolate")
O2_acs_lst = O2_interpolator(wavelengths)
O2_threshold_wavelength = 314.6879444
# noinspection PyUnresolvedReferences
O2_acs_lst[wavelengths > O2_threshold_wavelength] = 0
# print(O2_acs_lst)

O3_interpolator = interp1d(data_O3_acs_wl, data_O3_acs_acs_special, kind='linear', bounds_error=False,
                           fill_value="extrapolate")
O3_acs_lst_special = O3_interpolator(wavelengths)
O3_threshold_wavelength = 448.5001484
# noinspection PyUnresolvedReferences
O3_acs_lst_special[wavelengths > O3_threshold_wavelength] = 0

O3_sigma_interpolator = interp1d(data_O3_sigma_O2, data_O3_sigma_sigma, kind='linear', bounds_error=False,
                                 # TODO: See if it's ok to have linear
                                 fill_value="extrapolate")
# endregion

# region Perform spline interpolation on the log-transformed absorption cross-section data -> Getting tau list
def tau_value(t_age):
    if t_age < 2.203:  # Articles about CO2 and pressure
        pressure_past = 23000  # Pa
        percentage_CO2_past = 0.70
        pressure_CO2_past = pressure_past * percentage_CO2_past  # Pa
        percentage_N2_past = 0.30
        pressure_N2_past = pressure_past * percentage_N2_past  # Pa
        CO2_tau_lst = CO2_acs_lst * 10 ** (-4) * pressure_CO2_past / (mass_CO2 * g_Earth)
        N2_tau_lst = N2_acs_lst * 10 ** (-4) * pressure_N2_past / (mass_N2 * g_Earth)
        total_tau_lst = N2_tau_lst + CO2_tau_lst
        # for w in range(len(total_tau_lst)):
        #     if total_tau_lst[w] < 0:
        #         print(f"Negative tau value at {t_age}: {total_tau_lst[w]}")
    # elif 2.203 <= t_age < 4.003:
    #     O2_time_interpolator = interp1d(data_O2_time_time, data_O2_time_pressure, kind='linear', bounds_error=False,
    #                                     fill_value='extrapolate')
    #     pressure_O2_past = O2_time_interpolator(4.603 - t_age)
    #     # print(pressure_O2_past)
    #     pressure_interpolation = interp1d([2.203, 4.603], [23000, 101300], kind='linear', bounds_error=False,
    #                                       fill_value='extrapolate')
    #     pressure_past = pressure_interpolation(t_age)
    #     percentage_N2_interpolation = interp1d([2.203, 4.603], [0.30, 0.7808], kind='linear',
    #                                            bounds_error=False, fill_value='extrapolate')
    #     percentage_N2_past = percentage_N2_interpolation(t_age)
    #     pressure_N2_past = pressure_past * percentage_N2_past  # Pa
    #     # percentage_CO2_interpolation = interp1d([2.166031101, 4.603], [0.30, 0.7808], kind='linear',
    #     #                                        bounds_error=False, fill_value='extrapolate')
    #     # percentage_CO2_past = percentage_CO2_interpolation(t_age)
    #     # pressure_CO2_past = pressure_past * percentage_CO2_past
    #     percentage_CO2_interpolation = interp1d([2.203, 4.603], [0.70, 0.00042], kind='linear',
    #                                             bounds_error=False, fill_value='extrapolate')
    #     percentage_CO2_past = percentage_CO2_interpolation(t_age)
    #     pressure_CO2_past = pressure_past * percentage_CO2_past
    #     # pressure_CO2_past = pressure_past - pressure_N2_past - pressure_O2_past - 0.001
    #     # Subtracted -0.001 as an estimate of total partial pressure of other gases
    #     CO2_tau_lst = CO2_acs_lst * 10 ** (-4) * pressure_CO2_past / (mass_CO2 * g_Earth)
    #     N2_tau_lst = N2_acs_lst * 10 ** (-4) * pressure_N2_past / (mass_N2 * g_Earth)
    #     O2_tau_lst = O2_acs_lst * 10 ** (-4) * pressure_O2_past / (mass_O2 * g_Earth)
    #     total_tau_lst = N2_tau_lst + CO2_tau_lst + O2_tau_lst
    # for w in range(len(total_tau_lst)):
    #     if total_tau_lst[w] < 0:
    #         print(f"Negative tau value at {t_age}: {total_tau_lst[w]}")
    else:  # Articles / general data about pressure and gases
        # print(pressure_O2_past)
        # pressure_interpolation = interp1d([2.166031101, 4.603], [23000, 101300], kind='linear',bounds_error=False,
        #                                 fill_value='extrapolate')
        # pressure_past = pressure_interpolation(t_age)
        if t_age >= 4.603:
            t_age = 4.603
        else:
            t_age = t_age
        O2_time_interpolator = interp1d(data_O2_time_time, data_O2_time_pressure,
                                        kind='linear')  # TODO: Define perhaps as a function
        pressure_O2_past = O2_time_interpolator(4.603 - t_age) * 21222.35
        print(f"Past pressure at {4.603 - t_age} Ga: {pressure_O2_past}")
        pressure_interpolation = interp1d([2.203, 4.603], [23000, 101300], kind='linear', bounds_error=False,
                                          fill_value='extrapolate')
        pressure_past = pressure_interpolation(t_age)
        percentage_N2_interpolation = interp1d([2.203, 4.603], [0.30, 0.7808], kind='linear',
                                               bounds_error=False, fill_value='extrapolate')
        percentage_N2_past = percentage_N2_interpolation(t_age)
        pressure_N2_past = pressure_past * percentage_N2_past
        percentage_CO2_interpolation = interp1d([2.203, 4.603], [0.70, 0.00042], kind='linear',
                                                bounds_error=False, fill_value='extrapolate')
        percentage_CO2_past = percentage_CO2_interpolation(t_age)
        pressure_CO2_past = pressure_past * percentage_CO2_past
        CO2_tau_lst = CO2_acs_lst * 10 ** (-4) * pressure_CO2_past / (mass_CO2 * g_Earth)
        N2_tau_lst = N2_acs_lst * 10 ** (-4) * pressure_N2_past / (mass_N2 * g_Earth)
        O2_tau_lst = O2_acs_lst * 10 ** (-4) * pressure_O2_past / (mass_O2 * g_Earth)
        O3_tau_lst = O3_acs_lst_special * 10 ** (-4) * O3_sigma_interpolator(
            O2_time_interpolator(4.603 - t_age)) / mass_03
        total_tau_lst = N2_tau_lst + CO2_tau_lst + O2_tau_lst + O3_tau_lst
        # print(f"Total tau list: {total_tau_lst}")
        # for w in range(len(O2_acs_lst)):
        #     if O2_acs_lst[w] < 0:
        #         print(f"Negative tau value at O2 ACS list at {t_age}: {O2_acs_lst[w]}")
        # if pressure_O2_past < 0:
        #     print(f"Negative tau value at O2 Pressure list at {t_age}: {pressure_O2_past}")
        # print(O3_tau_lst)
        # print(total_tau_lst)
    return total_tau_lst


def N2_pressure(t_age):
    if t_age < 2.203:  # Articles about CO2 and pressure
        pressure_past = 0.230  # bar
        percentage_N2_past = 0.30
        pressure_N2_past = pressure_past * percentage_N2_past
        return pressure_N2_past
    elif 2.203 <= t_age < 4.003:
        pressure_interpolation = interp1d([2.203, 4.603], [0.230, 1.013], kind='linear', bounds_error=False,
                                          fill_value='extrapolate')
        pressure_past = pressure_interpolation(t_age)
        percentage_N2_interpolation = interp1d([2.203, 4.603], [0.30, 0.7808], kind='linear',
                                               bounds_error=False, fill_value='extrapolate')
        percentage_N2_past = percentage_N2_interpolation(t_age)
        pressure_N2_past = pressure_past * percentage_N2_past  # Pa
        return pressure_N2_past
    else:
        if t_age >= 4.603:
            t_age = 4.603
        else:
            t_age = t_age
        pressure_interpolation = interp1d([2.203, 4.603], [0.230, 1.013], kind='linear', bounds_error=False,
                                          fill_value='extrapolate')
        pressure_past = pressure_interpolation(t_age)
        percentage_N2_interpolation = interp1d([2.203, 4.603], [0.30, 0.7808], kind='linear',
                                               bounds_error=False, fill_value='extrapolate')
        percentage_N2_past = percentage_N2_interpolation(t_age)
        pressure_N2_past = pressure_past * percentage_N2_past
        return pressure_N2_past


def CO2_pressure(t_age):
    if t_age < 2.203:  # Articles about CO2 and pressure
        pressure_past = 0.230  # bar
        percentage_CO2_past = 0.70
        pressure_CO2_past = pressure_past * percentage_CO2_past  # Pa
        return pressure_CO2_past
    elif 2.203 <= t_age < 4.003:
        pressure_interpolation = interp1d([2.203, 4.603], [0.230, 1.013], kind='linear', bounds_error=False,
                                          fill_value='extrapolate')
        pressure_past = pressure_interpolation(t_age)
        percentage_CO2_interpolation = interp1d([2.203, 4.603], [0.70, 0.00042], kind='linear',
                                                bounds_error=False, fill_value='extrapolate')
        percentage_CO2_past = percentage_CO2_interpolation(t_age)
        pressure_CO2_past = pressure_past * percentage_CO2_past
        return pressure_CO2_past
    else:
        if t_age >= 4.603:
            t_age = 4.603
        else:
            t_age = t_age
        pressure_interpolation = interp1d([2.203, 4.603], [0.230, 1.013], kind='linear', bounds_error=False,
                                          fill_value='extrapolate')
        pressure_past = pressure_interpolation(t_age)
        percentage_CO2_interpolation = interp1d([2.203, 4.603], [0.70, 0.00042], kind='linear',
                                                bounds_error=False, fill_value='extrapolate')
        percentage_CO2_past = percentage_CO2_interpolation(t_age)
        pressure_CO2_past = pressure_past * percentage_CO2_past
        return pressure_CO2_past


def O2_pressure(t_age):
    if t_age < 2.203:  # Articles about CO2 and pressure
        pressure_O2_past = 0
        return pressure_O2_past
    else:
        if t_age >= 4.603:
            t_age = 4.603
        else:
            t_age = t_age
        O2_time_interpolator = interp1d(data_O2_time_time, data_O2_time_pressure, kind='linear', bounds_error=False,
                                        fill_value='extrapolate')
        pressure_O2_past = O2_time_interpolator(4.603 - t_age) * 21222.35 / 100000  # bar
        return pressure_O2_past


def O3_presence(t_age):  # TODO: fix graphing with new O3
    if t_age < 2.203:
        return 0
    else:
        if t_age >= 4.603:
            t_age = 4.603
        O2_time_interpolator = interp1d(data_O2_time_time, data_O2_time_pressure, kind='linear', bounds_error=False,
                                        fill_value='extrapolate')
        print(f"O2 compared from {t_age} to today: {O2_time_interpolator(4.603 - t_age)}")
        print(f"O3 layer: {O3_sigma_interpolator(O2_time_interpolator(4.603 - t_age))}")
        return O3_sigma_interpolator(O2_time_interpolator(4.603 - t_age))


def total_pressure(t_age):
    if t_age < 2.203:  # Articles about CO2 and pressure
        pressure_past = 0.23  # Pa
        return pressure_past
    elif 2.203 <= t_age < 4.003:
        pressure_interpolation = interp1d([2.203, 4.603], [0.230, 1.013], kind='linear', bounds_error=False,
                                          fill_value='extrapolate')
        pressure_past = pressure_interpolation(t_age)
        return pressure_past
    else:
        if t_age >= 4.603:
            t_age = 4.603
        else:
            t_age = t_age
        pressure_interpolation = interp1d([2.203, 4.603], [0.23000, 1.013], kind='linear', bounds_error=False,
                                          fill_value='extrapolate')
        pressure_past = pressure_interpolation(t_age)
        return pressure_past
# endregion

# Calculate flux post atmosphere
flux_post_atmosphere = solar_to_atmosphere_flux(flux_needed, tau_value(target_age))
flux_post_atmosphere_photons = solar_to_atmosphere_flux(flux_needed_photons, tau_value(target_age))
# endregion


# region Ocean Water Modeling
file_path_water_1 = r'Water_Data\Water_Data_1.csv'
file_path_water_2 = r'Water_Data\Water_Data_2.txt'
file_path_water_3 = r'Water_Data\Water_Data_3.txt'  # Change to 1, 2, 3, 4 as needed (csv for 1 and 4, txt for 2 and 3)

data_water_1 = pd.read_csv(file_path_water_1, comment='#', names=['wavelength', 'coefficient'],
                           nrows=157)  # For Water_Data_1
data_water_2 = pd.read_csv(file_path_water_2, sep='\\s+', comment='#', names=['wavelength', 'coefficient'],
                           skiprows=4, nrows=21)  # For Water_Data_2
data_water_3 = pd.read_csv(file_path_water_3, sep='\\s+', comment='#', names=['wavelength', 'coefficient'],
                           skiprows=4, nrows=9)  # For Water_Data_3

water_wavelengths_1 = data_water_1['wavelength'].astype(float)
water_wavelengths_2 = data_water_2['wavelength'].astype(float)
water_wavelengths_3 = data_water_3['wavelength'].astype(float)

water_coefficients_1 = data_water_1['coefficient'].astype(float)
water_coefficients_2 = data_water_2['coefficient'].astype(float)
water_coefficients_3 = data_water_3['coefficient'].astype(float)

water_coefficients_2 = water_coefficients_2 * 100  # For 2, 3, 4
water_coefficients_3 = water_coefficients_3 * 100

# water_interpolator_1 = UnivariateSpline(water_wavelengths_1, water_coefficients_1)
water_interpolator_1 = interp1d(water_wavelengths_1, water_coefficients_1, kind='linear', bounds_error=False,
                                fill_value="extrapolate")
water_coefficients_lst_1 = np.array(water_interpolator_1(wavelengths))
# water_interpolator_2 = UnivariateSpline(water_wavelengths_2, water_coefficients_2)
water_interpolator_2 = interp1d(water_wavelengths_2, water_coefficients_2, kind='linear', bounds_error=False,
                                fill_value="extrapolate")
water_coefficients_lst_2 = np.array(water_interpolator_2(wavelengths))
# water_interpolator_3 = UnivariateSpline(water_wavelengths_3, water_coefficients_3)
water_interpolator_3 = interp1d(water_wavelengths_3, water_coefficients_3, kind='linear', bounds_error=False,
                                fill_value="extrapolate")
water_coefficients_lst_3 = np.array(water_interpolator_3(wavelengths))

main_water_coefficient_lst = water_coefficients_lst_1
flux_post_water_specific = atmosphere_to_water_flux(flux_post_atmosphere, main_water_coefficient_lst, water_depth)
flux_post_water_photons = atmosphere_to_water_flux(flux_post_atmosphere_photons, main_water_coefficient_lst,
                                                   water_depth)
flux_post_water_mixed = atmosphere_to_mixed_water_flux(flux_post_atmosphere, main_water_coefficient_lst, mixed_layer_depth)
# endregion

# region Defining Masks
wavelength_mask = (wavelengths >= wavelength_integration_start) & (wavelengths <= wavelength_integration_end)
wavelength_interval = wavelengths[wavelength_mask]
water_coefficient_lst_need = main_water_coefficient_lst[wavelength_mask]

organic_wavelength_mask = (wavelengths >= 100) & (wavelengths <= single_CN_bond_wavelength)
organic_wavelength_interval = wavelengths[organic_wavelength_mask]
water_coefficient_lst_org = main_water_coefficient_lst[organic_wavelength_mask]
# endregion

# region Average Flux in Mixed Layer
depth_steps = 0.1
depth_list = [depth_steps]
depth_change = 0.1
max_depth = 100

integrated_flux_water_specific = simpson(y=atmosphere_to_water_flux(flux_post_atmosphere, main_water_coefficient_lst,
                                                                depth_steps)[wavelength_mask], x=wavelength_interval)
integrated_flux_lst_water_specific = [integrated_flux_water_specific]

integrated_flux_water_mixed = simpson(y=atmosphere_to_mixed_water_flux(flux_post_atmosphere, main_water_coefficient_lst,
                                                                   depth_steps)[wavelength_mask], x=wavelength_interval)
integrated_flux_lst_water_mixed = [integrated_flux_water_mixed]

while depth_steps < max_depth:
    depth_steps += depth_change
    depth_list.append(depth_steps)

    integrated_flux_water_specific = simpson(y=atmosphere_to_water_flux(flux_post_atmosphere,
                                                                    main_water_coefficient_lst, depth_steps)[wavelength_mask],
                                           x=wavelength_interval)
    integrated_flux_lst_water_specific.append(integrated_flux_water_specific)

    integrated_flux_water_mixed = simpson(y=atmosphere_to_mixed_water_flux(flux_post_atmosphere, main_water_coefficient_lst,
                                                                       depth_steps)[wavelength_mask], x=wavelength_interval)
    integrated_flux_lst_water_mixed.append(integrated_flux_water_mixed)
# endregion

# region Making time and flux lists for plotting
time_steps = 0.6  # in Gyr
time_list = [time_steps]
time_change = 0.1  # in Gyr
max_time = 5.6  # in Gyr
flux_bol = flux_bolometric_calc(time_steps)
flux_bol_list = [flux_bol]

flux_interval = smoothing_spline_interpolation_solar(time_steps, ages_of_Sun, flux_data,
                                                     smoothing_factor)[wavelength_mask]
integrated_flux_solar = simpson(y=flux_interval, x=wavelength_interval)
integrated_flux_lst_solar = [integrated_flux_solar]
integrated_flux_atmosphere = simpson(y=solar_to_atmosphere_flux(flux_interval, tau_value(time_steps)[wavelength_mask]),
                                   x=wavelength_interval)
integrated_flux_lst_atmosphere = [integrated_flux_atmosphere]
integrated_flux_water = simpson(y=atmosphere_to_water_flux(solar_to_atmosphere_flux(flux_interval,
                                                                                tau_value(time_steps)[wavelength_mask]),
                                                       water_coefficient_lst_need, water_depth), x=wavelength_interval)
integrated_flux_lst_water = [integrated_flux_water]

integrated_flux_water_mixed_time = simpson(y=atmosphere_to_mixed_water_flux(solar_to_atmosphere_flux(flux_interval,
                                                                                    tau_value(time_steps)[
                                                                                        wavelength_mask]),
                                                           water_coefficient_lst_need, mixed_layer_depth),
                                  x=wavelength_interval)
integrated_flux_lst_water_mixed_time = [integrated_flux_water_mixed_time]

flux_interval_photons = smoothing_spline_interpolation_solar(time_steps, ages_of_Sun, flux_data_photons,
                                                             smoothing_factor)[wavelength_mask]
integrated_flux_solar_photons = simpson(y=flux_interval_photons, x=wavelength_interval)
integrated_flux_lst_solar_photons = [integrated_flux_solar_photons]
integrated_flux_atmosphere_photons = simpson(y=solar_to_atmosphere_flux(flux_interval_photons,
                                                                    tau_value(time_steps)[wavelength_mask]),
                                           x=wavelength_interval)
integrated_flux_lst_atmosphere_photons = [integrated_flux_atmosphere_photons]
integrated_flux_water_photons = simpson(y=atmosphere_to_water_flux_photons(solar_to_atmosphere_flux_photons(
    flux_interval_photons, tau_value(time_steps)[wavelength_mask]), water_coefficient_lst_need, water_depth),
    x=wavelength_interval)
integrated_flux_lst_water_photons = [integrated_flux_water_photons]

organic_flux_interval_photons = smoothing_spline_interpolation_solar(time_steps, ages_of_Sun, flux_data_photons,
                                                                     smoothing_factor)[organic_wavelength_mask]
organic_integrated_flux_atmosphere_photons = simpson(y=solar_to_atmosphere_flux_photons(organic_flux_interval_photons,
                                                                                    tau_value(time_steps)
                                                                                    [organic_wavelength_mask]),
                                                   x=organic_wavelength_interval)
organic_integrated_flux_atmosphere_grams = 110 * 300 * 1.66054 * 10 ** (-24) * \
                                           organic_integrated_flux_atmosphere_photons / 300  # TODO fix numbers
organic_integrated_flux_lst_atmosphere_grams = [organic_integrated_flux_atmosphere_grams]

organic_integrated_flux_water_photons = simpson(y=atmosphere_to_water_flux_photons(solar_to_atmosphere_flux_photons(
    organic_flux_interval_photons, tau_value(time_steps)[organic_wavelength_mask]), water_coefficient_lst_org,
    water_depth),
    x=organic_wavelength_interval)
organic_integrated_flux_water_grams = 110 * 300 * 1.66054 * 10 ** (-24) * \
                                      organic_integrated_flux_water_photons / 300  # TODO fix numbers
organic_integrated_flux_lst_water_grams = [organic_integrated_flux_water_grams]

organic_integrated_flux_mixed_water_photons = simpson(y=atmosphere_to_mixed_water_flux_photons(solar_to_atmosphere_flux_photons(
    organic_flux_interval_photons, tau_value(time_steps)[organic_wavelength_mask]), water_coefficient_lst_org,
    mixed_layer_depth),
    x=organic_wavelength_interval)
organic_integrated_flux_mixed_water_grams = 110 * 300 * 1.66054 * 10 ** (-24) * \
                                      organic_integrated_flux_mixed_water_photons / 300  # TODO fix numbers
organic_integrated_flux_lst_mixed_water_grams = [organic_integrated_flux_mixed_water_grams]

N2_pressure_lst = [N2_pressure(time_steps)]
CO2_pressure_lst = [CO2_pressure(time_steps)]
O2_pressure_lst = [O2_pressure(time_steps)]
O3_presence_lst = [O3_presence(time_steps)]
total_pressure_lst = [total_pressure(time_steps)]

# region Specific UV code
# uva_mask = (wavelengths >= only_uva[0]) & (wavelengths <= only_uva[1])
# uva_wavelengths = wavelengths[uva_mask]
# uva_flux = smoothing_spline_interpolation(time_steps, ages_of_Sun, flux_data, smoothing_factor)[uva_mask]
# uvb_mask = (wavelengths >= only_uvb[0]) & (wavelengths <= only_uvb[1])
# uvb_wavelengths = wavelengths[uvb_mask]
# uvb_flux = smoothing_spline_interpolation(time_steps, ages_of_Sun, flux_data, smoothing_factor)[uvb_mask]
# uvc_mask = (wavelengths >= only_uvc[0]) & (wavelengths <= only_uvc[1])
# uvc_wavelengths = wavelengths[uvc_mask]
# uvc_flux = smoothing_spline_interpolation(time_steps, ages_of_Sun, flux_data, smoothing_factor)[uvc_mask]
#
# integrated_uva_flux = simpson(y=uva_flux, x=uva_wavelengths)
# integrated_uva_flux_lst = [integrated_uva_flux]
#
# integrated_uvb_flux = simpson(y=uvb_flux, x=uvb_wavelengths)
# integrated_uvb_flux_lst = [integrated_uvb_flux]
#
# integrated_uvc_flux = simpson(y=uvc_flux, x=uvc_wavelengths)
# integrated_uvc_flux_lst = [integrated_uvc_flux]
# endregion

while time_steps < max_time:
    time_steps += time_change
    time_list.append(time_steps)

    flux_bol = flux_bolometric_calc(time_steps)
    flux_bol_list.append(flux_bol)

    flux_interval = smoothing_spline_interpolation_solar(time_steps, ages_of_Sun, flux_data,
                                                         smoothing_factor)[wavelength_mask]
    integrated_flux_solar = simpson(y=flux_interval, x=wavelength_interval)
    integrated_flux_lst_solar.append(integrated_flux_solar)

    integrated_flux_atmosphere = simpson(y=solar_to_atmosphere_flux(flux_interval, tau_value(time_steps)[wavelength_mask]),
                                       x=wavelength_interval)
    integrated_flux_lst_atmosphere.append(integrated_flux_atmosphere)
    integrated_flux_water = simpson(y=atmosphere_to_water_flux(solar_to_atmosphere_flux(flux_interval,
                                                                                    tau_value(time_steps)[
                                                                                        wavelength_mask]),
                                                           water_coefficient_lst_need, water_depth),
                                  x=wavelength_interval)
    integrated_flux_lst_water.append(integrated_flux_water)

    integrated_flux_water_mixed_time = simpson(y=atmosphere_to_mixed_water_flux(solar_to_atmosphere_flux(flux_interval,
                                                                                    tau_value(time_steps)[
                                                                                        wavelength_mask]),
                                                           water_coefficient_lst_need, mixed_layer_depth),
                                  x=wavelength_interval)
    integrated_flux_lst_water_mixed_time.append(integrated_flux_water_mixed_time)

    flux_interval_photons = smoothing_spline_interpolation_solar(time_steps, ages_of_Sun, flux_data_photons,
                                                                 smoothing_factor)[wavelength_mask]
    integrated_flux_solar_photons = simpson(y=flux_interval_photons, x=wavelength_interval)
    integrated_flux_lst_solar_photons.append(integrated_flux_solar_photons)

    integrated_flux_atmosphere_photons = simpson(y=solar_to_atmosphere_flux(flux_interval_photons,
                                                                        tau_value(time_steps)[wavelength_mask]),
                                               x=wavelength_interval)
    integrated_flux_lst_atmosphere_photons.append(integrated_flux_atmosphere_photons)
    integrated_flux_water_photons = simpson(y=atmosphere_to_water_flux_photons(solar_to_atmosphere_flux_photons(
        flux_interval_photons, tau_value(time_steps)[wavelength_mask]), water_coefficient_lst_need, water_depth),
        x=wavelength_interval)
    integrated_flux_lst_water_photons.append(integrated_flux_water_photons)

    organic_flux_interval_photons = smoothing_spline_interpolation_solar(time_steps, ages_of_Sun, flux_data_photons,
                                                                         smoothing_factor)[organic_wavelength_mask]
    organic_integrated_flux_atmosphere_photons = simpson(y=solar_to_atmosphere_flux(organic_flux_interval_photons,
                                                                                tau_value(time_steps)
                                                                                [organic_wavelength_mask]),
                                                       x=organic_wavelength_interval)
    organic_integrated_flux_atmosphere_grams = 110 * 300 * 1.66054 * 10 ** (-24) * \
                                               organic_integrated_flux_atmosphere_photons / 300
    organic_integrated_flux_lst_atmosphere_grams.append(organic_integrated_flux_atmosphere_grams)

    organic_integrated_flux_water_photons = simpson(y=atmosphere_to_water_flux_photons(solar_to_atmosphere_flux_photons(
        organic_flux_interval_photons, tau_value(time_steps)[organic_wavelength_mask]), water_coefficient_lst_org,
        water_depth),
        x=organic_wavelength_interval)
    organic_integrated_flux_water_grams = 110 * 300 * 1.66054 * 10 ** (-24) * \
                                          organic_integrated_flux_water_photons / 300
    organic_integrated_flux_lst_water_grams.append(organic_integrated_flux_water_grams)

    organic_integrated_flux_mixed_water_photons = simpson(
        y=atmosphere_to_mixed_water_flux_photons(solar_to_atmosphere_flux_photons(
            organic_flux_interval_photons, tau_value(time_steps)[organic_wavelength_mask]), water_coefficient_lst_org,
            mixed_layer_depth),
        x=organic_wavelength_interval)
    organic_integrated_flux_mixed_water_grams = 110 * 300 * 1.66054 * 10 ** (-24) * \
                                                organic_integrated_flux_mixed_water_photons / 300  # TODO fix numbers
    organic_integrated_flux_lst_mixed_water_grams.append(organic_integrated_flux_mixed_water_grams)

    N2_pressure_lst.append(N2_pressure(time_steps))
    CO2_pressure_lst.append(CO2_pressure(time_steps))
    O2_pressure_lst.append(O2_pressure(time_steps))
    O3_presence_lst.append(O3_presence(time_steps))
    total_pressure_lst.append(total_pressure(time_steps))

    # region Specific UV code
    # uva_flux = smoothing_spline_interpolation(time_steps, ages_of_Sun, flux_data, smoothing_factor)[uva_mask]
    # integrated_uva_flux = simpson(y=uva_flux, x=uva_wavelengths)
    # integrated_uva_flux_lst.append(integrated_uva_flux)
    #
    # uvb_flux = smoothing_spline_interpolation(time_steps, ages_of_Sun, flux_data, smoothing_factor)[uvb_mask]
    # integrated_uvb_flux = simpson(y=uvb_flux, x=uvb_wavelengths)
    # integrated_uvb_flux_lst.append(integrated_uvb_flux)
    #
    # uvc_flux = smoothing_spline_interpolation(time_steps, ages_of_Sun, flux_data, smoothing_factor)[uvc_mask]
    # integrated_uvc_flux = simpson(y=uvc_flux, x=uvc_wavelengths)
    # integrated_uvc_flux_lst.append(integrated_uvc_flux)
    # endregion
# endregion

# region Creating points to plot on Flux vs Time
# region Specific UV code
# point_x1a = 4.603
# point_y1a = simpson(y=smoothing_spline_interpolation(4.603, ages_of_Sun, flux_data, smoothing_factor)[uva_mask],
#                   x=uva_wavelengths)
# point_x2a = target_age
# point_y2a = simpson(y=smoothing_spline_interpolation(target_age, ages_of_Sun, flux_data, smoothing_factor)[uva_mask],
#                   x=uva_wavelengths)
# point_x1b = 4.603
# point_y1b = simpson(y=smoothing_spline_interpolation(4.603, ages_of_Sun, flux_data, smoothing_factor)[uvb_mask],
#                   x=uvb_wavelengths)
# point_x2b = target_age
# point_y2b = simpson(y=smoothing_spline_interpolation(target_age, ages_of_Sun, flux_data, smoothing_factor)[uvb_mask],
#                   x=uvb_wavelengths)
# point_x1c = 4.603
# point_y1c = simpson(y=smoothing_spline_interpolation(4.603, ages_of_Sun, flux_data, smoothing_factor)[uvc_mask],
#                   x=uvc_wavelengths)
# point_x2c = target_age
# point_y2c = simpson(y=smoothing_spline_interpolation(target_age, ages_of_Sun, flux_data, smoothing_factor)[uvc_mask],
#                   x=uvc_wavelengths)
# endregion

# point_x1 = 4.603
# point_y1 = simpson(y=smoothing_spline_interpolation_solar(4.603, ages_of_Sun, flux_data, smoothing_factor)[wavelength_mask],
#                  x=wavelength_interval)
# point_z1 = simpson(y=smoothing_spline_interpolation_solar(4.603, ages_of_Sun, flux_data_photons,
#                                                       smoothing_factor)[wavelength_mask], x=wavelength_interval)
#
# point_x2 = target_age
# point_y2 = simpson(y=smoothing_spline_interpolation_solar(target_age, ages_of_Sun, flux_data,
#                                                       smoothing_factor)[wavelength_mask], x=wavelength_interval)
# point_z2 = simpson(y=smoothing_spline_interpolation_solar(target_age, ages_of_Sun, flux_data_photons,
#                                                       smoothing_factor)[wavelength_mask], x=wavelength_interval)
#
# point_x1a = 4.603
# point_y1a = simpson(y=solar_to_atmosphere_flux(flux_00[wavelength_mask], tau_value(4.603)[wavelength_mask]),
#                   x=wavelength_interval)
# point_z1a = simpson(y=solar_to_atmosphere_flux(flux_00_photons[wavelength_mask], tau_value(4.603)[wavelength_mask]),
#                   x=wavelength_interval)
# point_N2_1 = N2_pressure(4.603)
# point_CO2_1 = CO2_pressure(4.603)
# point_O2_1 = O2_pressure(4.603)
# point_O3_1 = O3_presence(4.603)
# point_atm_1 = total_pressure(4.603)
# point_w1 = simpson(y=atmosphere_to_water_flux(solar_to_atmosphere_flux(flux_00[wavelength_mask],
#                                                                    tau_value(4.603)[wavelength_mask]),
#                                           water_coefficient_lst_need, water_depth), x=wavelength_interval)
# point_wp1 = simpson(y=atmosphere_to_water_flux(solar_to_atmosphere_flux_photons(flux_00_photons[wavelength_mask],
#                                                                             tau_value(4.603)[wavelength_mask]),
#                                            water_coefficient_lst_need, water_depth), x=wavelength_interval)
# point_alpha1 = 110 * 300 * 1.66054 * 10 ** (-24) * simpson(y=solar_to_atmosphere_flux_photons(
#     flux_00_photons[organic_wavelength_mask], tau_value(4.603)[organic_wavelength_mask]),
#     x=organic_wavelength_interval) / 300
# point_beta1 = 110 * 300 * 1.66054 * 10 ** (-24) * simpson(
#     y=atmosphere_to_water_flux_photons(solar_to_atmosphere_flux_photons(
#         flux_00_photons[organic_wavelength_mask], tau_value(4.603)[organic_wavelength_mask]), water_coefficient_lst_org,
#         water_depth),
#     x=organic_wavelength_interval) / 300
#
# point_x2a = target_age
# point_y2a = simpson(y=solar_to_atmosphere_flux(flux_needed[wavelength_mask], tau_value(target_age)[wavelength_mask]),
#                   x=wavelength_interval)
# point_z2a = simpson(y=solar_to_atmosphere_flux(flux_needed_photons[wavelength_mask],
#                                            tau_value(target_age)[wavelength_mask]), x=wavelength_interval)
# point_N2_2 = N2_pressure(target_age)
# point_CO2_2 = CO2_pressure(target_age)
# point_O2_2 = O2_pressure(target_age)
# point_O3_2 = O3_presence(target_age)
# point_atm_2 = total_pressure(target_age)
# point_w2 = simpson(y=atmosphere_to_water_flux(solar_to_atmosphere_flux(flux_needed[wavelength_mask],
#                                                                    tau_value(target_age)[wavelength_mask]),
#                                           water_coefficient_lst_need, water_depth), x=wavelength_interval)
# point_wp2 = simpson(y=atmosphere_to_water_flux(solar_to_atmosphere_flux_photons(flux_needed_photons[wavelength_mask],
#                                                                             tau_value(target_age)[wavelength_mask]),
#                                            water_coefficient_lst_need, water_depth), x=wavelength_interval)
# point_alpha2 = 110 * 300 * 1.66054 * 10 ** (-24) * simpson(y=solar_to_atmosphere_flux_photons(
#     flux_needed_photons[organic_wavelength_mask], tau_value(target_age)[organic_wavelength_mask]),
#     x=organic_wavelength_interval) / 300
# point_beta2 = 110 * 300 * 1.66054 * 10 ** (-24) * simpson(
#     y=atmosphere_to_water_flux_photons(solar_to_atmosphere_flux_photons(
#         flux_needed_photons[organic_wavelength_mask], tau_value(target_age)[organic_wavelength_mask]),
#         water_coefficient_lst_org,
#         water_depth),
#     x=organic_wavelength_interval) / 300

# endregion

directory = r'Results'

# Create the directory if it doesn't exist
if not os.path.exists(directory):
    os.makedirs(directory)

# region Plotting
plt.figure()
plt.plot(depth_list, integrated_flux_lst_water_specific, linewidth=1)
plt.title("Flux Under Water (at Specific Depth) vs. Specific Depth")
plt.xlabel("Depth, m")
plt.ylabel("Flux, W/m^2")
plt.tight_layout()
timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
file_path_save = os.path.join(directory, f"Specific_Flux_vs_Specific_Depth_{timestamp}.pdf")
plt.savefig(file_path_save, format='pdf')
plt.show()


if wavelength_integration_start == 280 and wavelength_integration_end == 315:
    plt.figure()
    plt.plot(integrated_flux_lst_water_specific, depth_list, linewidth=1)
    plt.title("Cockell Comparison: Depth vs. Flux at Depth (UVB only)")
    plt.xlabel("Flux, W/m^2")
    plt.ylabel("Depth, m")
    plt.xscale('log')
    plt.ylim(40, 0)
    plt.tight_layout()
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    file_path_save = os.path.join(directory, f"Cockell_Comparison_{timestamp}.pdf")
    plt.savefig(file_path_save, format='pdf')
    plt.show()

    plt.figure()
    plt.plot(depth_list, integrated_flux_lst_water_mixed, linewidth=1)
    plt.title("Flux in Mixed Layer vs. Mixed Layer Depth")
    plt.xlabel("Mixed Layer Depth, m")
    plt.ylabel("Average Flux, W/m^2")
    plt.tight_layout()
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    file_path_save = os.path.join(directory, f"Average_Flux_vs_Mixed_Layer_Depth_{mixed_layer_depth}m_{timestamp}.pdf")
    plt.savefig(file_path_save, format='pdf')
    plt.show()


# plt.figure()
# plt.plot(wavelengths[(200 <= wavelengths)], water_coefficients_lst_1[(200 <= wavelengths)], label=f'H2O Coeff. Data 1',
#          color='red')
# plt.plot(wavelengths[(200 <= wavelengths)], water_coefficients_lst_2[(200 <= wavelengths)], label=f'H2O Coeff. Data 2',
#          color='limegreen')
# plt.plot(wavelengths[(200 <= wavelengths)], water_coefficients_lst_3[(200 <= wavelengths)], label=f'H2O Coeff. Data 3',
#          color='blue')
# # plt.plot(wavelengths, water_coefficients_lst_4, label=f'H2O Coeff. Data 4', color='deeppink', linewidth=0.5)
# plt.title("Water Coefficients (from different data) vs. Wavelengths")
# plt.xlabel("Wavelength, nm")
# plt.ylabel("Water Coefficients, 1/m")
# plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.25), ncol=2)
# plt.tight_layout()
# file_path_save = os.path.join(directory, f"Water_Coefficients.pdf")
# plt.savefig(file_path_save, format='pdf')
# plt.show()

plt.figure()
plt.plot(wavelengths, flux_needed, label=f'Top of Atmosphere', color='#D55E00', linewidth=1)
plt.plot(wavelengths, flux_post_atmosphere, label=f"At Earth's Surface", color='#0072B2', linewidth=1)
plt.plot(wavelengths, flux_post_water_specific, label=f"Under {water_depth} m of water", color='black', linewidth=1)
plt.plot(wavelengths, flux_post_water_mixed, label=f"Ocean mixed layer of {mixed_layer_depth} m", color='purple', linewidth=1)
print(f'Flux after water: {flux_post_water_specific}')
plt.title(f'Irradiance vs. Wavelength (Sun age: {target_age} Gyr)')
plt.xlabel('Wavelength, nm')
plt.ylabel('Irradiance, W * m^-2 * nm^-1')
plt.yscale('log')

# region Bond-Breaking Wavelengths Lines
plt.axvline(x=single_CN_bond_wavelength, color='dimgrey', linewidth=1.75)
plt.text(single_CN_bond_wavelength, 10 ** (-1.5), "C-N", rotation=90, verticalalignment='center',
         horizontalalignment='right', color='black')

plt.axvline(x=single_PO_bond_wavelength, color='dimgrey', linewidth=1.75)
plt.text(single_PO_bond_wavelength, 10 ** (-2.5), "P-O", rotation=90, verticalalignment='center',
         horizontalalignment='left', color='black')

plt.axvline(x=single_CC_bond_wavelength, color='dimgrey', linewidth=1.75)
plt.text(single_CC_bond_wavelength, 10 ** (-2.5), "C-C", rotation=90, verticalalignment='center',
         horizontalalignment='left', color='black')

plt.axvline(x=single_CH_bond_wavelength, color='dimgrey', linewidth=1.75)
plt.text(single_CH_bond_wavelength, 10 ** (-1.5), "C-H", rotation=90, verticalalignment='center',
         horizontalalignment='left', color='black')

plt.axvline(x=single_CO_bond_wavelength, color='dimgrey', linewidth=1.75)
plt.text(single_CO_bond_wavelength, 10 ** (-2.5), "C-O", rotation=90, verticalalignment='center',
         horizontalalignment='left', color='black')

# plt.axvline(x=single_CP_bond_wavelength, color='dimgrey', linewidth=1.75)
# plt.text(single_CP_bond_wavelength, 10 ** (-0.5), "C-P", rotation=90, verticalalignment='center',
#          horizontalalignment='left', color='black')

plt.axvline(x=double_CC_bond_wavelength, color='dimgrey', linewidth=1.75)
plt.text(double_CC_bond_wavelength, 10 ** (-0.5), "C=C", rotation=90, verticalalignment='center',
         horizontalalignment='left', color='black')

plt.axvline(x=double_CO_bond_wavelength, color='dimgrey', linewidth=1.75)
plt.text(double_CO_bond_wavelength, 10 ** (-0.5), "C=O", rotation=90, verticalalignment='center',
         horizontalalignment='left', color='black')

if bond_break == "C-N":
    plt.axvline(x=single_CN_bond_wavelength, color='rebeccapurple', linewidth=4, alpha=0.4)
elif bond_break == "P-O":
    plt.axvline(x=single_PO_bond_wavelength, color='rebeccapurple', linewidth=4, alpha=0.4)
elif bond_break == "C-C":
    plt.axvline(x=single_CC_bond_wavelength, color='rebeccapurple', linewidth=4, alpha=0.4)
elif bond_break == "C-H":
    plt.axvline(x=single_CH_bond_wavelength, color='rebeccapurple', linewidth=4, alpha=0.4)
elif bond_break == "C-O":
    plt.axvline(x=single_CO_bond_wavelength, color='rebeccapurple', linewidth=4, alpha=0.4)
# elif bond_break == "C-P":
#     plt.axvline(x=single_CP_bond_wavelength, color='rebeccapurple', linewidth=4, alpha=0.4)
elif bond_break == "C=C":
    plt.axvline(x=double_CC_bond_wavelength, color='rebeccapurple', linewidth=4, alpha=0.4)
elif bond_break == "C=O":
    plt.axvline(x=double_CO_bond_wavelength, color='rebeccapurple', linewidth=4, alpha=0.4)
# endregion

# region Specific UV code
plt.fill_between(wavelengths, flux_needed, where=((wavelengths >= only_uva[0]) & (wavelengths <= only_uva[1])),
                 color='#F0E442', alpha=0.35, label="UVA (315nm-400nm)")
plt.fill_between(wavelengths, flux_needed, where=((wavelengths >= only_uvb[0]) & (wavelengths <= only_uvb[1])),
                 color='#009E73', alpha=0.35, label="UVB (280nm-315nm)")
plt.fill_between(wavelengths, flux_needed, where=((wavelengths >= only_uvc[0]) & (wavelengths <= only_uvc[1])),
                 color='#56B4E9', alpha=0.35, label="UVC (100nm-280nm)")
# endregion

plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.25), ncol=1)
plt.tight_layout()
timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
file_path_save = os.path.join(directory, f"Flux_vs_Wavelength_{timestamp}.pdf")
plt.savefig(file_path_save, format='pdf')
plt.show()

flux_vs_wavelength_data = dict()
flux_vs_wavelength_data['Wavelength, nm'] = wavelengths
flux_vs_wavelength_data[f'Flux at top of atmosphere, W/(m^2 * nm)'] = flux_needed
flux_vs_wavelength_data[f'Flux at Earth surface, W/(m^2 * nm)'] = flux_post_atmosphere
flux_vs_wavelength_data[f'Flux under {water_depth} m of water, W/(m^2 * nm)'] = flux_post_water_specific
flux_vs_wavelength_data[f'Flux in ocean mixed layer of {mixed_layer_depth} m, W/(m^2 * nm)'] = flux_post_water_mixed
flux_vs_wavelength_df = pd.DataFrame(flux_vs_wavelength_data)
flux_vs_wavelength_df.to_csv(rf'Results\Flux_vs_Wavelength_Data_{timestamp}.csv', index=False)


fig, (ax1, ax2) = plt.subplots(2, 1, sharex="all")
ax1.plot(time_list, integrated_flux_lst_solar, label=f"Top of atm", color='#7E2954', linewidth=1)
ax1.plot(time_list, integrated_flux_lst_atmosphere, label=f"At Earth's Surface", color='#337538', linewidth=1,
         linestyle='dashed')
ax1.plot(time_list, integrated_flux_lst_water, label=f"Under {water_depth} m of water", color='#2E2585',
         linewidth=1, linestyle='dashdot')
ax1.plot(time_list, integrated_flux_lst_water_mixed_time, label=f"Ocean mixed layer of {mixed_layer_depth} m",
         color='green', linewidth=1, linestyle='dotted')
time_list = np.array(time_list)
integrated_flux_lst_atmosphere = np.array(integrated_flux_lst_atmosphere)
integrated_flux_lst_water = np.array(integrated_flux_lst_water)
pre_goe_region = (time_list <= 2.203)
pre_ozone_region = (time_list >= 2.203) & (time_list <= 4.003)
modern_region = (time_list >= 4.003)
ax1.axvspan(0, 2.203, label=f"Pre GOE", color='#C26A77', alpha=0.35)
ax1.axvspan(2.203, 4.003, label=f"Post GOE", color='#DCCD7D', alpha=0.35)
ax1.axvspan(4.003, 5.6, label=f"Post GOE, modern O3", color='#94CBEC', alpha=0.35)
if bond_break != "":
    if 100 <= wavelength_integration_start and wavelength_integration_end <= 280:
        ax1.set_title(f'({bond_break}) UVC Flux ({wavelength_integration_start} nm to '
                      f'{wavelength_integration_end} nm) vs. Time', fontsize=12)
    elif 280 <= wavelength_integration_start and wavelength_integration_end <= 315:
        ax1.set_title(f'({bond_break}) UVB Flux ({wavelength_integration_start} nm to '
                      f'{wavelength_integration_end} nm) vs. Time', fontsize=12)
    elif 315 <= wavelength_integration_start and wavelength_integration_end <= 400:
        ax1.set_title(f'({bond_break}) UVA Flux ({wavelength_integration_start} nm to '
                      f'{wavelength_integration_end} nm) vs. Time', fontsize=12)
    else:
        ax1.set_title(f"({bond_break})  Flux ({wavelength_integration_start} to "
                      f"{wavelength_integration_end} nm) vs. Time", fontsize=12)
else:
    if 100 <= wavelength_integration_start and wavelength_integration_end <= 280:
        ax1.set_title(f'UVC Flux ({wavelength_integration_start} nm to {wavelength_integration_end} nm) vs. Time',
                      fontsize=12)
    elif 280 <= wavelength_integration_start and wavelength_integration_end <= 315:
        ax1.set_title(f'UVB Flux ({wavelength_integration_start} nm to {wavelength_integration_end} nm) vs. Time',
                      fontsize=12)
    elif 315 <= wavelength_integration_start and wavelength_integration_end <= 400:
        ax1.set_title(f'UVA Flux ({wavelength_integration_start} nm to {wavelength_integration_end} nm) vs. Time',
                      fontsize=12)
    else:
        ax1.set_title(f"Flux ({wavelength_integration_start} to {wavelength_integration_end} nm) vs. Time",
                      fontsize=12)
ax1.set_ylabel("Flux (in W * m^-2)")
ax1.legend(loc='center left', bbox_to_anchor=(1.05, 0.5), fontsize='small')
ax2.plot(time_list, N2_pressure_lst, label=f"N2 pressure", color='#FFB000')
ax2.plot(time_list, CO2_pressure_lst, label=f"CO2 pressure", color='#FE6100')
ax2.plot(time_list, O2_pressure_lst, label=f"O2 pressure", color='#648FFF')
ax2.plot(time_list, total_pressure_lst, label=f"Total atm pressure", color="#DC267F")
ax2.set_xlabel("Time (in Gyr)")
ax2.set_ylabel("Pressure (in bar)", fontsize='small')
ax2.set_title("Atmosphere Gas Pressures vs. Time")
ax2.legend(loc='center left', bbox_to_anchor=(1.2, 0.5), fontsize='small')
ax3 = ax2.twinx()
ax3.plot(time_list, O3_presence_lst, label="O3 P/g (plug)", color='#785EF0')
ax3.set_ylabel("O3 P/g (plug), kg/m^2", fontsize='small')
ax3.set_ylim(0, 0.007)
ax3.legend(loc='center left', bbox_to_anchor=(1.2, 0.2), fontsize='small')

plt.tight_layout()
timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
file_path_save = os.path.join(directory, f"Flux_vs_Time_{timestamp}.pdf")
plt.savefig(file_path_save, format='pdf')
plt.show()


plt.figure()
plt.plot(time_list, integrated_flux_lst_solar, label=f"Top of atm", color='#7E2954', linewidth=1.5)
plt.plot(time_list, integrated_flux_lst_atmosphere, label=f"At Earth's Surface", color='#337538', linewidth=1.5,
         linestyle='dashed')
plt.plot(time_list, integrated_flux_lst_water, label=f"Under {water_depth} m of water", color='#2E2585',
         linewidth=1.5, linestyle='dashdot')
plt.plot(time_list, integrated_flux_lst_water_mixed_time, label=f"Ocean mixed layer of {mixed_layer_depth} m",
         color='green', linewidth=1.5, linestyle='dotted')
time_list = np.array(time_list)
integrated_flux_lst_atmosphere = np.array(integrated_flux_lst_atmosphere)
integrated_flux_lst_water = np.array(integrated_flux_lst_water)
pre_goe_region = (time_list <= 2.203)
pre_ozone_region = (time_list >= 2.203) & (time_list <= 4.003)
modern_region = (time_list >= 4.003)
plt.axvspan(0, 2.203, label=f"Pre GOE", color='#C26A77', alpha=0.35)
plt.axvspan(2.203, 4.003, label=f"Post GOE", color='#DCCD7D', alpha=0.35)
plt.axvspan(4.003, 5.6, label=f"Post GOE, modern O3", color='#94CBEC', alpha=0.35)
if bond_break != "":
    if 100 <= wavelength_integration_start and wavelength_integration_end <= 280:
        plt.title(f'({bond_break}) UVC Flux ({wavelength_integration_start} nm to '
                      f'{wavelength_integration_end} nm) vs. Time', fontsize=12)
    elif 280 <= wavelength_integration_start and wavelength_integration_end <= 315:
        plt.title(f'({bond_break}) UVB Flux ({wavelength_integration_start} nm to '
                      f'{wavelength_integration_end} nm) vs. Time', fontsize=12)
    elif 315 <= wavelength_integration_start and wavelength_integration_end <= 400:
        plt.title(f'({bond_break}) UVA Flux ({wavelength_integration_start} nm to '
                      f'{wavelength_integration_end} nm) vs. Time', fontsize=12)
    else:
        plt.title(f"({bond_break})  Flux ({wavelength_integration_start} to "
                      f"{wavelength_integration_end} nm) vs. Time", fontsize=12)
else:
    if 100 <= wavelength_integration_start and wavelength_integration_end <= 280:
        plt.title(f'UVC Flux ({wavelength_integration_start} nm to {wavelength_integration_end} nm) vs. Time',
                      fontsize=12)
    elif 280 <= wavelength_integration_start and wavelength_integration_end <= 315:
        plt.title(f'UVB Flux ({wavelength_integration_start} nm to {wavelength_integration_end} nm) vs. Time',
                      fontsize=12)
    elif 315 <= wavelength_integration_start and wavelength_integration_end <= 400:
        plt.title(f'UVA Flux ({wavelength_integration_start} nm to {wavelength_integration_end} nm) vs. Time',
                      fontsize=12)
    else:
        plt.title(f"Flux ({wavelength_integration_start} to {wavelength_integration_end} nm) vs. Time",
                      fontsize=12)
plt.ylabel("Flux (in W * m^-2)")
plt.xlabel("Time (in Gyr)")
plt.legend(loc='center left', bbox_to_anchor=(1.05, 0.5), fontsize='small')

plt.tight_layout()
file_path_save = os.path.join(directory, f"Flux_vs_Time_no_atmosphere_{timestamp}.pdf")
plt.savefig(file_path_save, format='pdf')
plt.show()

flux_vs_time_data = dict()
flux_vs_time_data['Age of Sun and Earth, Gyr'] = time_list
flux_vs_time_data[f'Flux at top of atmosphere, W/(m^2)'] = integrated_flux_lst_solar
flux_vs_time_data[f'Flux at Earth surface, W/(m^2)'] = integrated_flux_lst_atmosphere
flux_vs_time_data[f'Flux under {water_depth} m of water, W/(m^2)'] = integrated_flux_lst_water
flux_vs_time_data[f'Flux in ocean mixed layer of {mixed_layer_depth} m, W/(m^2)'] = integrated_flux_lst_water_mixed_time
flux_vs_time_df = pd.DataFrame(flux_vs_time_data)
flux_vs_time_df.to_csv(rf'Results\Flux_vs_Time_Data_from_{wavelength_integration_start}_nm_to_{wavelength_integration_end}_nm_{timestamp}.csv', index=False)


plt.figure()
plt.plot(time_list, organic_integrated_flux_lst_atmosphere_grams,
         label=f"Mass of amino acids (or other organics) broken down (breaking C-N) at surface", color='black', linewidth=1)
plt.plot(time_list, organic_integrated_flux_lst_water_grams,
         label=f"Mass of amino acids (or other organics) under {water_depth} m of water",
         color='#0072B2', linewidth=1, linestyle="dashed")
plt.plot(time_list, organic_integrated_flux_lst_mixed_water_grams,
         label=f"Mass of amino acids (or other organics) in ocean mixed layer of {mixed_layer_depth} m",
         color='green', linewidth=1, linestyle="dashdot")  #TODO Fix some colouring for averages
time_list = np.array(time_list)
organic_integrated_flux_lst_atmosphere_grams = np.array(organic_integrated_flux_lst_atmosphere_grams)
organic_integrated_flux_lst_water_grams = np.array(organic_integrated_flux_lst_water_grams)
pre_goe_region = (time_list <= 2.203)
pre_ozone_region = (time_list >= 2.203) & (time_list <= 4.003)
modern_region = (time_list >= 4.003)
plt.axvspan(0, 2.203, label=f"Pre GOE", color='#C26A77', alpha=0.35)
plt.axvspan(2.203, 4.003, label=f"Post GOE", color='#DCCD7D', alpha=0.35)
plt.axvspan(4.003, 5.6, label=f"Post GOE, modern O3", color='#94CBEC', alpha=0.35)
plt.xlabel("Time (in Gyr)")
plt.ylabel("Broken organics flux (in g * s^-1 * m^-2)")
plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.25), ncol=1)
# plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.25), ncol=1)
# plt.scatter(point_x1, point_alpha1, color='navy')
# plt.text(point_x1, point_alpha1, f'({point_x1}, {point_alpha1:.2e})', color='black', size=10)
# plt.scatter(point_x2, point_alpha2, color='navy')
# plt.text(point_x2, point_alpha2, f'({point_x2}, {point_alpha2:.2e})', color='black', size=10)
# plt.scatter(point_x1, point_beta1, color='navy')
# plt.text(point_x1, point_beta1, f'({point_x1}, {point_beta1:.2e})', color='black', size=10)
# plt.scatter(point_x2, point_beta2, color='navy')
# plt.text(point_x2, point_beta2, f'({point_x2}, {point_beta2:.2e})', color='black', size=10)
plt.title("Flux of broken organic matter vs Time")
plt.tight_layout()
timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
file_path_save = os.path.join(directory, f"Organic_Flux_vs_Time_{timestamp}.pdf")
plt.savefig(file_path_save, format='pdf')
plt.show()
# endregion

# plt.figure()
# depths = np.array([0])
# for p in range(100):
#     depths = np.append(depths, depths[p]+0.1)
# uvb_integrated_flux = []
# uvb_mask = (280 <= wavelengths) & (wavelengths <= 315)
# uvb_wavelengths = wavelengths[uvb_mask]
# for p in range(len(depths)):
#     flux_water = np.array(atmosphere_to_water_flux(flux_post_atmosphere, main_water_coefficient_lst, depths[p]))
#     uvb_integrated_flux.append(simpson(y=flux_water[uvb_mask], x=uvb_wavelengths))
# plt.plot(depths, uvb_integrated_flux, color='magenta')
# plt.title("UVB Flux vs. Depth in Water (Sanity Check)")
# plt.xlabel("Depth, m")
# plt.ylabel("Flux, W / m^2")
# plt.tight_layout()
# file_path_save = os.path.join(directory, f"UVB_vs_Depth.pdf")
# plt.savefig(file_path_save, format='pdf')
# plt.show()
#
# plt.figure()
# depths = np.array([0])
# for p in range(700):
#     depths = np.append(depths, depths[p]+0.1)
# uva_integrated_flux = []
# uva_mask = (315 <= wavelengths) & (wavelengths <= 400)
# uva_wavelengths = wavelengths[uva_mask]
# for p in range(len(depths)):
#     flux_water = np.array(atmosphere_to_water_flux(flux_post_atmosphere, main_water_coefficient_lst, depths[p]))
#     uva_integrated_flux.append(simpson(y=flux_water[uva_mask], x=uva_wavelengths))
# plt.plot(depths, uva_integrated_flux, color='blue')
# plt.title("UVA Flux vs. Depth in Water (Sanity Check)")
# plt.xlabel("Depth, m")
# plt.ylabel("Flux, W / m^2")
# plt.tight_layout()
# file_path_save = os.path.join(directory, f"UVA_vs_Depth.pdf")
# plt.savefig(file_path_save, format='pdf')
# plt.show()

fig, ax1 = plt.subplots(1, 1)
ax1.plot(time_list, O2_pressure_lst, label=f"O2 pressure", color='#648FFF')
ax1.set_ylabel("Pressure of O2, bar")
ax1.set_xlabel("Time, Gyr")
ax1.set_title("Oxygen-related Gas Check")
ax2 = ax1.twinx()
ax2.plot(time_list, O3_presence_lst, label="O3 P/g (plug)", color='#785EF0')
ax2.set_ylabel("O3 P/g (plug), kg/m^2")
ax2.set_ylim(0, 0.007)
plt.tight_layout()
timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
file_path_save = os.path.join(directory, f"Oxygen_Amounts_and_Ozone_Plug_vs_Time_{timestamp}.pdf")
plt.savefig(file_path_save, format='pdf')
plt.show()

# region Final prints
print(f"Bolometric flux for {target_age} Gyr: {flux_bolometric_calc(target_age)} W/m^2.")
end_time = time.time()
print(f"Elapsed time: {end_time - start_time} seconds")
# endregion
