import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from scipy.interpolate import UnivariateSpline, interp1d
from scipy.constants import h, c
from math import e
import time
import os
import unicodeit
from datetime import datetime


plt.rcParams["font.family"] = ["Arial", "DejaVu Sans"]


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


# region Intro Print, Input Statements, Time Start
print("Note: UVC is 100nm to 280nm; UVB is 280nm to 315nm; UVA is 315nm to 400nm.")
print("Disclaimer: For photons units, numbers less than 1000 will be considered negligible (rounded to 0).")
print("Disclaimer: The bond energy of the C-N bond is for an amide (peptide) bond, not a normal C-N bond.")
print("Note: Taking average mass of 1 amino acid as 110 Da.")
print("Note: Taking an average of 300 amino acids per each chain.")

double_CO_bond_wavelength = 161.0  # in nm
double_CC_bond_wavelength = 199.0  # in nm
single_CH_bond_wavelength = 291.0  # in nm
single_CO_bond_wavelength = 334.0  # in nm
single_CC_bond_wavelength = 346.0  # in nm
single_PO_bond_wavelength = 357.0  # in nm
single_CN_bond_wavelength = 392.0  # in nm

# target_age_1 = float(input("Age of the Sun (in Gyr) #1 : "))
# target_age_2 = float(input("Age of the Sun (in Gyr) #2 : "))
# target_age_3 = float(input("Age of the Sun (in Gyr) #3 : "))
target_age_1 = 2
target_age_2 = 3.5
target_age_3 = 4.603
mixed_layer_depth = 50
start_time = time.time()
# endregion

# region Reading data files for Sun
file_path_00 = r'C:\Users\User\Documents\McGill Internship 2024\Claire_Data\SunModern.txt'
data_00 = pd.read_csv(file_path_00, sep='\\s+', comment='#', names=['wavelength', 'flux'],
                      skiprows=103, nrows=5726)

file_path_06 = r'C:\Users\User\Documents\McGill Internship 2024\Claire_Data\Sun0.6Ga.txt'
data_06 = pd.read_csv(file_path_06, sep='\\s+', comment='#', names=['wavelength', 'flux'],
                      skiprows=103, nrows=5726)

file_path_18 = r'C:\Users\User\Documents\McGill Internship 2024\Claire_Data\Sun1.8Ga.txt'
data_18 = pd.read_csv(file_path_18, sep='\\s+', comment='#', names=['wavelength', 'flux'],
                      skiprows=103, nrows=5726)

file_path_24 = r'C:\Users\User\Documents\McGill Internship 2024\Claire_Data\Sun2.4Ga.txt'
data_24 = pd.read_csv(file_path_24, sep='\\s+', comment='#', names=['wavelength', 'flux'],
                      skiprows=103, nrows=5726)

file_path_38 = r'C:\Users\User\Documents\McGill Internship 2024\Claire_Data\Sun3.8Ga.txt'
data_38 = pd.read_csv(file_path_38, sep='\\s+', comment='#', names=['wavelength', 'flux'],
                      skiprows=103, nrows=5726)

file_path_44 = r'C:\Users\User\Documents\McGill Internship 2024\Claire_Data\Sun4.4Ga.txt'
data_44 = pd.read_csv(file_path_44, sep='\\s+', comment='#', names=['wavelength', 'flux'],
                      skiprows=103, nrows=5726)
# endregion

# region Creating wavelengths and fluxes lists from data
wavelengths = np.array(data_00['wavelength'])
wavelengths_conv = 10 ** (-9) * np.array(data_06['wavelength'])
flux_00 = h * c * 10 ** 4 * np.array(data_00['flux']) / wavelengths_conv
flux_06 = h * c * 10 ** 4 * np.array(data_06['flux']) / wavelengths_conv
flux_18 = h * c * 10 ** 4 * np.array(data_18['flux']) / wavelengths_conv
flux_24 = h * c * 10 ** 4 * np.array(data_24['flux']) / wavelengths_conv
flux_38 = h * c * 10 ** 4 * np.array(data_38['flux']) / wavelengths_conv
flux_44 = h * c * 10 ** 4 * np.array(data_44['flux']) / wavelengths_conv

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
ages_of_Sun = np.array([0.203, 0.803, 2.203, 2.803, 4.003, 4.603])
flux_data = np.vstack([flux_44, flux_38, flux_24, flux_18, flux_06, flux_00])
flux_data_photons = np.vstack([flux_44_photons, flux_38_photons, flux_24_photons, flux_18_photons,
                               flux_06_photons, flux_00_photons])


def top_atmosphere_flux(target_age):
    flux_needed = []
    smoothing_factor = 0
    if target_age in ages_of_Sun:
        if target_age == 0.203:
            flux_needed = flux_44
        if target_age == 0.803:
            flux_needed = flux_38
        if target_age == 2.203:
            flux_needed = flux_24
        if target_age == 2.803:
            flux_needed = flux_18
        if target_age == 4.003:
            flux_needed = flux_06
        if target_age == 4.603:
            flux_needed = flux_00
    else:
        flux_needed = smoothing_spline_interpolation_solar(target_age, ages_of_Sun, flux_data, smoothing_factor)
    return flux_needed


flux_needed_1 = top_atmosphere_flux(target_age_1)
flux_needed_2 = top_atmosphere_flux(target_age_2)
flux_needed_3 = top_atmosphere_flux(target_age_3)
# endregion

# region Basic Atmosphere Data
mass_CO2 = 7.308032317 * 10 ** (-26)  # kg
mass_N2 = 4.652830391 * 10 ** (-26)
mass_O2 = 5.313724929 * 10 ** (-26)
mass_03 = 7.970587394 * 10 ** (-26)
g_Earth = 9.8  # N/kg
# endregion

# region Reading Atmosphere files
file_path_CO2_acs = r'C:\Users\User\Documents\McGill Internship 2024\ACS_Data\CO2_ACS.csv'
data_CO2_acs = pd.read_csv(file_path_CO2_acs, comment='#', names=['wavelength', 'acs'], nrows=322)
# data_CO2_acs = data_CO2_acs.drop_duplicates(subset='wavelength')
data_CO2_acs_wl = data_CO2_acs['wavelength'].astype(float)
data_CO2_acs_acs = data_CO2_acs['acs'].astype(float)

file_path_N2_acs = r'C:\Users\User\Documents\McGill Internship 2024\ACS_Data\N2_ACS.csv'
data_N2_acs = pd.read_csv(file_path_N2_acs, comment='#', names=['wavelength', 'acs'], nrows=12)
data_N2_acs_wl = data_N2_acs['wavelength'].astype(float)
data_N2_acs_acs = data_N2_acs['acs'].astype(float)

file_path_O2_acs = r'C:\Users\User\Documents\McGill Internship 2024\ACS_Data\O2_ACS.csv'
data_O2_acs = pd.read_csv(file_path_O2_acs, comment='#', names=['wavelength', 'acs'], nrows=322)
data_O2_acs_wl = data_O2_acs['wavelength'].astype(float)
data_O2_acs_acs = data_O2_acs['acs'].astype(float)

file_path_O2_time_lower = r'C:\Users\User\Documents\McGill Internship 2024\ACS_Data\O2_vs_Time_Lower.csv'
data_O2_time_lower = pd.read_csv(file_path_O2_time_lower, delimiter=",", comment='#', names=['time', 'pressure'],
                                 nrows=180)
data_O2_time_time_lower = data_O2_time_lower['time'].astype(float)
data_O2_time_pressure_lower = data_O2_time_lower['pressure'].astype(float)

file_path_O2_time_upper = r'C:\Users\User\Documents\McGill Internship 2024\ACS_Data\O2_vs_Time_Upper.csv'
data_O2_time_upper = pd.read_csv(file_path_O2_time_upper, delimiter=",", comment='#', names=['time', 'pressure'],
                                 nrows=174)
data_O2_time_time_upper = data_O2_time_upper['time'].astype(float)
data_O2_time_pressure_upper = data_O2_time_upper['pressure'].astype(float)

file_path_O2_time_average = r'C:\Users\User\Documents\McGill Internship 2024\ACS_Data\O2_vs_Time.csv'
data_O2_time_average = pd.read_csv(file_path_O2_time_average, delimiter=",", comment='#', names=['time', 'pressure'],
                                   nrows=177)
data_O2_time_time_average = data_O2_time_average['time'].astype(float)
data_O2_time_pressure_average = data_O2_time_average['pressure'].astype(float)

file_path_O3_acs = r'C:\Users\User\Documents\McGill Internship 2024\ACS_Data\O3_ACS.csv'
data_O3_acs = pd.read_csv(file_path_O3_acs, comment='#', names=['wavelength', 'acs'], nrows=562)
data_O3_acs_wl = data_O3_acs['wavelength'].astype(float)
data_O3_acs_acs_special = data_O3_acs['acs'].astype(float)

file_path_O3_sigma = r'C:\Users\User\Documents\McGill Internship 2024\ACS_Data\O3_Sigma.csv'
data_O3_sigma = pd.read_csv(file_path_O3_sigma, comment="#", names=['O2_amount', 'O3_sigma'])
data_O3_sigma_O2 = data_O3_sigma['O2_amount'].astype(float)
data_O3_sigma_sigma = data_O3_sigma['O3_sigma'].astype(float)
data_O3_sigma_sigma = data_O3_sigma_sigma * 2.687 * 10 ** 20 * mass_03
# endregion

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

# region Perform spline interpolation on the log-transformed absorption cross-section data -> Getting tau list and pressures


def tau_value(t_age, O2_history):
    if "l" in O2_history:
        data_O2_time_time = data_O2_time_time_lower
        data_O2_time_pressure = data_O2_time_pressure_lower
    elif "u" in O2_history:
        data_O2_time_time = data_O2_time_time_upper
        data_O2_time_pressure = data_O2_time_pressure_upper
    else:
        data_O2_time_time = data_O2_time_time_average
        data_O2_time_pressure = data_O2_time_pressure_average
    if t_age < 2.203:  # Articles about CO2 and pressure
        pressure_past = 23000  # Pa
        percentage_CO2_past = 0.70
        pressure_CO2_past = pressure_past * percentage_CO2_past  # Pa
        percentage_N2_past = 0.30
        pressure_N2_past = pressure_past * percentage_N2_past  # Pa
        CO2_tau_lst = CO2_acs_lst * 10 ** (-4) * pressure_CO2_past / (mass_CO2 * g_Earth)
        N2_tau_lst = N2_acs_lst * 10 ** (-4) * pressure_N2_past / (mass_N2 * g_Earth)
        total_tau_lst = N2_tau_lst + CO2_tau_lst
    else:  # Articles / general data about pressure and gases
        if t_age >= 4.603:
            t_age = 4.603
        else:
            t_age = t_age
        O2_time_interpolator = interp1d(data_O2_time_time, data_O2_time_pressure,
                                        kind='linear')  # TODO: Define perhaps as a function
        pressure_O2_past = O2_time_interpolator(4.603 - t_age) * 21222.35
        # print(f"Past pressure at {4.603 - t_age} Ga: {pressure_O2_past}")
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
    return total_tau_lst


tau_lst_1_lower = tau_value(target_age_1, "lower")
tau_lst_1_average = tau_value(target_age_1, "average")
tau_lst_1_upper = tau_value(target_age_1, "upper")
tau_lst_2_lower = tau_value(target_age_2, "lower")
tau_lst_2_average = tau_value(target_age_2, "average")
tau_lst_2_upper = tau_value(target_age_2, "upper")
tau_lst_3_lower = tau_value(target_age_3, "lower")
tau_lst_3_average = tau_value(target_age_3, "average")
tau_lst_3_upper = tau_value(target_age_3, "upper")
# endregion

# region Find flux at Earth surface
flux_post_atmosphere_1_lower = solar_to_atmosphere_flux(flux_needed_1, tau_lst_1_lower)
flux_post_atmosphere_1_average = solar_to_atmosphere_flux(flux_needed_1, tau_lst_1_average)
flux_post_atmosphere_1_upper = solar_to_atmosphere_flux(flux_needed_1, tau_lst_1_upper)
flux_post_atmosphere_2_lower = solar_to_atmosphere_flux(flux_needed_2, tau_lst_2_lower)
flux_post_atmosphere_2_average = solar_to_atmosphere_flux(flux_needed_2, tau_lst_2_average)
flux_post_atmosphere_2_upper = solar_to_atmosphere_flux(flux_needed_2, tau_lst_2_upper)
flux_post_atmosphere_3_lower = solar_to_atmosphere_flux(flux_needed_3, tau_lst_3_lower)
flux_post_atmosphere_3_average = solar_to_atmosphere_flux(flux_needed_3, tau_lst_3_average)
flux_post_atmosphere_3_upper = solar_to_atmosphere_flux(flux_needed_3, tau_lst_3_upper)
# endregion

# region Ocean Water Modeling
file_path_water_1 = r'C:\Users\User\Documents\McGill Internship 2024\Water_Data\Water_Data_1.csv'

data_water_1 = pd.read_csv(file_path_water_1, comment='#', names=['wavelength', 'coefficient'],
                           nrows=157)  # For Water_Data_1

water_wavelengths_1 = data_water_1['wavelength'].astype(float)

water_coefficients_1 = data_water_1['coefficient'].astype(float)

water_interpolator_1 = interp1d(water_wavelengths_1, water_coefficients_1, kind='linear', bounds_error=False,
                                fill_value="extrapolate")
water_coefficients_lst_1 = np.array(water_interpolator_1(wavelengths))

flux_post_water_mixed_1_lower = atmosphere_to_mixed_water_flux(flux_post_atmosphere_1_lower, water_coefficients_lst_1, mixed_layer_depth)
flux_post_water_mixed_1_average = atmosphere_to_mixed_water_flux(flux_post_atmosphere_1_average, water_coefficients_lst_1, mixed_layer_depth)
flux_post_water_mixed_1_upper = atmosphere_to_mixed_water_flux(flux_post_atmosphere_1_upper, water_coefficients_lst_1, mixed_layer_depth)
flux_post_water_mixed_2_lower = atmosphere_to_mixed_water_flux(flux_post_atmosphere_2_lower, water_coefficients_lst_1, mixed_layer_depth)
flux_post_water_mixed_2_average = atmosphere_to_mixed_water_flux(flux_post_atmosphere_2_average, water_coefficients_lst_1, mixed_layer_depth)
flux_post_water_mixed_2_upper = atmosphere_to_mixed_water_flux(flux_post_atmosphere_2_upper, water_coefficients_lst_1, mixed_layer_depth)
flux_post_water_mixed_3_lower = atmosphere_to_mixed_water_flux(flux_post_atmosphere_3_lower, water_coefficients_lst_1, mixed_layer_depth)
flux_post_water_mixed_3_average = atmosphere_to_mixed_water_flux(flux_post_atmosphere_3_average, water_coefficients_lst_1, mixed_layer_depth)
flux_post_water_mixed_3_upper = atmosphere_to_mixed_water_flux(flux_post_atmosphere_3_upper, water_coefficients_lst_1, mixed_layer_depth)
# endregion

# region Directory making
directory = r'C:\Users\User\Documents\McGill Internship 2024\PDF_Plots'

if not os.path.exists(directory):
    os.makedirs(directory)
# endregion

# region Plotting itself
plt.figure(figsize=(10, 9))
plt.plot(wavelengths, flux_post_water_mixed_1_lower, color='#D55E00', linewidth=1, linestyle="dotted")
plt.plot(wavelengths, flux_post_water_mixed_1_average, color='#D55E00', linewidth=1)
plt.plot(wavelengths, flux_post_water_mixed_1_average, color='#D55E00', linewidth=8, alpha=0.5)
plt.plot(wavelengths, flux_post_water_mixed_1_upper, color='#D55E00', linewidth=1, linestyle='dashed')
plt.plot(wavelengths, flux_post_water_mixed_2_lower, color='#009E73', linewidth=1, linestyle="dotted")
plt.plot(wavelengths, flux_post_water_mixed_2_average, color='#009E73', linewidth=1)
plt.plot(wavelengths, flux_post_water_mixed_2_average, color='#009E73', linewidth=8, alpha=0.5)
plt.plot(wavelengths, flux_post_water_mixed_2_upper, color='#009E73', linewidth=1, linestyle='dashed')
plt.plot(wavelengths, flux_post_water_mixed_3_lower, color='#7E2954', linewidth=1, linestyle="dotted")
plt.plot(wavelengths, flux_post_water_mixed_3_average, color='#7E2954', linewidth=1)
plt.plot(wavelengths, flux_post_water_mixed_3_average, color='#7E2954', linewidth=8, alpha=0.5)
plt.plot(wavelengths, flux_post_water_mixed_3_upper, color='#7E2954', linewidth=1, linestyle='dashed')
plt.xlabel('Wavelength (nm)', fontsize=17)
plt.ylabel(f'Irradiance ({unicodeit.replace("W m^{-2} nm^{-1}")})', fontsize=17)
plt.yscale('log')
plt.tick_params(axis='both', labelsize=16)
plt.plot([], [], label="Archean", color='#D55E00', linewidth=3)
plt.plot([], [], label="Proterozoic", color='#009E73', linewidth=3)
plt.plot([], [], label="Phanerozoic", color='#7E2954', linewidth=3)
plt.plot([], [], label=f"Min {unicodeit.replace("O_2")} pressure values", color='black', linewidth=3, linestyle="dotted")
plt.plot([], [], label=f"Max {unicodeit.replace("O_2")} pressure values", color='black', linewidth=3, linestyle='dashed')
plt.plot([], [], label=f"Average {unicodeit.replace("O_2")} pressure values", color='black', linewidth=3)
plt.legend(loc='best', fontsize=14, frameon=False)

# plt.axvline(x=single_CN_bond_wavelength, color='#0072B2', linewidth=1.5)
# plt.text(single_CN_bond_wavelength, 10 ** (-1.5), "C-N", rotation=90, verticalalignment='center',
#          horizontalalignment='right', color='black', fontsize=27)
# plt.axvline(x=single_CN_bond_wavelength, color='#0072B2', linewidth=10, alpha=0.5)
#
# plt.axvline(x=single_CC_bond_wavelength, color='red', linewidth=1.5)
# plt.text(single_CC_bond_wavelength, 10 ** (-3.5), "C-C", rotation=90, verticalalignment='center',
#          horizontalalignment='right', color='black', fontsize=27)
# plt.axvline(x=single_CC_bond_wavelength, color='red', linewidth=10, alpha=0.5)
#
# plt.axvline(x=single_PO_bond_wavelength, color='black', linewidth=1.5)
# plt.text(single_PO_bond_wavelength, 10 ** (-4.5), "P-O", rotation=90, verticalalignment='center',
#          horizontalalignment='left', color='black', fontsize=27)
#
# plt.axvline(x=single_CH_bond_wavelength, color='black', linewidth=1.5)
# plt.text(single_CH_bond_wavelength, 10 ** (-0.5), "C-H", rotation=90, verticalalignment='center',
#          horizontalalignment='left', color='black', fontsize=27)
#
# plt.axvline(x=single_CO_bond_wavelength, color='black', linewidth=1.5)
# plt.text(single_CO_bond_wavelength, 10 ** (-2.5), "C-O", rotation=90, verticalalignment='center',
#          horizontalalignment='right', color='black', fontsize=27)
#
# plt.axvline(x=double_CC_bond_wavelength, color='black', linewidth=1.5)
# plt.text(double_CC_bond_wavelength, 10 ** (-0.5), "C=C", rotation=90, verticalalignment='center',
#          horizontalalignment='left', color='black', fontsize=27)
#
# plt.axvline(x=double_CO_bond_wavelength, color='black', linewidth=1.5)
# plt.text(double_CO_bond_wavelength, 10 ** (-0.5), "C=O", rotation=90, verticalalignment='center',
#          horizontalalignment='left', color='black', fontsize=27)

# plt.axvspan(only_uvc[0], only_uvc[1], label=f"UVC (100nm-280nm)", color='#7F66FF', alpha=0.35)
# plt.axvspan(only_uvb[0], only_uvb[1], label=f"UVB (280nm-315nm)", color='#E3006A', alpha=0.35)
# plt.axvspan(only_uva[0], only_uva[1], label=f"UVA (315nm-400nm), modern O3", color='#FFA500', alpha=0.35)

ax = plt.gca()

y_bottom = 0
y_top = 0.5 * 10**(-5)
height = y_top - y_bottom

uvc_patch = patches.Rectangle((100, y_bottom), 180, height, edgecolor="none", facecolor="dimgrey", alpha=0.6)
ax.add_patch(uvc_patch)
ax.text(190, 10**(-5.48), "UVC", ha='center', va='center', fontsize=11, fontweight='bold')

uvb_patch = patches.Rectangle((280, y_bottom), 35, height, edgecolor="none", facecolor="silver", alpha=0.8)
ax.add_patch(uvb_patch)
ax.text(297.5, 10**(-5.48), "UVB", ha='center', va='center', fontsize=11, fontweight='bold')

uva_patch = patches.Rectangle((315, y_bottom), 85, height, edgecolor="none", facecolor="gainsboro", alpha=0.6)
ax.add_patch(uva_patch)
ax.text(357.5, 10**(-5.48), "UVA", ha='center', va='center', fontsize=11, fontweight='bold')

plt.tight_layout()
timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
file_path_save = os.path.join(directory, f"Irradiance_vs_Wavelength_Mixed_Water_{timestamp}.pdf")
plt.savefig(file_path_save, format='pdf')
file_path_save = os.path.join(directory, f"Irradiance_vs_Wavelength_Mixed_Water_{timestamp}.png")
plt.savefig(file_path_save, format='png', dpi=1200)
plt.show()
# endregion

end_time = time.time()
print(f"Elapsed time: {end_time - start_time} seconds")
