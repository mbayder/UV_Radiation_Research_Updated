import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline, interp1d
from scipy.constants import h, c
from math import e
from scipy.integrate import simpson
import time
import os
import unicodeit
from datetime import datetime


def smoothing_spline_interpolation_solar(age, age_lst, flux_lst, s):
    interpolated_flux = np.zeros(flux_lst.shape[1])
    for i in range(flux_lst.shape[1]):
        spline = UnivariateSpline(age_lst, flux_lst[:, i], s=s)
        interpolated_flux[i] = spline(age)
    return interpolated_flux


def solar_to_atmosphere_flux_photons(flux_given, tau_lst):
    flux_post_atm = []
    for m in range(len(flux_given)):
        flux_post_atm_item = flux_given[m] * e ** ((-1) * tau_lst[m])
        if flux_post_atm_item >= 1000:
            flux_post_atm.append(flux_post_atm_item)
        else:
            flux_post_atm.append(0.0)
    return flux_post_atm


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
print("Note: Assuming 0.45 g of carbon per g of amino acid.")

double_CO_bond_wavelength = 161.0  # in nm
double_CC_bond_wavelength = 199.0  # in nm
single_CH_bond_wavelength = 291.0  # in nm
single_CO_bond_wavelength = 334.0  # in nm
single_CC_bond_wavelength = 346.0  # in nm
single_PO_bond_wavelength = 357.0  # in nm
single_CN_bond_wavelength = 392.0  # in nm

print("Wavelengths needed to break each bond below")
print(" C=O: 161 nm \nC=C: 195 nm \nC-H: 290 nm \nC-O: 334 nm \nC-C: 345 nm \nP-O: 357 nm \nC-N: 392 nm")

mixed_layer_depth = 50
organic_mask_start_1 = 100
organic_mask_end_1 = single_CC_bond_wavelength
organic_mask_start_2 = 100
organic_mask_end_2 = single_CN_bond_wavelength
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
smoothing_factor = 0
# endregion

# region Basic data for atmosphere
mass_CO2 = 7.308032317 * 10 ** (-26)  # kg
mass_N2 = 4.652830391 * 10 ** (-26)
mass_O2 = 5.313724929 * 10 ** (-26)
mass_03 = 7.970587394 * 10 ** (-26)
g_Earth = 9.8  # N/kg
# endregion

# region Reading Atmosphere files
file_path_CO2_acs = r'C:\Users\User\Documents\McGill Internship 2024\ACS_Data\CO2_ACS.csv'
data_CO2_acs = pd.read_csv(file_path_CO2_acs, comment='#', names=['wavelength', 'acs'], nrows=322)
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
    return total_tau_lst


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
# endregion

# region Wavelength masks
organic_wavelength_mask_1 = (wavelengths >= organic_mask_start_1) & (wavelengths <= organic_mask_end_1)
organic_wavelength_interval_1 = wavelengths[organic_wavelength_mask_1]
water_coefficient_lst_org_1 = water_coefficients_lst_1[organic_wavelength_mask_1]

organic_wavelength_mask_2 = (wavelengths >= organic_mask_start_2) & (wavelengths <= organic_mask_end_2)
organic_wavelength_interval_2 = wavelengths[organic_wavelength_mask_2]
water_coefficient_lst_org_2 = water_coefficients_lst_1[organic_wavelength_mask_2]
# endregion


def organic_carbon_flux_mixedlayer_O2history(mixed_layer_depth, O2_history, organic_wavelength_mask):
    time_st = 0.603  # in Gyr
    time_lst = [time_st]
    time_change = 0.05  # in Gyr
    max_time = 4.603  # in Gyr

    organic_flux_interval_photons = smoothing_spline_interpolation_solar(time_st, ages_of_Sun, flux_data_photons,
                                                                         smoothing_factor)[organic_wavelength_mask]

    organic_integrated_flux_mixed_water_photons = simpson(
        y=atmosphere_to_mixed_water_flux_photons(solar_to_atmosphere_flux_photons(
            organic_flux_interval_photons, tau_value(time_st, O2_history)[organic_wavelength_mask]),
            water_coefficients_lst_1[organic_wavelength_mask],
            mixed_layer_depth),
        x=wavelengths[organic_wavelength_mask])
    organic_carbon_integrated_flux_mixed_water_grams = 0.45 * 110 * 300 * 1.66054 * 10 ** (-24) * \
                                                       organic_integrated_flux_mixed_water_photons / 300  # TODO fix numbers
    organic_carbon_integrated_flux_lst_mixed_water_grams = [organic_carbon_integrated_flux_mixed_water_grams]

    while time_st < max_time:
        time_st += time_change
        time_lst.append(time_st)

        organic_flux_interval_photons = smoothing_spline_interpolation_solar(time_st, ages_of_Sun, flux_data_photons,
                                                                             smoothing_factor)[organic_wavelength_mask]

        organic_integrated_flux_mixed_water_photons = simpson(
            y=atmosphere_to_mixed_water_flux_photons(solar_to_atmosphere_flux_photons(
                organic_flux_interval_photons, tau_value(time_st, O2_history)[organic_wavelength_mask]),
                water_coefficients_lst_1[organic_wavelength_mask],
                mixed_layer_depth),
            x=wavelengths[organic_wavelength_mask])
        organic_carbon_integrated_flux_mixed_water_grams = 0.45 * 110 * 300 * 1.66054 * 10 ** (-24) * \
                                                           organic_integrated_flux_mixed_water_photons / 300  # TODO fix numbers
        organic_carbon_integrated_flux_lst_mixed_water_grams.append(organic_carbon_integrated_flux_mixed_water_grams)

    results = {"time_lst": time_lst,
               "organic_carbon_integrated_flux_mixed_water_grams": organic_carbon_integrated_flux_lst_mixed_water_grams}
    return results


# region Calculating flux values (not normalized)
time_list = organic_carbon_flux_mixedlayer_O2history(mixed_layer_depth, "average", organic_wavelength_mask_1)["time_lst"]

organic_carbon_flux_organic1_lowerO2 = organic_carbon_flux_mixedlayer_O2history(mixed_layer_depth, "lower", organic_wavelength_mask_1)["organic_carbon_integrated_flux_mixed_water_grams"]
organic_carbon_flux_organic1_averageO2 = organic_carbon_flux_mixedlayer_O2history(mixed_layer_depth, "average", organic_wavelength_mask_1)["organic_carbon_integrated_flux_mixed_water_grams"]
organic_carbon_flux_organic1_upperO2 = organic_carbon_flux_mixedlayer_O2history(mixed_layer_depth, "upper", organic_wavelength_mask_1)["organic_carbon_integrated_flux_mixed_water_grams"]

organic_carbon_flux_organic2_lowerO2 = organic_carbon_flux_mixedlayer_O2history(mixed_layer_depth, "lower", organic_wavelength_mask_2)["organic_carbon_integrated_flux_mixed_water_grams"]
organic_carbon_flux_organic2_averageO2 = organic_carbon_flux_mixedlayer_O2history(mixed_layer_depth, "average", organic_wavelength_mask_2)["organic_carbon_integrated_flux_mixed_water_grams"]
organic_carbon_flux_organic2_upperO2 = organic_carbon_flux_mixedlayer_O2history(mixed_layer_depth, "upper", organic_wavelength_mask_2)["organic_carbon_integrated_flux_mixed_water_grams"]
# endregion

# region Creating directory
directory = r'C:\Users\User\Documents\McGill Internship 2024\PDF_Plots'
if not os.path.exists(directory):
    os.makedirs(directory)
# endregion

# region Plotting itself
plt.figure(figsize=(10, 9))
plt.plot([4.603-x for x in time_list], organic_carbon_flux_organic1_lowerO2, color='red', linewidth=3, linestyle="dotted")
plt.plot([4.603-x for x in time_list], organic_carbon_flux_organic1_averageO2, color='red', linewidth=3)
plt.plot([4.603-x for x in time_list], organic_carbon_flux_organic1_averageO2, color='red', linewidth=10, alpha=0.5)
plt.plot([4.603-x for x in time_list], organic_carbon_flux_organic1_upperO2, color='red', linewidth=3, linestyle="dashed")

plt.plot([4.603-x for x in time_list], organic_carbon_flux_organic2_lowerO2, color='#0072B2', linewidth=3, linestyle="dotted")
plt.plot([4.603-x for x in time_list], organic_carbon_flux_organic2_averageO2, color='#0072B2', linewidth=3)
plt.plot([4.603-x for x in time_list], organic_carbon_flux_organic2_averageO2, color='#0072B2', linewidth=10, alpha=0.5)
plt.plot([4.603-x for x in time_list], organic_carbon_flux_organic2_upperO2, color='#0072B2', linewidth=3, linestyle="dashed")

# plt.axvline(x=4.603-2, color='#D55E00', linewidth=5)
# plt.axvline(x=4.603-3.5, color='#009E73', linewidth=5)
# plt.axvline(x=4.603-4.603, color='#7E2954', linewidth=5)

plt.axvspan(4.603-0, 4.603-2.203, color='#D55E00', alpha=0.35)
plt.axvspan(4.603-2.203, 4.603-4.003, color='#009E73', alpha=0.35)
plt.axvspan(4.603-4.003, -0.1, color='#7E2954', alpha=0.35)

plt.gca().set_xlim(abs(min(time_list)-4.603), -0.1)

plt.xlabel("Time, Ga", fontsize=25)
plt.ylabel(f"Organic carbon flux, {unicodeit.replace("g s^{-1} m^{-2}")}", fontsize=25)
plt.tick_params(axis='both', labelsize=25)

plt.tight_layout()
timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
file_path_save = os.path.join(directory, f"Flux_vs_Time_Organic_{timestamp}.pdf")
plt.savefig(file_path_save, format='pdf')
plt.show()
# endregion

end_time = time.time()
print(f"Elapsed time: {end_time - start_time} seconds or {(end_time - start_time) / 60} minutes.")
