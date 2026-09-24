import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from scipy.interpolate import UnivariateSpline
from scipy.constants import h, c
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
start_time = time.time()

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

only_uva = [315, 400]
until_uva = [100, 400]
only_uvb = [280, 315]
until_uvb = [100, 315]
only_uvc = [100, 280]
until_uvc = [100, 280]

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

directory = r'Results'

if not os.path.exists(directory):
    os.makedirs(directory)

plt.figure(figsize=(10, 5))
plt.plot(wavelengths, flux_needed_1, color='#D55E00', linewidth=1)
plt.plot(wavelengths, flux_needed_2, color='#009E73', linewidth=1)
plt.plot(wavelengths, flux_needed_3, color='#7E2954', linewidth=1)
plt.xlabel('Wavelength (nm)', fontsize=17)
plt.ylabel(f'Irradiance ({unicodeit.replace("W m^{-2} nm^{-1}")})', fontsize=17)
plt.yscale('log')
plt.tick_params(axis='both', labelsize=16)
plt.plot([], [], label="Archean", color='#D55E00', linewidth=3)
plt.plot([], [], label="Proterozoic", color='#009E73', linewidth=3)
plt.plot([], [], label="Phanerozoic", color='#7E2954', linewidth=3)
plt.legend(loc='best', fontsize=14, frameon=False)

# plt.axvline(x=single_CN_bond_wavelength, color='#0072B2', linewidth=1.5)
# plt.text(single_CN_bond_wavelength, 10 ** (-1.5), "C-N", rotation=90, verticalalignment='center',
#          horizontalalignment='right', color='black', fontsize=20)
# plt.axvline(x=single_CN_bond_wavelength, color='#0072B2', linewidth=10, alpha=0.5)
#
# plt.axvline(x=single_CC_bond_wavelength, color='red', linewidth=1.5)
# plt.text(single_CC_bond_wavelength, 10 ** (-3.5), "C-C", rotation=90, verticalalignment='center',
#          horizontalalignment='right', color='black', fontsize=20)
# plt.axvline(x=single_CC_bond_wavelength, color='red', linewidth=10, alpha=0.5)
#
# plt.axvline(x=single_PO_bond_wavelength, color='black', linewidth=1.5)
# plt.text(single_PO_bond_wavelength, 10 ** (-4.5), "P-O", rotation=90, verticalalignment='center',
#          horizontalalignment='left', color='black', fontsize=20)
#
# plt.axvline(x=single_CH_bond_wavelength, color='black', linewidth=1.5)
# plt.text(single_CH_bond_wavelength, 10 ** (-1.5), "C-H", rotation=90, verticalalignment='center',
#          horizontalalignment='left', color='black', fontsize=20)
#
# plt.axvline(x=single_CO_bond_wavelength, color='black', linewidth=1.5)
# plt.text(single_CO_bond_wavelength, 10 ** (-2.5), "C-O", rotation=90, verticalalignment='center',
#          horizontalalignment='right', color='black', fontsize=20)
#
# plt.axvline(x=double_CC_bond_wavelength, color='black', linewidth=1.5)
# plt.text(double_CC_bond_wavelength, 10 ** (-0.5), "C=C", rotation=90, verticalalignment='center',
#          horizontalalignment='left', color='black', fontsize=20)
#
# plt.axvline(x=double_CO_bond_wavelength, color='black', linewidth=1.5)
# plt.text(double_CO_bond_wavelength, 10 ** (-0.5), "C=O", rotation=90, verticalalignment='center',
#          horizontalalignment='left', color='black', fontsize=20)

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
file_path_save_pdf = os.path.join(directory, f"Irradiance_vs_Wavelength_Top_of_Atmosphere_{timestamp}.pdf")
plt.savefig(file_path_save_pdf, format='pdf')
file_path_save_png = os.path.join(directory, f"Irradiance_vs_Wavelength_Top_of_Atmosphere_{timestamp}.png")
plt.savefig(file_path_save_png, format='png', dpi=1200)
plt.show()

end_time = time.time()
print(f"Elapsed time: {end_time - start_time} seconds")
