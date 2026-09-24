import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from scipy.interpolate import interp1d
import time
import os
import unicodeit
from datetime import datetime


plt.rcParams["font.family"] = ["Arial", "DejaVu Sans"]


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

target_age_1 = 2
target_age_2 = 3.5
target_age_3 = 4.603
start_time = time.time()
# endregion

file_path_00 = r'Claire_Data\SunModern.txt'
data_00 = pd.read_csv(file_path_00, sep='\\s+', comment='#', names=['wavelength', 'flux'],
                      skiprows=103, nrows=5726)
wavelengths = np.array(data_00['wavelength'])
wavelengths_conv = 10 ** (-9) * np.array(data_00['wavelength'])

only_uva = [315, 400]
until_uva = [100, 400]
only_uvb = [280, 315]
until_uvb = [100, 315]
only_uvc = [100, 280]
until_uvc = [100, 280]

file_path_water_1 = r'Water_Data\Water_Data_1.csv'
data_water_1 = pd.read_csv(file_path_water_1, comment='#', names=['wavelength', 'coefficient'], nrows=157)
water_wavelengths_1 = data_water_1['wavelength'].astype(float)
water_coefficients_1 = data_water_1['coefficient'].astype(float)
water_interpolator_1 = interp1d(water_wavelengths_1, water_coefficients_1, kind='linear', bounds_error=False,
                                fill_value="extrapolate")
water_coefficients_lst_1 = np.array(water_interpolator_1(wavelengths))

directory = r'Results'

if not os.path.exists(directory):
    os.makedirs(directory)

plt.figure(figsize=(10, 9))
plt.plot(wavelengths[(200 <= wavelengths)], water_coefficients_lst_1[(200 <= wavelengths)], color='navy', linewidth=3,
         label=f"{unicodeit.replace("H_2O")} absorption coefficients")
plt.xlabel("Wavelength (nm)", fontsize=17)
plt.ylabel(f"Attenuation in {unicodeit.replace("H_2O")} ({unicodeit.replace("m^{-1}")})", fontsize=17)
plt.tick_params(axis='both', labelsize=16)
# plt.legend(loc='best', fontsize=14, frameon=False)
plt.ylim(-0.1)

# plt.axvline(x=single_CN_bond_wavelength, color='#0072B2', linewidth=1.5)
# plt.text(single_CN_bond_wavelength, 2, "C-N", rotation=90, verticalalignment='center',
#          horizontalalignment='right', color='black', fontsize=27)
# plt.axvline(x=single_CN_bond_wavelength, color='#0072B2', linewidth=10, alpha=0.5)
#
# plt.axvline(x=single_CC_bond_wavelength, color='red', linewidth=1.5)
# plt.text(single_CC_bond_wavelength, 1.5, "C-C", rotation=90, verticalalignment='center',
#          horizontalalignment='right', color='black', fontsize=27)
# plt.axvline(x=single_CC_bond_wavelength, color='red', linewidth=10, alpha=0.5)
#
# plt.axvline(x=single_PO_bond_wavelength, color='black', linewidth=1.5)
# plt.text(single_PO_bond_wavelength, 1, "P-O", rotation=90, verticalalignment='center',
#          horizontalalignment='left', color='black', fontsize=27)
#
# plt.axvline(x=single_CH_bond_wavelength, color='black', linewidth=1.5)
# plt.text(single_CH_bond_wavelength, 2, "C-H", rotation=90, verticalalignment='center',
#          horizontalalignment='left', color='black', fontsize=27)
#
# plt.axvline(x=single_CO_bond_wavelength, color='black', linewidth=1.5)
# plt.text(single_CO_bond_wavelength, 2, "C-O", rotation=90, verticalalignment='center',
#          horizontalalignment='right', color='black', fontsize=27)

# plt.axvspan(200, only_uvc[1], label=f"UVC (100nm-280nm)", color='#7F66FF', alpha=0.35)
# plt.axvspan(only_uvb[0], only_uvb[1], label=f"UVB (280nm-315nm)", color='#E3006A', alpha=0.35)
# plt.axvspan(only_uva[0], only_uva[1], label=f"UVA (315nm-400nm), modern O3", color='#FFA500', alpha=0.35)

ax = plt.gca()

y_bottom = -0.1
y_top = 0.03
height = y_top - y_bottom

uvc_patch = patches.Rectangle((200, y_bottom), 80, height, edgecolor="none", facecolor="dimgrey", alpha=0.6)
ax.add_patch(uvc_patch)
ax.text(240, -0.05, "UVC", ha='center', va='center', fontsize=11, fontweight='bold')

uvb_patch = patches.Rectangle((280, y_bottom), 35, height, edgecolor="none", facecolor="silver", alpha=0.8)
ax.add_patch(uvb_patch)
ax.text(297.5, -0.05, "UVB", ha='center', va='center', fontsize=11, fontweight='bold')

uva_patch = patches.Rectangle((315, y_bottom), 85, height, edgecolor="none", facecolor="gainsboro", alpha=0.6)
ax.add_patch(uva_patch)
ax.text(357.5, -0.05, "UVA", ha='center', va='center', fontsize=11, fontweight='bold')

plt.tight_layout()
timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
file_path_save = os.path.join(directory, f"Water_Coefficients_vs_Wavelength_{timestamp}.pdf")
plt.savefig(file_path_save, format='pdf')
file_path_save = os.path.join(directory, f"Water_Coefficients_vs_Wavelength_{timestamp}.png")
plt.savefig(file_path_save, format='png', dpi=1200)
plt.show()

end_time = time.time()
print(f"Elapsed time: {end_time - start_time} seconds")
# endregion
