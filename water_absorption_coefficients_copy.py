import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import time
import os
import unicodeit
from datetime import datetime


# region Intro Print, Input Statements, Time Start
print("Note: UVC is 100nm to 280nm; UVB is 280nm to 315nm; UVA is 315nm to 400nm.")
print("Disclaimer: For photons units, numbers less than 1000 will be considered negligible (rounded to 0).")
print("Disclaimer: The bond energy of the C-N bond is for an amide (peptide) bond, not a normal C-N bond.")
print("Note: Taking average mass of 1 amino acid as 110 Da.")
print("Note: Taking an average of 300 amino acids per each chain.")

# double_CO_bond_wavelength = 161.0  # in nm
# double_CC_bond_wavelength = 199.0  # in nm
single_CP_bond_wavelength = 338.0  # in nm  #OLD DATA: 233.0nm
single_CH_bond_wavelength = 291.0  # in nm
single_CO_bond_wavelength = 334.0  # in nm
single_CC_bond_wavelength = 346.0  # in nm
single_PO_bond_wavelength = 266.0  # in nm   #OLD DATA: 357.0nm
single_CN_bond_wavelength = 392.0  # in nm

target_age_1 = 2
target_age_2 = 3.5
target_age_3 = 4.603
start_time = time.time()
# endregion

file_path_00 = r'C:\Users\User\Documents\McGill Internship 2024\Claire_Data\SunModern.txt'
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

file_path_water_1 = r'C:\Users\User\Documents\McGill Internship 2024\Water_Data\Water_Data_1.csv'
data_water_1 = pd.read_csv(file_path_water_1, comment='#', names=['wavelength', 'coefficient'], nrows=157)
water_wavelengths_1 = data_water_1['wavelength'].astype(float)
water_coefficients_1 = data_water_1['coefficient'].astype(float)
water_interpolator_1 = interp1d(water_wavelengths_1, water_coefficients_1, kind='linear', bounds_error=False,
                                fill_value="extrapolate")
water_coefficients_lst_1 = np.array(water_interpolator_1(wavelengths))

directory = r'C:\Users\User\Documents\McGill Internship 2024\PDF_Plots'

if not os.path.exists(directory):
    os.makedirs(directory)

plt.figure(figsize=(10, 9))
plt.plot(wavelengths[(200 <= wavelengths)], water_coefficients_lst_1[(200 <= wavelengths)], color='navy', linewidth=5)
plt.xlabel("Wavelength, nm", fontsize=40)
plt.ylabel(f"Attenuation in {unicodeit.replace("H_2O")}, {unicodeit.replace("m^{-1}")}", fontsize=40)
plt.tick_params(axis='both', labelsize=30)

plt.axvline(x=single_CN_bond_wavelength, color='black', linewidth=1.5)
plt.text(single_CN_bond_wavelength, 2, "C-N", rotation=90, verticalalignment='center',
         horizontalalignment='right', color='black', fontsize=35)

plt.axvline(x=single_CC_bond_wavelength, color='black', linewidth=1.5)
plt.text(single_CC_bond_wavelength, 1.5, "C-C", rotation=90, verticalalignment='center',
         horizontalalignment='left', color='black', fontsize=35)

plt.axvline(x=single_PO_bond_wavelength, color='black', linewidth=1.5)
plt.text(single_PO_bond_wavelength, 1, "P-O", rotation=90, verticalalignment='center',
         horizontalalignment='left', color='black', fontsize=35)

plt.axvline(x=single_CH_bond_wavelength, color='black', linewidth=1.5)
plt.text(single_CH_bond_wavelength, 2, "C-H", rotation=90, verticalalignment='center',
         horizontalalignment='left', color='black', fontsize=35)

plt.axvline(x=single_CO_bond_wavelength, color='black', linewidth=1.5)
plt.text(single_CO_bond_wavelength, 2, "C-O", rotation=90, verticalalignment='center',
         horizontalalignment='right', color='black', fontsize=35)

plt.axvline(x=single_CP_bond_wavelength, color='black', linewidth=1.5)
plt.text(single_CO_bond_wavelength, 1, "C-P", rotation=90, verticalalignment='center',
         horizontalalignment='left', color='black', fontsize=35)

plt.axvspan(200, only_uvc[1], label=f"UVC (100nm-280nm)", color='#7F66FF', alpha=0.25)
plt.axvspan(only_uvb[0], only_uvb[1], label=f"UVB (280nm-315nm)", color='#E3006A', alpha=0.25)
plt.axvspan(only_uva[0], only_uva[1], label=f"UVA (315nm-400nm), modern O3", color='#FF6600', alpha=0.25)

plt.tight_layout()
timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
file_path_save = os.path.join(directory, f"Water_Coefficients_vs_Wavelength_{timestamp}.pdf")
plt.savefig(file_path_save, format='pdf')
plt.show()

end_time = time.time()
print(f"Elapsed time: {end_time - start_time} seconds")
# endregion
