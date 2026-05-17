import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import time
import os
from datetime import datetime
import unicodeit


plt.rcParams["font.family"] = ["Arial", "DejaVu Sans"]


start_time = time.time()

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
O3_sigma_interpolator = interp1d(data_O3_sigma_O2, data_O3_sigma_sigma, kind='linear', bounds_error=False,
                                 # TODO: See if it's ok to have linear
                                 fill_value="extrapolate")
# endregion

# region Time functions for gases in atmosphere
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


def O2_pressure(t_age, O2_history):
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


def O3_presence(t_age, O2_history):
    if "l" in O2_history:
        data_O2_time_time = data_O2_time_time_lower
        data_O2_time_pressure = data_O2_time_pressure_lower
    elif "u" in O2_history:
        data_O2_time_time = data_O2_time_time_upper
        data_O2_time_pressure = data_O2_time_pressure_upper
    else:
        data_O2_time_time = data_O2_time_time_average
        data_O2_time_pressure = data_O2_time_pressure_average
    if t_age < 2.203:
        return 0
    else:
        if t_age >= 4.603:
            t_age = 4.603
        O2_time_interpolator = interp1d(data_O2_time_time, data_O2_time_pressure, kind='linear', bounds_error=False,
                                        fill_value='extrapolate')
        # print(f"O2 compared from {t_age} to today: {O2_time_interpolator(4.603 - t_age)}")
        # print(f"O3 layer: {O3_sigma_interpolator(O2_time_interpolator(4.603 - t_age))}")
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

# region Making lists for plotting
time_steps = 0.6  # in Gyr
time_list = [time_steps]
time_change = 0.001  # in Gyr
max_time = 4.603

N2_pressure_lst = [N2_pressure(time_steps)]
CO2_pressure_lst = [CO2_pressure(time_steps)]
O2_pressure_lower_lst = [O2_pressure(time_steps, "lower")]
O2_pressure_average_lst = [O2_pressure(time_steps, "average")]
O2_pressure_upper_lst = [O2_pressure(time_steps, "upper")]
O3_presence_lower_lst = [O3_presence(time_steps, "lower")]
O3_presence_average_lst = [O3_presence(time_steps, "average")]
O3_presence_upper_lst = [O3_presence(time_steps, "upper")]
total_pressure_lst = [total_pressure(time_steps)]

while time_steps < max_time:
    time_steps += time_change
    time_list.append(time_steps)
    N2_pressure_lst.append(N2_pressure(time_steps))
    CO2_pressure_lst.append(CO2_pressure(time_steps))
    O2_pressure_lower_lst.append(O2_pressure(time_steps, "lower"))
    O2_pressure_average_lst.append(O2_pressure(time_steps, "average"))
    O2_pressure_upper_lst.append(O2_pressure(time_steps, "upper"))
    O3_presence_lower_lst.append(O3_presence(time_steps, "lower"))
    O3_presence_average_lst.append(O3_presence(time_steps, "average"))
    O3_presence_upper_lst.append(O3_presence(time_steps, "upper"))
    total_pressure_lst.append(total_pressure(time_steps))
# endregion

# region Directory making
directory = r'C:\Users\User\Documents\McGill Internship 2024\PDF_Plots'
if not os.path.exists(directory):
    os.makedirs(directory)
# endregion

# region Plotting itself
fig, ax1 = plt.subplots(1, 1, figsize=(10, 5))
ax1.plot([4.603-x for x in time_list], N2_pressure_lst, color='crimson', linewidth=3, label=unicodeit.replace("N_2"))
ax1.plot([4.603-x for x in time_list], N2_pressure_lst, color='crimson', linewidth=8, alpha=0.5)
ax1.plot([4.603-x for x in time_list], CO2_pressure_lst, color='darkgreen', linewidth=3, label=unicodeit.replace("CO_2"))
ax1.plot([4.603-x for x in time_list], CO2_pressure_lst, color='darkgreen', linewidth=8, alpha=0.5)
ax1.plot([4.603-x for x in time_list], O2_pressure_lower_lst,color='#28CC7A', linestyle='dotted', linewidth=3)
ax1.plot([4.603-x for x in time_list], O2_pressure_average_lst,color='#28CC7A', linewidth=3, label=unicodeit.replace("O_2"))
ax1.plot([4.603-x for x in time_list], O2_pressure_average_lst,color='#28CC7A', linewidth=8, alpha=0.5)
ax1.plot([4.603-x for x in time_list], O2_pressure_upper_lst,color='#28CC7A', linestyle='dashed', linewidth=3)
ax1.plot([4.603-x for x in time_list], total_pressure_lst, color="#4D0000", linewidth=3)
ax1.plot([4.603-x for x in time_list], total_pressure_lst, color="#4D0000", linewidth=8, alpha=0.5)
ax1.set_xlabel("Time (Ga)", fontsize=17)
ax1.set_ylabel("Partial Pressure (bar)", fontsize=17)
ax1.tick_params(axis='both', labelsize=16)
ax2 = ax1.twinx()
ax2.plot([4.603-x for x in time_list], O3_presence_lower_lst, color='#3333FF', linestyle='dotted', linewidth=3)
ax2.plot([4.603-x for x in time_list], O3_presence_average_lst, label=unicodeit.replace("O_3"), color='#3333FF', linewidth=3)
ax2.plot([4.603-x for x in time_list], O3_presence_average_lst, color='#3333FF', linewidth=8, alpha=0.5)
ax2.plot([4.603-x for x in time_list], O3_presence_upper_lst, color='#3333FF', linestyle='dashed', linewidth=3)
ax2.set_ylabel(f"{unicodeit.replace("O_{3}")} ({unicodeit.replace("kg/m^{2}")})", fontsize=17, labelpad=8)
ax2.plot([], [], label=unicodeit.replace("N_2"), color='crimson', linewidth=3)
ax2.plot([], [], label=unicodeit.replace("CO_2"), color='darkgreen', linewidth=3)
ax2.plot([], [], label=unicodeit.replace("O_2"), color='#28CC7A', linewidth=3)
ax2.plot([], [], label=f"Min {unicodeit.replace("O_2")} pressure values", color='black', linewidth=3, linestyle="dotted")
ax2.plot([], [], label=f"Max {unicodeit.replace("O_2")} pressure values", color='black', linewidth=3, linestyle='dashed')
ax2.plot([], [], label=f"Average {unicodeit.replace("O_2")} pressure values", color='black', linewidth=3)
ax2.tick_params(axis='y', labelsize=17)
ax2.set_ylim(0, 0.007)
ax2.legend(loc='best', fontsize=14, frameon=False)

# plt.axvline(x=4.603-2, color='#D55E00', linewidth=5)
# plt.axvline(x=4.603-3.5, color='#009E73', linewidth=5)
# plt.axvline(x=4.603-4.603, color='#7E2954', linewidth=5)

# plt.axvspan(4.603-0, 4.603-2.203, color='#D55E00', alpha=0.35)
# plt.axvspan(4.603-2.203, 4.603-4.003, color='#009E73', alpha=0.35)
# plt.axvspan(4.603-4.003, -0.1, color='#7E2954', alpha=0.35)

plt.gca().set_xlim(abs(min(time_list)-4.603), -0.1)

plt.tight_layout()
timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
file_path_save = os.path.join(directory, f"Atmosphere_Gases_vs_Time_{timestamp}.pdf")
plt.savefig(file_path_save, format='pdf')
file_path_save = os.path.join(directory, f"Atmosphere_Gases_vs_Time_{timestamp}.png")
plt.savefig(file_path_save, format='png', dpi=1200)
plt.show()
# endregion

end_time = time.time()
print(f"Elapsed time: {end_time - start_time} seconds")
