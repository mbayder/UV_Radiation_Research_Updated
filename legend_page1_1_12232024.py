import matplotlib.pyplot as plt
import os
import unicodeit
from datetime import datetime


plt.rcParams["font.family"] = ["Arial", "DejaVu Sans"]


directory = r'Results'
if not os.path.exists(directory):
    os.makedirs(directory)

plt.figure()
plt.plot([], [], label="Before GOE (2.6 Ga)", color='#D55E00', linewidth=3)
plt.plot([], [], label="Before NOE, after GOE (1.1 Ga)", color='#009E73', linewidth=3)
plt.plot([], [], label="After NOE (Today)", color='#7E2954', linewidth=3)
plt.plot([], [], label="C-N bond", color='#0072B2', linewidth=3)
plt.plot([], [], label="C-C bond", color='red', linewidth=3)
plt.plot([], [], label="Other bonds", color='black', linewidth=3)
plt.plot([], [], label="UVA", color='#FFA500', linewidth=3)
plt.plot([], [], label="UVB", color='#E3006A', linewidth=3)
plt.plot([], [], label="UVC", color='#7F66FF', linewidth=3)
plt.legend(loc='center', fontsize=18, frameon=False)
plt.axis('off')
timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
file_path_save = os.path.join(directory, f"legend_page1_part1_{timestamp}.pdf")
plt.savefig(file_path_save, format='pdf')
plt.show()

plt.figure()
plt.plot([], [], label="Total atmosphere", color='#4D0000', linewidth=3)
plt.plot([], [], label=unicodeit.replace("N_2"), color='crimson', linewidth=3)
plt.plot([], [], label=unicodeit.replace("CO_2"), color='darkgreen', linewidth=3)
plt.plot([], [], label=unicodeit.replace("O_2"), color='#28CC7A', linewidth=3)
plt.plot([], [], label=unicodeit.replace("O_3"), color='#3333FF', linewidth=3)
plt.plot([], [], label=f"Min {unicodeit.replace("O_2")} pressure values", color='black', linewidth=3, linestyle="dotted")
plt.plot([], [], label=f"Max {unicodeit.replace("O_2")} pressure values", color='black', linewidth=3, linestyle='dashed')
plt.plot([], [], label=f"Average {unicodeit.replace("O_2")} pressure values", color='black', linewidth=3)
plt.legend(loc='center', fontsize=18, frameon=False)
plt.axis('off')
timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
file_path_save = os.path.join(directory, f"legend_page1_part2_{timestamp}.pdf")
plt.savefig(file_path_save, format='pdf')
plt.show()

# plt.figure()
# plt.plot([], [], label="Before GOE (2.6 Ga)", color='#D55E00', linewidth=3)
# plt.plot([], [], label="Before NOE, after GOE (1.1 Ga)", color='#009E73', linewidth=3)
# plt.plot([], [], label="After NOE (Today)", color='#7E2954', linewidth=3)
# plt.plot([], [], label="C-N bond", color='#0072B2', linewidth=3)
# plt.plot([], [], label="C-C bond", color='red', linewidth=3)
# plt.plot([], [], label="Other bonds", color='black', linewidth=3)
# plt.plot([], [], label="UVA", color='#FFA500', linewidth=3)
# plt.plot([], [], label="UVB", color='#E3006A', linewidth=3)
# plt.plot([], [], label="UVC", color='#7F66FF', linewidth=3)
# plt.plot([], [], label="Total atmosphere", color='#4D0000', linewidth=3)
# plt.plot([], [], label=unicodeit.replace("N_2"), color='crimson', linewidth=3)
# plt.plot([], [], label=unicodeit.replace("CO_2"), color='darkgreen', linewidth=3)
# plt.plot([], [], label=unicodeit.replace("O_2"), color='#28CC7A', linewidth=3)
# plt.plot([], [], label=unicodeit.replace("O_3"), color='#3333FF', linewidth=3)
# plt.plot([], [], label=f"Min {unicodeit.replace("O_2")} pressure values", color='black', linewidth=3, linestyle="dotted")
# plt.plot([], [], label=f"Max {unicodeit.replace("O_2")} pressure values", color='black', linewidth=3, linestyle='dashed')
# plt.plot([], [], label=f"Average {unicodeit.replace("O_2")} pressure values", color='black', linewidth=3)
# plt.legend(loc='center', fontsize=16, frameon=False)
# plt.axis('off')
# timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
# file_path_save = os.path.join(directory, f"legend_page1_full_{timestamp}.pdf")
# plt.savefig(file_path_save, format='pdf')
# plt.show()
