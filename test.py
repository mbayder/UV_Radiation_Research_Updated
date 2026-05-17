import matplotlib.pyplot as plt
import os
from datetime import datetime

plt.plot([], [], label="Legend Entry 1", color='blue')
plt.plot([], [], label="Legend Entry 2", color='orange')
plt.legend(loc='center', fontsize=30, frameon=False)  # Position the legend in the center (adjust as needed)
plt.axis('off')  # Turn off the axis
directory = r'C:\Users\User\Documents\McGill Internship 2024\PDF_Plots'
if not os.path.exists(directory):
    os.makedirs(directory)
timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
file_path_save = os.path.join(directory, f"test_legend_{timestamp}.pdf")
plt.savefig(file_path_save, format='pdf')
plt.show()
