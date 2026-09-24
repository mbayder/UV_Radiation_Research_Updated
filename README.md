# Light for a Living: Ultraviolet Radiation in Surface Waters Over Geological Time

### **Part of project and paper "Microbial exoenzymes promoted the transition to an oxygenated Earth" by Sanger et al.**

**Code Author:** Maria Bayder

**Code Author Email:** maria.bayder@mail.mcgill.ca

**Supervisors:** Prof. Nicolas Cowan and Prof. Nagissa Mahmoudi

**Code Use Instructions:** If all folders, files, and programs were installed properly, they should all work 
right away! Here are some helpful tips and descriptions: 
- Some code may require user input, so please be sure to be able to type them in somewhere on your screen.
- The outputs should include some printed information, data in CSV format, and figures in PNG and PDF formats. 
- If you wish, you may use your own local directories for outputs; simply replace 'Results' directory with the path to 
your local directory. 
- Please do not touch the relative paths to data files. 
- To avoid potential library conflicts, please use the following package versions: contourpy 1.3.3 ;
cycler 0.12.1 ; fonttools 4.63.0 ; kiwisolver 1.5.0 ; matplotlib 3.10.9 ; numpy 2.4.5 ; packaging 26.2 ; pandas 3.0.3 ; 
pillow 12.2.0 ; pyparsing 3.3.2 ; python-dateutil 2.9.0.post0 ; pytz 2026.2 ; scipy 1.17.1 ; six 1.17.0 ; 
tzdata 2026.2 ; unicodeit 0.7.5 

Please feel free to email the code author if you have any issues, questions, comments, suggestions, etc. 

**Scientific Background of the Code:** 

We simulate the ultraviolet (UV) flux available for photodegradation in the surface layer of Earth's ocean throughout 
its geological history. Our radiative transfer model accounts for the secular changes in the spectrum of the Sun, 
the evolving composition of the Earth's atmosphere, and the mixed layer depth of the ocean. In the following three 
paragraphs, we summarize each of the three principle modules: the Sun, Earth's atmosphere, and Earth's oceans. 
In the fourth paragraph, we explain how to obtain the flux of photodegrading UV light and the rate of organic matter 
degradation. The final paragraph outlines a sensitivity analysis of important variables.

We use the Sun's irradiance spectrum throughout geological time from Claire et al. (2012), and we interpolate in time 
using univariate spline interpolation to estimate the Sun's irradiance at multiple epochs. 

When modeling the atmosphere, we consider only carbon dioxide, nitrogen, oxygen, and ozone (Cnosssen et al., 2007). 
We assume a total pressure of 0.23 bar on 2.7 Ga Archean Earth (Som et al., 2016) with 70\% carbon dioxide 
(Lehmer et al., 2020) and 30\% nitrogen. Our model assumes that the total atmospheric pressure and percentages of carbon
dioxide and nitrogen change linearly from 2.7Ga Archean values to modern values of 1.013 bar, 0.042\% carbon dioxide, 
and 78\% nitrogen. We use the evolution of oxygen over time from Lyons et al. (2014), and we estimate the thickness of 
the ozone layer over time based on its dependency on oxygen levels, as described by Cooke et al. (2022). Using the 
calculated amounts of gases for a specific time and their absorption cross-sections from Cnossen et al. (2007), our 
model calculates irradiance spectra at Earth’s surface at different epochs. 

We model attenuation of UV light in ocean water following Hale and Querry (1973) and Cockell (2000), and we assume that 
ocean water optical properties remain constant throughout Earth's history. Our model calculates the irradiance spectrum 
as a function of the depth in the ocean. Parcels of water in the mixed layer circulate between different depths, so they
are exposed to varying levels of UV radiation. Our model therefore calculates the average irradiance spectrum in the 
ocean's mixed layer at different geological epochs. 

Afterwards, we choose a limiting wavelength range within the UV range of 100 nm to 400 nm, where the maximum corresponds
to the longest wavelength radiation that can photodegrade amino acids. Our model integrates the average irradiance 
spectra in the ocean mixed layer at different epochs over the chosen wavelength range to obtain the average UV fluxes. 
Our model then converts the UV flux to the rate of organic matter breakdown. According to Tiessen et al. (2012), there 
are an average of 300 amino acids per chain, and each amino acid has a mass of approximately 110 Da. From the UV fluxes,
we calculate the fluxes of photons with enough energy to break peptide bonds in the ocean's mixed layer for different 
geological epochs. The program converts the number of photons to the mass of broken-down organic matter. Thus, we obtain
the rate of organic matter breakdown due to photodegradation over geological time.

To assess the robustness of our results, we examined how the limiting wavelength range, depth of the ocean's mixed 
layer, and atmospheric composition during different epochs affected our results. We found that the factor that impacts 
our results the most is the chosen limiting wavelength range. If the limiting wavelength range is between 100 nm and 
315 nm, then we find that the rate of organic matter breakdown decreases significantly over billions of years. 
However, if we choose a limiting wavelength range beyond 315 nm, the rate of breakdown decreases insignificantly or even
increases over billions of years. This effect occurs because the absorption of UV light wavelengths longer than 315 nm 
is almost unaffected by the changes in atmospheric composition. When we used different models of oxygen evolution in the
atmosphere, these different evolutions did not affect the general trend of the organic matter breakdown rate over time, 
but they affected how much and when increases or decreases in the rate occur. Changing the depth of the ocean's mixed 
layer, the total atmospheric pressure, and the proportion of carbon dioxide to nitrogen affected the magnitude of the 
rate of organic matter breakdown and minimally affected how much and when the rate increases or decreases.

To see figures and bibliography, please refer to "Microbial exoenzymes promoted the transition to an oxygenated 
Earth" by Sanger et al.
