# Marine mammal recovery is associated with the resurgence of a nematode parasite
---

Dataset includes anisakid (*Contracaecum* and *Anisakis* spp.) data collected from natural history specimens of five species of fish hosts collected from 1880-2018 in Puget Sound, Washington. These files include Pacific hake (hake_2022.02.21.csv), Pacific herring (herring_2022.02.23.csv), walleye pollock (pollock_2022.02.21.csv), copper rockfish (rockfish_2022.02.21.csv), and surf smelt (smelt_2022.02.21.csv). There is also harbor seal abundance data from Chasco et al. 2017 (harborSeals.csv), pollutant data from Brandenberger et al. 2008 (brandenberger_pollutants.csv), and temperature data from Race Rocks lighthouse (48.2980°N, 123.5315°W) British Columbia Lightstation Sea-Surface Temperature and Salinity Data (RR_temp.csv).


## Description of the Data and file structure

Each fish species' anisakid count and metadata, the temperature data, the harbor seal abundance data, and the pollutant data are all housed in independent .csv files in the data folder. For analysis, the fish data (hake, herring, pollock, rockfish, and smelt) have been compiled into one dataset in the compiled_data.RDS file. The code to load the data and compile the fish, seal, temperature, bird, and pollutant data into one dataset is in load_format_data.R. 

# RR_temp.csv
YEAR: year of temperature reading
temp_na_rm: average temperature, calculated with NAs removed

# brandenberger_pollutants.csv
year: year of pollutant reading from Tacoma, WA and Seattle, WA
Pb:	average lead reading for that year in μg/g
As:	average arsenic reading for that year in μg/g
Zn:	average zinc reading for that year in μg/g
Ni:	average nickel reading for that year in μg/g
V: average vanadium reading for that year in μg/g
Cr:	average chromium reading for that year in μg/g
Cu:	average copper reading for that year in μg/g
Ba:	average barium reading for that year in μg/g
Be:	average beryllium reading for that year in μg/g
Sig8_lignin:	average lignin reading for that year in mg/g
Lamb8: average soil biomarker concentration reading for that year
Bd.V_soil_biomarker: average soil biomarker concentration reading for that year

# harborSeals.csv



## Sharing/access Information

Links to other publicly accessible locations of the data: https://github.com/wood-lab/Mastick_et_al_Ecology

Data derived from the following sources:
Wood CL, Welicky RL, Preisser WC, et al. 2023 A reconstruction of parasite burden reveals one century of climate-associated parasite decline. Proc Nat Acad Sci 120: e2211903120.

Brandenberger JM, Crecelius EA, Louchouarn P. 2008 Historical inputs and natural recovery rates for heavy metals and organic biomarkers in Puget Sound during the 20th century. Environ Sci Technol 42: 6786–90.

Chasco B, Kaplan IC, Thomas A et al. 2017 Estimates of Chinook salmon consumption in Washington State inland waters by four marine mammal predators from 1970 to 2015. Can J Fish Aquat Sci 74: 1173–94.

British Columbia Lightstation Sea-Surface Temperature and Salinity Data (Pacific), 1914-present. https://open.canada.ca/data/en/dataset/719955f2-bf8e-44f7-bc26-6bd623e82884.

National Audubon Society (2020). The Christmas Bird Count Historical Results [Online]. Available http://www.christmasbirdcount.org [October 10, 2022]

Strom A, Francis RC, Mantua NJ, Miles EL, Peterson DL. 2004 North Pacific climate recorded in growth rings of geoduck clams: a new tool for paleoenvironmental reconstruction. Geophysical Research Letters, 31(6).

## Code/software
Scripts were composed in R, version 4.2.3. For analysis, the fish data (hake, herring, pollock, rockfish, and smelt) have been compiled into one dataset in the compiled_data.RDS file. The code to load the data and compile the fish, seal, temperature, bird, and pollutant data into one dataset is in load_format_data.R. 
