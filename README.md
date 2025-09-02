# Marine mammal recovery is associated with the resurgence of a nematode parasite
---

Dataset includes anisakid (*Contracaecum* and *Anisakis* spp.) data collected from natural history specimens of five species of fish hosts collected from 1880-2018 in Puget Sound, Washington. These files include Pacific hake (hake_2022.02.21.csv), Pacific herring (herring_2022.02.23.csv), walleye pollock (pollock_2022.02.21.csv), copper rockfish (rockfish_2022.02.21.csv), and surf smelt (smelt_2022.02.21.csv). There is also harbor seal abundance data from Chasco et al. 2017 (harborSeals.csv), pollutant data from Brandenberger et al. 2008 (brandenberger_pollutants.csv), and temperature data from Race Rocks lighthouse (48.2980°N, 123.5315°W) British Columbia Lightstation Sea-Surface Temperature and Salinity Data (RR_temp.csv).


## Description of the data and file structure

Each fish species' anisakid count and metadata, the temperature data, the harbor seal abundance data, and the pollutant data are all housed in independent .csv files in the data folder. For analysis, the fish data (hake, herring, pollock, rockfish, and smelt) have been compiled into one dataset in the compiled_data.RDS file. The code to load the data and compile the fish, seal, temperature, bird, and pollutant data into one dataset is in load_format_data.R. 

### RR_temp.csv
YEAR: year of temperature reading
temp_na_rm: average temperature, calculated with NAs removed

### brandenberger_pollutants.csv
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

### harborSeals.csv
Year: year of the estimate
HS.WA.StraitJuanDeFuca: harbor seal population estimate in the Strait of Juan de Fuca	
HS.WA.SanJuanIslands: harbor seal population estimate in the San Juan Islands, WA
HS.WA.EasternBays: harbor seal population estimate in the Eastern Bays of Washington	
HS.WA.PugetSound: harbor seal population estimate from Puget Sound, WA
HS.WA.HoodCanal: harbor seal population estimate from Hood Canal, WA	
HS.WA.Georgia.Strait: harbor seal population estimate from Georgia Strait, WA	
HS.WA.OlympicPeninsula: harbor seal population estimate from the Olympic Peninsula in Washington
HS.WA.CoastalEstuaries: harbor seal population estimate in the Coastal Estuaries of Washington	
HS.BC.Inland: harbor seal population estimate in inland British Columbia

#### hake_2022.02.21.csv , herring_2022.02.23.csv , pollock_2022.02.21.csv , rockfish_2022.02.21.csv , smelt_2022.02.21.csv
Year.Collected: the year the fish was collected from the wild and preserved in the Burke Museum's Ichthyology collection
Columns named in the format XX_XX_XX: the count of each genus of parasites found in each fish, the suffix indicating the initials of the scientific name of that species in some cases, or "all" in others, where it was found in more than one species. For example, Anisakis_sp_all indicates Anisakis sp. nematodes found throughout the body of the particular host. 
latjitt	site: the collection site latitude, but the coordinate has been "jittered" by +/- a small amount to reduce overlap in mapping and analysis.	
Fish.ID: the catalog number of each fish dissected
Standard.Length..mm.: standard length of the fish in mm
Total.Length..mm.: total length of the fish in mm
long: longitude the fish was collected
lat: latitude the fish was collected	
fish_density: an estimate of the density of the fish in Puget Sound from that year	
temp_na_rm: data compiled from the RR_temp.csv dataset, temperature with NAs removed	
sTemp: temperature standardized to a value between 0-1. 
Pb	As	Zn	Ni	V	Cr	Cu	Ba	Be	Sig8_lignin	Lamb8	Bd.V_soil_biomarker: pollutant values compiled from brandenberger_pollutants.csv

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
