# Marine mammal recovery is associated with the resurgence of a nematode parasite
---

Dataset includes anisakid (*Contracaecum* and *Anisakis* spp.) data collected from natural history specimens of five species of fish hosts collected from 1880-2018 in Puget Sound, Washington. There is also harbor seal abundance data from Chasco et al. 2017, pollutant data from Brandenberger et al. 2008, temperature data from Race Rocks lighthouse (48.2980°N, 123.5315°W) British Columbia Lightstation Sea-Surface Temperature and Salinity Data, and bird abundance data from the Audubon Christmas Day Bird Count. 


## Description of the Data and file structure

Each fish species' anisakid count and metadata, the temperature data, the harbor seal abundance data, the bird abundance data, and the pollutant data are all housed in independent .csv files. The fish data have been compiled in the compiled_data.RDS file. The code to load the data and compile the fish, seal, temperature, bird, and pollutant data into one dataset is in load_format_data.R. 

## Sharing/access Information

Links to other publicly accessible locations of the data: https://github.com/wood-lab/Mastick_et_al_Ecology

Was data derived from another source?
If yes, list source(s): 
Wood CL, Welicky RL, Preisser WC, et al. 2023 A reconstruction of parasite burden reveals one century of climate-associated parasite decline. Proc Nat Acad Sci 120: e2211903120.

Brandenberger JM, Crecelius EA, Louchouarn P. 2008 Historical inputs and natural recovery rates for heavy metals and organic biomarkers in Puget Sound during the 20th century. Environ Sci Technol 42: 6786–90.

Chasco B, Kaplan IC, Thomas A et al. 2017 Estimates of Chinook salmon consumption in Washington State inland waters by four marine mammal predators from 1970 to 2015. Can J Fish Aquat Sci 74: 1173–94.

British Columbia Lightstation Sea-Surface Temperature and Salinity Data (Pacific), 1914-present. https://open.canada.ca/data/en/dataset/719955f2-bf8e-44f7-bc26-6bd623e82884.

National Audubon Society (2020). The Christmas Bird Count Historical Results [Online]. Available http://www.christmasbirdcount.org [October 10, 2022]

Strom A, Francis RC, Mantua NJ, Miles EL, Peterson DL. 2004 North Pacific climate recorded in growth rings of geoduck clams: a new tool for paleoenvironmental reconstruction. Geophysical Research Letters, 31(6).
