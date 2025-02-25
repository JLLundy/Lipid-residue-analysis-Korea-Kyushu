# Data and R scripts for the manuscript "Lipid residue analysis reveals divergent culinary practices in Japan and Korea at the dawn of intensive agriculture"

This repository contains data and scripts used in the following paper:

Craig, O.E., Lundy, J., Bondetti, M., Nicholson-Lailey, S., Murakami,N., Szuki ,M ., Stevens, C.,  Lucquin, A., Talbot, H..M, Son, J.H. Fujio, S., Sakamoto, M., Yamashita, Y., Kobayashi, K., Crema, E.R., Shoda, S. (2025), Lipid residue analysis reveals divergent culinary practices in Japan and Korea at the dawn of intensive agriculture

The repository is organised into the following directories:

* _Data_ ... contains all raw data required for analyses
* _isotopesandmap_ ... contains R scripts for generating figures pertaining human stable isotope data
* _organicresidues_ ... contains R scripts and resulting outputs for analysing organic residues from ceramics
* _c14dates_ ... contains R scripts and resulting outputs for the analyses of radiocarbon dates.
* _figures_ ... contains figure outputs used in the manuscript.

Each analysis directory (_Isotopesandmap_,_organicresidues_, and _c14dates_) is self-contained and can be executed independently.  

### Data 

The data directory contains CSV files required for all analyses and directly referenced in the manuscript (e.g. 'Dataset S1'):

| Filename        | Name       | Description                                                                      |
|-----------------|------------|----------------------------------------------------------------------------------|
| S1.csv          | Dataset S1 | Organic residue analysis datafile with details of molecular and isotopic results |
| S1_metadata.csv | Dataset S1 | Metadata for Dataset S1.                                                         |
| S2.csv          | Dataset S2 | Reference isotopic data for authentic food products.                             |
| S3.csv          | Dataset S3 | Reference fatty acid concentrations for modern food products.                    |
| S4.csv          | Dataset S4 | Stable Human Isotope Data                                                        |
| S5.csv          | Dataset S5 | 14C dated rice and millet remains                                                |
| S6.csv          | Dataset S6 | Output from mixing model                                                         |

### c14dates 
The directory is organised into a series of R script and R image files with the prefix number representing the execution order. Estimates of the arrival date of rice and millet can be computed by running first `01_prepare_arrival_analysis.R` (requiring the input data `S5.csv`), which will produce the R image files `02_millets_run.RData` and `02_rice_data.RData`. The two files are required for the execution of the core Bayesian analyses (`03_run_milletsmodel.R` and `03_run_ricemodel.R`), each producing an R image file (`04_milletsarrival_results.RData` and `04_ricearrival_results.RData`). Bayesian analyses of population growth rates estimated from the time-frequency of radiocarbon dates can be executed with the script `05_korea_growthrates.R` producing the R image file `06_koreagrowth_results.RData`. The script automatically downloads the required data from the web. Finally, the script `07_plotresults.R` combines the output of the two analyses on a single figure (`c14_analysis.pdf`), corresponding to figure 3 in tha main text. 


### isotopesandmap
The directory contains an R script file (`isotopes_and_maps.R`) that generates panels B and C of Figure 1 in the manuscript. The script requires the files `S1.csv` and `S4.csv` from the _data_ directory.

### organicresidues
The directory contains a Rmarkdown file (`organicresidues.Rmd`) required to run the Bayesian mixing model of the organic residue data. The script requires the input files `S1.csv`, ` S2.csv`, and `S3.csv` from the _data_ directory, and generates a series of files inside the directory during its execution (not contained in the repository). The final output of the model is stored in the file `S6.csv` in the _data_ directory.

### figures

# Funding
This research was funded by the ERC grant _Demography, Cultural Change, and the Diffusion of Rice and Millets during the Jomon-Yayoi transition in prehistoric Japan (ENCOUNTER)_ (Project N. 801953, PI: Enrico Crema).

# Licence
CC-BY 3.0
