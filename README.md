# 2020-Universal-PIP
Files for universal, four-exponential positron implantation profile (PIP) paper
1. GEANT4: GEANT4 simulation input files
The GEANT4 code outputs histograms (1 um bins) of the annihilation positions of all input positrons in x/y/z along, labeled by the integer number of the material being run. Within the main file (PositronImplantationProfile.cc), one specifies the number of materials to be run in an outer for loop (the DetectorConstruction object gets passed the current material integer). The DetectorConstruction has the list of all materials, one of which is run. The current working version has two lists, one for the running of all elements at 4g/cc and another for running the 270 material database. These are selected in the DetectorConstruction.cc file. Note, to input the 22Na spectrum for the incident positrons, the file "Na22_EnergyHist.txt" must be contained within the Build directory.

2. MATLAB: MATLAB fitting files
These files take the direct output of the GEANT histogram and construct an exponential fit to the data of the number of terms specified in the file name. The coefficient values and the fit goodness parameter files are written as output.
    
All other folders contain material information, GEANT4 simulation output files, MATLAB exponential fit coefficients, and MATLAB exponential model goodness values.
