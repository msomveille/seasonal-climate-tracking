This repository contains the code and data required for the analyses in the manuscript: Somveille M, Bossu CM, DeSaix MG, Alvarado AH, Gómez S, Rodríguez Otero G, Hernández-Baños BE, Smith TB & Ruegg KC. Broad-scale seasonal climate tracking is a consequence, not a driver, of avian migratory connectivity.

This repository contains three folders: scripts (which contains the code to run the analysis), resources (which contains the data), and results (where the analysis output is saved).

In the scripts folder, the file seasonal_climate_tracking_analysis.R contains the R code used for the analysis, and the file orsim.py contains the Python code used to run the ORSIM model.

The resources folder contains the file data_for_analysis_final.csv, which is a table with the following information for each individual used in the analysis: species, season of sample collection (i.e. breeding or wintering), population assignment, and geographical coordinates of sample collection). This data was obtained from previously conducted genoscape analyses (see main text and references in the manuscript), whose raw genetic data is accessible via the links provided in the Data Accessibility Statement in the manuscript. 

# Running the analysis
To run the analysis, you need to first download the ecoregion polygons (shapefile format), climate layers (raster format) and species seasonal relative abundance surfaces (raster format), from freely available sources (see manuscript), and add these files to the resources folder. The code in the file seasonal_climate_tracking_analysis.R can then be run until line 168, which then requires to run the Python script orsim.py. Once ORSIM finished running, the rest of the code in the file seasonal_climate_tracking_analysis.R can be run to output the results presented in the manuscript.

