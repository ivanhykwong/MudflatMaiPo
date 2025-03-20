# Topography of intertidal mudflats in Deep Bay (Hong Kong) from 1991 to 2021

Supplementary materials used in the following presentation:

**Kwong, I. H. Y., Lai, D. Y. F., Wong, F. K. K., & Fung, T. (2024). Topographic changes of intertidal mudflat in Deep Bay (Hong Kong) using Landsat and PlanetScope satellite imagery. Poster Presentation at American Geophysical Union (AGU) Annual Meeting 2024, 9-13 December 2024.**

- Link to AGU24 Session: https://agu.confex.com/agu/agu24/meetingapp.cgi/Paper/1680250
- ESS Open Archive: https://doi.org/10.22541/essoar.173445471.18731500/v1
- iPoster: https://agu24.ipostersessions.com/?s=75-45-52-44-23-92-B0-B7-78-D7-E9-8B-AF-D2-2C-00

---

GIS data produced from this study:

*   **MudflatElevation_DeepBayHK_yyyy-yyyy.tif**: Raster data (GeoTiff format) showing the elevation of intertidal mudflats and the coverage of other features in the study area in four periods, including (i) 1991-2000, (ii) 2001-2010, (iii) 2011-2020, and (iv) 2020-2021. Refer to the poster for the detailed methodology.
  
*   **MudflatElevation_DeepBayHK_symbology.csv**: Explanatory notes and a suggested symbology for different pixel values in the GeoTiff files.

*   **MudflatElevation_DeepBayHK_ArcGISsymbology.lyrx**: Used to apply the suggested symbology in ArcGIS Pro, as shown in the figure below.

---

R code used in this study:

*   **DeepBayHK_CreateMudflatFromLandsat.R**: Create elevation raster of intertidal mudflats from Landsat imagery for the first three periods (1991-2020)
  
*   **DeepBayHK_CreateMudflatFromPlanetScope.R**: Create elevation raster of intertidal mudflats from PlanetScope imagery for the last period (2020-2021) 

*   **DeepBayHK_MudflatAnalysis.R**: Create figures in the poster based on the elevation rasters 

---

![Fig6_Result](https://github.com/user-attachments/assets/59f62ec5-f89d-4a8f-a5ee-b037c5e86173)

*Last updated in March 2025*
