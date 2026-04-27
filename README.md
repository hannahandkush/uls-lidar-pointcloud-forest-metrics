# UAV and LiDAR Data Processing for Forest Metrics

This repository is dedicated to the processing and analysis of Unmanned Aerial Vehicle (UAV) and Light Detection and Ranging (LiDAR) data. The primary focus is on extracting high-precision forest metrics from point clouds, specifically for Eucalyptus forest management and research.

## Project Scope
- **UAV-LiDAR Integration**: Handling high-density point clouds captured from UAV platforms.
- **Forest Inventory**: Automated extraction of individual tree parameters and stand-level metrics.
- **Spatial Modeling**: Generation of high-resolution Digital Terrain Models (DTM), Canopy Height Models (CHM), and Above-Ground Biomass (AGB) maps.

## Current Projects
### [Eucalyptus Project](./Eucalyptus_Project)
A complete workflow for processing ULS LiDAR data in Eucalyptus stands, including:
- Point cloud cleaning and noise removal.
- Ground classification and normalization.
- Area-Based Approach (ABA) metrics extraction.
- Spatial mapping of forest attributes (p95, VCC, AGB).

## Tools and Technologies
- **Language**: R
- **Key Libraries**: `lidR`, `raster`, `sf`, `ggplot2`
- **Data Formats**: LAS/LAZ (Point Clouds), GeoTIFF (Rasters), SHP (Vector Data)

---
*Maintained by Hannah Nathanson*
