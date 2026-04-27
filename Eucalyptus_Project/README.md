# Eucalyptus Forest Metrics from ULS LiDAR Point Clouds

## 1. Overview
This project focuses on processing Unmanned Laser Scanning (ULS) LiDAR point clouds to extract critical forest metrics for Eucalyptus stands. The workflow covers noise removal, ground classification, normalization, and the generation of Area-Based Approach (ABA) metrics, including Canopy Height Models (CHM) and Above-Ground Biomass (AGB) maps.

## 2. Project Structure
The project is organized into functional directories following the processing pipeline:

```
.
├── 01_data/           # Raw LiDAR data (.las, .lax) and site shapefiles.
├── 02_scripts/        # R scripts for DTM/CHM generation and metrics extraction.
├── 03_processing/     # Intermediate processing stages (Noise removal, Normalization).
├── 04_documents/      # Methodology presentations and technical notes.
├── 05_outputs/        # Final tiff maps (AGB, VCC, p95) and ABA metrics CSVs.
└── G3/                # Group-specific working files and supplementary data.
```

## 3. Key Processing Steps
1.  **Noise Removal**: Cleaning raw point clouds to remove outliers and atmospheric interference.
2.  **Ground Classification**: Identifying ground points to establish a digital terrain reference.
3.  **Normalization**: Converting absolute altitudes to heights above ground level.
4.  **Metrics Extraction**: 
    *   **CHM**: Canopy Height Model at 1m resolution.
    *   **ABA Metrics**: Statistical summaries of height distributions (e.g., p95).
    *   **AGB/VCC**: Spatial mapping of Above-Ground Biomass and Vegetation Canopy Cover.

## 4. Usage
The primary analysis is conducted using the `liDR` package in R.
*   The main processing logic is located in: `02_scripts/Practical_Lab1_liDR_DEM_CHM.R`.
*   Outputs such as `ABA_metrics.tif` and `AGB_map.png` provide the final spatial insights.

## 5. Requirements
*   **Software**: R (version 4.0+)
*   **Key Libraries**: `lidR`, `raster`, `sf`, `ggplot2`.

---
*This project was developed as part of the ULS LiDAR Point Cloud Forest Metrics course.*
