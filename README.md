# Urban-Structure-LST-Grid-sensitivity-
# STURLA Grid-Based Urban Landscape Analysis

This repository contains data and code for:

**Kremer, P., Weaver, D., & Stewart, J.D. (in review). Sensitivity of Urban Structure-Temperature Relationships to Grid Parameterization. Submitted to Ecological Informatics.**

## Repository Contents

### Code

- **Grid_STURLA_ST_for_publication_25Oct13.R** - Grid generation and analysis script that creates grids of varying sizes and orientations to analyze urban landscape structure and surface temperature relationships.

- **STURLA_Sensitivity_Analysis_for_publication_25Oct13.R** - Statistical analysis and figure generation script. Produces all manuscript figures (Figures 2-8) and supplementary materials.

### Data

- **sturla_organized_counts.csv** - Count of grid cells for each STURLA class at different resolutions and angles.

- **sturla_organized_ranking.csv** - Rank (1 = most abundant) of each STURLA class at different resolutions and angles.

- **sturla_organized_temps.csv** - Mean surface temperature (°C) for each STURLA class at different resolutions and angles.

## STURLA Classification System

The **S**patial **T**ypology of **U**rban **L**andscape **A**rrangements (STURLA) uses the following abbreviations:

- **T** = Tree canopy
- **G** = Grass/shrubs
- **B** = Bare soil
- **W** = Water
- **P** = Paved surfaces
- **L** = Low-rise buildings (1-3 stories)
- **M** = Mid-rise buildings (4-9 stories)
- **H** = High-rise buildings (>9 stories)

Classes are combinations of these elements (e.g., "TGPL" = Trees, Grass, Paved, Low-rise).

## Requirements

### R Packages

```r
# Required packages
library(terra)       # Raster data processing
library(sf)          # Spatial operations
library(tidyverse)   # Data manipulation and visualization
library(ggplot2)     # Plotting
library(viridis)     # Color palettes
library(gridExtra)   # Multiple plot arrangements
```

## Input Data Format (for Grid_STURLA_ST script)

Land cover/building height raster:
- Format: GeoTIFF
- Data type: Integer with values 1-10
- Required classes:
  - 1 = Background/No data
  - 2 = Trees (T)
  - 3 = Grass/Shrubs (G)
  - 4 = Bare Soil (B)
  - 5 = Water (W)
  - 6 = Paved
  - 7 = Low-rise Buildings (L)
  - 8 = Mid-rise Buildings (M)
  - 9 = High-rise Buildings (H)
- Coordinate system: Projected (not geographic)

## Surface temperature raster:
- Format: GeoTIFF
- Units: Kelvin × 10 (standard Landsat format)
- Same or compatible coordinate system as land cover

## Usage

### Running the Sensitivity Analysis

The sensitivity analysis script uses the provided CSV files:

```r
# Set working directory to repository location
setwd("path/to/repository")

# Source the script
source("STURLA_Sensitivity_Analysis_for_publication_25Oct13.R")
```

This will generate all manuscript figures.

### Running the Grid Generation (optional)

If you have your own input rasters:

```r
# Edit lines 131-133 in Grid_STURLA_ST_for_publication_25Oct13.R:
lulc_file <- "path/to/your/landcover.tif"
temp_file <- "path/to/your/temperature.tif"
output_file <- "path/to/output.csv"

# Source the script
source("Grid_STURLA_ST_for_publication_25Oct13.R")
```

## Figures Generated

The sensitivity analysis script produces:

- **Figure 2**: Number of unique classes by resolution and angle
- **Figure 3**: Ranking of top 5 STURLA classes (heatmaps)
- **Figure 4**: Rank dynamics of most abundant classes
- **Figure 5**: Class composition and diversity (stacked bar charts)
- **Figure 6**: Class stability analysis (3 panels)
- **Figure 7**: Mean surface temperature by STURLA class (heatmaps)
- **Figure 8**: Temperature prediction performance (R² analysis)
- **Figure S1**: Cumulative proportions (supplementary material)

## Citation

If you use this code or data, please cite:

Kremer, P., Weaver, D., & Stewart, J.D. (in review). Sensitivity of Urban Structure-Temperature Relationships to Grid Parameterization. Submitted to Ecological Informatics.

## Authors

- Peleg Kremer
- Dennis Weaver
- Justin D. Stewart

## Contact

For questions about this code or data, please contact Dr. Peleg Kremer peleg.kremer@villanova.edu

## License

This code and data are made available for research and educational purposes. Please cite the paper if you use these materials in your work.
