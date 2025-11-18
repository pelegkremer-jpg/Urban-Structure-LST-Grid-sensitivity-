#' ============================================================================
#' STURLA Grid-Based Analysis for Urban Landscape Classification
#' ============================================================================
#' 
#' This script generates grids of varying sizes and orientations to analyze
#' urban landscape structure and surface temperature relationships.
#' 
#' Citation:
#' Kremer, P., Weaver, D., & Stewart, J.D. (in review). 
#' Sensitivity of Urban Structure-Temperature Relationships to Grid Parameterization. 
#' Submitted to Ecological Informatics
#'
#' Authors: Peleg Kremer, Dennis Weaver, Justin D. Stewart
#' Date: October 2025
#' Version: 1.0 - Publication Release
#' 
#' ============================================================================
#' REQUIRED INPUT DATA FORMAT
#' ============================================================================
#' 
#' 1. Land cover/building height raster (lulc_file):
#'    - Format: GeoTIFF
#'    - Data type: Integer
#'    - Required classes (integer values):
#'      1 = Background/No data (excluded from analysis)
#'      2 = Trees (T)
#'      3 = Grass/Shrubs (G)
#'      4 = Bare Soil (B)
#'      5 = Water (W)
#'      6 = Roads
#'      7 = Other Paved surfaces
#'      8 = Low-rise Buildings (L: 1-3 stories)
#'      9 = Mid-rise Buildings (M: 4-9 stories)
#'      10 = High-rise Buildings (H: >9 stories)
#'    - Note: Roads (6) and Other Paved (7) are combined into Paved (P)
#'    - Coordinate system: Projected (not geographic lat/lon)
#' 
#' 2. Surface temperature raster (temp_file):
#'    - Format: GeoTIFF
#'    - Units: Kelvin * 10 (standard Landsat format)
#'      * If your data uses different units, modify conversion at line 330
#'    - Coordinate system: Same or compatible with land cover raster
#' 
#' 3. Output CSV file contains:
#'    - resolution: Grid cell size (meters)
#'    - angle: Grid orientation (degrees)
#'    - cell_id: Unique identifier for each grid cell
#'    - temperature: Mean surface temperature (°C)
#'    - pixel_count: Number of valid pixels in cell
#'    - prop_*: Proportion of each land cover type (0-1)
#'    - sturla_class: Combined classification string (e.g., "TGPL")
#'
#' NOTE FOR NYC CASE STUDY:
#' The analysis in Kremer et al. (in review) used:
#' - Land cover data from 2008 (MacFaden et al., 2012)
#' - Building heights from 2011 (NYC MapPLUTO, NYC Dept of City Planning)
#' - Landsat 7 surface temperature from summer 2008 (June-August)
#' 
#' ============================================================================

# Clear workspace
rm(list = ls())

# ============================================================================
# PACKAGE DEPENDENCIES
# ============================================================================

required_packages <- c("stringr", "exactextractr", "sp", "sf", "raster", "rgdal")

# Check if packages are installed
missing_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(missing_packages) > 0) {
  stop("Missing required packages: ", paste(missing_packages, collapse = ", "), 
       "\nPlease install using: install.packages(c('", 
       paste(missing_packages, collapse = "','"), "'))")
}

# Load libraries
library(stringr)
library(exactextractr)
library(sp)
library(sf)
library(raster)
library(rgdal)

# Print session info for reproducibility
cat("\n=== R SESSION INFORMATION ===\n")
cat("R version:", R.version.string, "\n")
cat("\nPackage versions:\n")
for(pkg in required_packages) {
  cat(sprintf("  %s: %s\n", pkg, as.character(packageVersion(pkg))))
}

#' ============================================================================
#' USER CONFIGURATION - MODIFY THESE PARAMETERS FOR YOUR DATA
#' ============================================================================

# Set base directory (modify to match your data location)
base_dir <- "."  # Current directory - change to your data folder path

# Input file paths
# USERS: Update these filenames to match your input data
lulc_file <- file.path(base_dir, "land_cover_with_buildings.tif")
temp_file <- file.path(base_dir, "surface_temperature.tif")

# Output directory
output_dir <- file.path(base_dir, "output")

# Grid parameters
# Full analysis (as used in manuscript): 21 resolutions × 19 angles = 399 combinations
resolutions <- c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 
                 200, 300, 400, 500, 600, 700, 800, 900, 
                 1000, 5000, 10000)  # in meters

angles <- seq(0, 180, by = 10)  # 0° to 180° in 10° increments

# For testing/demonstration with faster runtime, uncomment these lines:
# resolutions <- c(100, 500, 1000)
# angles <- c(0, 90, 180)

# Coordinate Reference System (CRS)
# USERS: Modify this to match your study area's projected coordinate system
# Example below is NYC State Plane (NAD83, US feet)
# Common alternatives:
#   - UTM Zone 18N (meters): '+proj=utm +zone=18 +datum=WGS84 +units=m +no_defs'
#   - State Plane meters: check spatialreference.org for your area
study_area_crs <- '+proj=lcc +lat_0=40.1666666666667 +lon_0=-74 +lat_1=41.0333333333333 +lat_2=40.6666666666667 +x_0=300000 +y_0=0 +datum=NAD83 +units=us-ft +no_defs'

# Expansion parameters for grid generation
# These ensure grids cover the entire study area when rotated
# Adjust based on your study area size
expansion_cells <- 7000  # Number of cells to add for buffer
cell_size_ft <- 9.8425   # Cell size in feet (for US feet CRS)
# If using meters, change to appropriate value

#' ============================================================================
#' SETUP AND VALIDATION
#' ============================================================================

cat("\n=== STARTING STURLA GRID ANALYSIS ===\n")
cat("Started at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("Created output directory:", output_dir, "\n")
}

# Validate input files exist
if (!file.exists(lulc_file)) {
  stop("Land cover file not found: ", lulc_file, 
       "\nPlease check the file path and name.")
}
if (!file.exists(temp_file)) {
  stop("Temperature file not found: ", temp_file,
       "\nPlease check the file path and name.")
}

# Display configuration
cat("Configuration:\n")
cat("  Input files:\n")
cat("    Land cover/buildings:", basename(lulc_file), "\n")
cat("    Surface temperature:", basename(temp_file), "\n")
cat("  Grid parameters:\n")
cat("    Resolutions:", length(resolutions), "levels from", min(resolutions), 
    "to", max(resolutions), "meters\n")
cat("    Angles:", length(angles), "angles from", min(angles), 
    "to", max(angles), "degrees\n")
cat("    Total combinations:", length(resolutions) * length(angles), "\n")
cat("  Estimated runtime: Several hours for full analysis\n\n")

#' ============================================================================
#' DATA LOADING AND PREPROCESSING
#' ============================================================================

cat("=== LOADING RASTER DATA ===\n")

# Load land cover/building height raster
lulc_raster <- raster(lulc_file)
cat("  Land cover raster:\n")
cat("    Dimensions:", nrow(lulc_raster), "rows ×", ncol(lulc_raster), "cols\n")
cat("    Resolution:", res(lulc_raster)[1], "×", res(lulc_raster)[2], "\n")
cat("    Extent:", paste(round(extent(lulc_raster)@xmin), 
                         round(extent(lulc_raster)@xmax), 
                         round(extent(lulc_raster)@ymin), 
                         round(extent(lulc_raster)@ymax), sep=", "), "\n")

# Load surface temperature raster
temp_raster_original <- raster(temp_file)
cat("  Temperature raster:\n")
cat("    Dimensions:", nrow(temp_raster_original), "rows ×", 
    ncol(temp_raster_original), "cols\n")
cat("    Resolution:", res(temp_raster_original)[1], "×", 
    res(temp_raster_original)[2], "\n")

# Reproject temperature raster to match land cover CRS if needed
cat("\n  Checking coordinate systems...\n")
if (!compareCRS(lulc_raster, temp_raster_original)) {
  cat("  Reprojecting temperature raster to match land cover CRS...\n")
  temp_raster <- projectRaster(temp_raster_original, crs = study_area_crs)
  cat("  Reprojection complete\n")
} else {
  cat("  Coordinate systems match - no reprojection needed\n")
  temp_raster <- temp_raster_original
}

# Create expanded extent for grid generation
# This ensures rotated grids still cover the entire study area
expansion_ft <- (expansion_cells / 2) * cell_size_ft

extent_expanded <- extent(
  xmin = lulc_raster@extent@xmin - expansion_ft,
  xmax = lulc_raster@extent@xmax + expansion_ft,
  ymin = lulc_raster@extent@ymin - expansion_ft,
  ymax = lulc_raster@extent@ymax + expansion_ft
)

# Create reference raster with expanded extent
reference_raster <- raster(
  nrows = 17256, 
  ncols = 17256,
  ext = extent_expanded,
  crs = study_area_crs
)

# Calculate and display study area dimensions
area_width_m <- (extent_expanded@xmax - extent_expanded@xmin) * 0.3048
area_height_m <- (extent_expanded@ymax - extent_expanded@ymin) * 0.3048
cat("\n  Study area dimensions:", round(area_width_m), "m ×", 
    round(area_height_m), "m\n")
cat("  (Expanded extent to accommodate rotated grids)\n\n")

#' ============================================================================
#' MAIN ANALYSIS LOOP
#' ============================================================================

cat("=== STARTING GRID ANALYSIS ===\n")
cat("Processing", length(resolutions) * length(angles), 
    "grid configurations...\n")
cat("This may take several hours depending on study area size.\n\n")

# Initialize storage
all_results <- list()
result_counter <- 1

# Progress tracking
total_iterations <- length(resolutions) * length(angles)
current_iteration <- 0
start_time <- Sys.time()

# Main nested loop: iterate through all combinations
for (resolution in resolutions) {
  
  # Calculate grid dimensions for this resolution
  n_cells_x <- round(area_width_m / resolution)
  n_cells_y <- round(area_height_m / resolution)
  
  for (angle_idx in 1:length(angles)) {
    
    angle <- angles[angle_idx]
    current_iteration <- current_iteration + 1
    
    # Progress update
    elapsed <- difftime(Sys.time(), start_time, units = "mins")
    cat(sprintf("[%d/%d] Resolution: %dm, Angle: %d° (%.1f%% complete, %.1f min elapsed)\n", 
                current_iteration, total_iterations, resolution, angle,
                100 * current_iteration / total_iterations, elapsed))
    
    # Create grid for this resolution
    grid <- st_make_grid(reference_raster, n = c(n_cells_x, n_cells_y))
    
    # Apply rotation if angle > 0
    if (angle > 0) {
      # Rotation matrix function
      rotation_matrix <- function(a) matrix(c(cos(a), sin(a), -sin(a), cos(a)), 2, 2)
      
      # Get grid centroid as rotation center
      grid_center <- st_centroid(st_union(grid))
      
      # Apply rotation
      grid_rotated <- (grid - grid_center) * rotation_matrix(angle * pi / 180) + grid_center
    } else {
      grid_rotated <- grid
    }
    
    # Extract land cover counts for each class (2-10, skip background class 1)
    # Using exactextractr for accurate zonal statistics
    lulc_counts <- lapply(2:10, function(class_id) {
      exact_extract(lulc_raster, grid_rotated, function(values, coverage) {
        sum(values == class_id, na.rm = TRUE)
      })
    })
    
    # Extract mean temperature for each grid cell
    temp_values <- exact_extract(temp_raster, grid_rotated, fun = 'mean')
    
    # Process results for each grid cell
    n_grid_cells <- length(grid_rotated)
    
    for (cell_id in 1:n_grid_cells) {
      
      # Get land cover counts for this cell
      counts <- sapply(lulc_counts, function(x) x[cell_id])
      total_pixels <- sum(counts)
      
      # Skip cells with no data
      if (total_pixels == 0) next
      
      # Calculate proportions for each land cover type
      proportions <- counts / total_pixels
      
      # Build STURLA classification string
      # Only include components present in this cell
      sturla_components <- character(0)
      if (proportions[1] > 0) sturla_components <- c(sturla_components, "T")  # Trees
      if (proportions[2] > 0) sturla_components <- c(sturla_components, "G")  # Grass
      if (proportions[3] > 0) sturla_components <- c(sturla_components, "B")  # Bare soil
      if (proportions[4] > 0) sturla_components <- c(sturla_components, "W")  # Water
      if (proportions[5] > 0 || proportions[6] > 0) {
        sturla_components <- c(sturla_components, "P")  # Paved (roads + other)
      }
      if (proportions[7] > 0) sturla_components <- c(sturla_components, "L")  # Low-rise
      if (proportions[8] > 0) sturla_components <- c(sturla_components, "M")  # Mid-rise
      if (proportions[9] > 0) sturla_components <- c(sturla_components, "H")  # High-rise
      
      # Combine components into class name
      sturla_class <- paste(sturla_components, collapse = "")
      
      # Convert temperature: Kelvin * 10 to Celsius
      # NOTE: If your temperature data uses different units, modify this formula
      temp_celsius <- temp_values[cell_id] * 0.1 - 273.15
      
      # Store results for this cell
      all_results[[result_counter]] <- data.frame(
        resolution = resolution,
        angle = angle,
        cell_id = cell_id,
        temperature = temp_celsius,
        pixel_count = total_pixels,
        prop_trees = proportions[1],
        prop_grass = proportions[2],
        prop_soil = proportions[3],
        prop_water = proportions[4],
        prop_paved = proportions[5] + proportions[6],  # Combine roads and other paved
        prop_lowrise = proportions[7],
        prop_midrise = proportions[8],
        prop_highrise = proportions[9],
        sturla_class = sturla_class,
        stringsAsFactors = FALSE
      )
      
      result_counter <- result_counter + 1
    }
  }
}

#' ============================================================================
#' COMPILE AND SAVE RESULTS
#' ============================================================================

cat("\n=== COMPILING RESULTS ===\n")

# Combine all results into single data frame
final_results <- do.call(rbind, all_results)

# Data quality control: remove invalid rows
cat("  Removing invalid data...\n")
n_before <- nrow(final_results)
final_results <- final_results[!is.na(final_results$temperature), ]
final_results <- final_results[final_results$pixel_count > 0, ]
n_after <- nrow(final_results)
cat("  Removed", n_before - n_after, "invalid rows\n")

# Sort by resolution, angle, and cell_id for organization
final_results <- final_results[order(final_results$resolution, 
                                     final_results$angle, 
                                     final_results$cell_id), ]

# Generate output filename with timestamp
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output_file <- file.path(output_dir, sprintf("STURLA_analysis_%s.csv", timestamp))

# Save results
write.csv(final_results, output_file, row.names = FALSE)
cat("\n  Results saved to:", output_file, "\n")

#' ============================================================================
#' SUMMARY STATISTICS
#' ============================================================================

cat("\n=== ANALYSIS SUMMARY ===\n")

# Processing time
total_time <- difftime(Sys.time(), start_time, units = "hours")
cat(sprintf("Total processing time: %.2f hours (%.1f minutes)\n", 
            total_time, total_time * 60))

# Data summary
cat(sprintf("Total grid cells processed: %d\n", nrow(final_results)))
cat(sprintf("Valid cells with data: %d\n", sum(final_results$pixel_count > 0)))
cat(sprintf("Unique STURLA classes found: %d\n", 
            length(unique(final_results$sturla_class))))

# Temperature statistics
cat(sprintf("Temperature range: %.2f to %.2f°C\n", 
            min(final_results$temperature, na.rm = TRUE), 
            max(final_results$temperature, na.rm = TRUE)))
cat(sprintf("Mean temperature: %.2f°C\n", 
            mean(final_results$temperature, na.rm = TRUE)))

# Class distribution by resolution
cat("\nSTURLA classes by resolution:\n")
resolution_summary <- aggregate(sturla_class ~ resolution, 
                                data = final_results, 
                                FUN = function(x) length(unique(x)))
names(resolution_summary)[2] <- "n_unique_classes"
print(resolution_summary)

# Most common classes overall
cat("\nTop 10 most frequent STURLA classes:\n")
class_freq <- sort(table(final_results$sturla_class), decreasing = TRUE)
print(head(class_freq, 10))

# Save session info for reproducibility
session_file <- file.path(output_dir, sprintf("session_info_%s.txt", timestamp))
sink(session_file)
cat("STURLA Grid Analysis - Session Information\n")
cat("==========================================\n\n")
cat("Analysis completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Processing time:", sprintf("%.2f hours\n\n", total_time))
cat("Input files:\n")
cat("  Land cover:", lulc_file, "\n")
cat("  Temperature:", temp_file, "\n\n")
cat("Grid parameters:\n")
cat("  Resolutions:", paste(resolutions, collapse=", "), "\n")
cat("  Angles:", paste(angles, collapse=", "), "\n\n")
sessionInfo()
sink()

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Session info saved to:", session_file, "\n")
cat("\nOutput file is ready for use with analysis script:\n")
cat("  ", output_file, "\n")
cat("\nNext step: Use STURLA_Sensitivity_Analysis_for_publication.R\n")
cat("to generate figures and statistical analyses.\n\n")