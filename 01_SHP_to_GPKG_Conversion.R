# ==============================================================================
#' 01 SHP_to_GPKG_Conversion.R
#' 
#' This script converts all shapefiles in a directory to GeoPackage format
#' and transforms them to EPSG:3005 coordinate system
#'  
#' Author: Eclipse Geomatics
#' Date: 2025-07-23
# ==============================================================================

# Loading required libraries
library(sf)
library(dplyr)

# Defining input and output directories
input_directory <- "G:/My Drive/01 Provincial Connectivity/SHP_TO_GPKG/Inputs"
output_directory <- "G:/My Drive/01 Provincial Connectivity/SHP_TO_GPKG/Outputs"

# Creating output directory if it doesn't exist
if (!dir.exists(output_directory)) {
  dir.create(output_directory, recursive = TRUE)
}

# Finding all shapefile (.shp) files in the input directory
shapefile_list <- list.files(input_directory, pattern = "\\.shp$", full.names = TRUE, ignore.case = TRUE)

# Checking if any shapefiles were found
if (length(shapefile_list) == 0) {
  stop("No SHP files found in the input directory.")
}

# Setting target coordinate reference system to EPSG:3005
target_coordinate_system <- st_crs(3005)

# Processing each shapefile
for (current_shapefile in shapefile_list) {
  
  # Using error handling to continue processing other files if one fails
  tryCatch({
    
    # Reading the shapefile with all attributes preserved
    spatial_data <- st_read(current_shapefile, stringsAsFactors = FALSE, quiet = TRUE)
    
    # Getting the current coordinate reference system
    current_coordinate_system <- st_crs(spatial_data)
    
    # Handling coordinate system transformation
    # If CRS is missing, assume it's in geographic coordinates (EPSG:4326)
    if (is.na(current_coordinate_system$input)) {
      st_crs(spatial_data) <- 4326
      current_coordinate_system <- st_crs(4326)
    }
    
    # Transforming to target CRS if needed
    if (current_coordinate_system != target_coordinate_system) {
      spatial_data <- st_transform(spatial_data, crs = 3005)
    }
    
    # Creating output filename by replacing .shp extension with .gpkg
    shapefile_basename <- tools::file_path_sans_ext(basename(current_shapefile))
    output_filename <- file.path(output_directory, paste0(shapefile_basename, ".gpkg"))
    
    # Writing to GeoPackage format
    # delete_dsn = TRUE will overwrite existing files
    st_write(spatial_data, output_filename, 
             delete_dsn = TRUE,
             quiet = TRUE)
    
    # Providing basic progress feedback
    cat("Converted:", basename(current_shapefile), "->", basename(output_filename), "\n")
    
  }, error = function(error_message) {
    # Printing error message and continuing with next file
    cat("ERROR processing", basename(current_shapefile), ":", error_message$message, "\n")
  })
}

cat("Conversion process completed.\n")