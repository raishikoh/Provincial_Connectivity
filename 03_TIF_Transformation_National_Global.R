# ==============================================================================
# 03 TIF_Transformation_National_Global.R
# 
# This script transforms national/global TIFFs to BC Albers 100m resolution 
# with standard grid alignment. Extends BC grid system up to 500km beyond 
# provincial boundary based on Common Raster Cell Origin, Shape and Sizes 
# Standard for BC.
#  
# Author: Eclipse Geomatics
# Date: 2025-07-27
# ==============================================================================

# Loading required libraries
library(terra)
library(sf)
library(fs)
library(dplyr)

# Defining file paths
input_path <- "G:/My Drive/01 Provincial Connectivity/National_Global_TIFs/Inputs"
output_folder <- "G:/My Drive/01 Provincial Connectivity/National_Global_TIFs/Outputs"
boundary_500km_path <- "G:/My Drive/01 Provincial Connectivity/AIO/BC_Outline_Buffered_Dissolved.gpkg"

# Creating output directory if it doesn't exist
dir_create(output_folder)

# Setting BC Albers parameters (EPSG:3005)
bc_albers_crs <- "EPSG:3005"
target_resolution <- 100
bc_origin_x <- 159587.5
bc_origin_y <- 173787.5

# Loading BC boundary with 500km buffer
boundary_500km <- st_read(boundary_500km_path)

# Transforming boundary to BC Albers if needed
if (st_crs(boundary_500km)$input != bc_albers_crs) {
  boundary_500km <- st_transform(boundary_500km, bc_albers_crs)
}

# Getting extended boundary extent
extended_extent <- ext(vect(boundary_500km))

# ==============================================================================
# Function to align extent to BC standard grid
# 
# @param extent_obj Spatial extent object to be aligned
# @param resolution Target resolution in meters
# @return Aligned extent object
# @examples
# aligned_ext <- align_to_bc_grid_extended(raster_extent, 100)
# ==============================================================================
align_to_bc_grid_extended <- function(extent_obj, resolution) {
  
  # Extracting extent values [xmin, xmax, ymin, ymax]
  xmin <- extent_obj[1]
  xmax <- extent_obj[2] 
  ymin <- extent_obj[3]
  ymax <- extent_obj[4]
  
  # Calculating grid alignment
  # BC grid lines are at: bc_origin + (n * resolution) where n is any integer
  
  # Finding grid line at or below xmin
  n_x_min <- floor((xmin - bc_origin_x) / resolution)
  aligned_xmin <- bc_origin_x + (n_x_min * resolution)
  
  # Finding grid line at or above xmax  
  n_x_max <- ceiling((xmax - bc_origin_x) / resolution)
  aligned_xmax <- bc_origin_x + (n_x_max * resolution)
  
  # Finding grid line at or below ymin
  n_y_min <- floor((ymin - bc_origin_y) / resolution)
  aligned_ymin <- bc_origin_y + (n_y_min * resolution)
  
  # Finding grid line at or above ymax
  n_y_max <- ceiling((ymax - bc_origin_y) / resolution)
  aligned_ymax <- bc_origin_y + (n_y_max * resolution)
  
  # Verifying perfect alignment
  x_check <- (aligned_xmin - bc_origin_x) %% resolution
  y_check <- (aligned_ymin - bc_origin_y) %% resolution
  
  # Checking for valid extent
  if (aligned_xmin >= aligned_xmax || aligned_ymin >= aligned_ymax) {
    return(extent_obj)
  }
  
  # Checking alignment accuracy
  if (abs(x_check) > 0.001 || abs(y_check) > 0.001) {
    return(extent_obj)
  }
  
  # Returning perfectly aligned extent
  return(ext(aligned_xmin, aligned_xmax, aligned_ymin, aligned_ymax))
}

# ==============================================================================
# Function to transform single raster to BC standard with extended coverage
# 
# @param input_path Path to input raster file
# @param output_path Path for output raster file
# @param clip_boundary Boundary for clipping (optional)
# @return Transformed raster object
# @examples
# result <- transform_to_bc_extended("input.tif", "output.tif", boundary)
# ==============================================================================
transform_to_bc_extended <- function(input_path, output_path, clip_boundary = NULL) {
  
  cat("Processing:", basename(input_path), "\n")
  
  # Loading raster
  input_raster <- rast(input_path)
  
  # Clipping to 500km buffer boundary first (in original projection)
  if (!is.null(clip_boundary)) {
    
    # Transforming boundary to match raster CRS for clipping
    boundary_raster_crs <- project(vect(clip_boundary), crs(input_raster))
    
    # Getting boundary extent in raster CRS
    boundary_extent <- ext(boundary_raster_crs)
    
    # Checking for overlap
    raster_extent <- ext(input_raster)
    has_overlap <- !(boundary_extent[2] < raster_extent[1] || 
                       boundary_extent[1] > raster_extent[2] ||
                       boundary_extent[4] < raster_extent[3] || 
                       boundary_extent[3] > raster_extent[4])
    
    if (!has_overlap) {
      cat("  WARNING: No overlap between raster and boundary. Skipping.\n")
      return(NULL)
    }
    
    # Cropping raster to boundary extent
    input_raster <- crop(input_raster, boundary_raster_crs)
  }
  
  # Checking if already in BC Albers, reprojecting if needed
  current_crs_info <- crs(input_raster, describe = TRUE)
  current_epsg <- ifelse(is.null(current_crs_info$code), NA, current_crs_info$code)
  
  if (is.na(current_epsg) || current_epsg != "3005") {
    input_raster <- project(input_raster, bc_albers_crs, method = "bilinear")
  }
  
  # Getting current extent after reprojection
  current_extent <- ext(input_raster)
  
  # Further refining with boundary if needed (now in BC Albers)
  if (!is.null(clip_boundary)) {
    clip_extent <- ext(vect(clip_boundary))
    
    # Final intersection in BC Albers
    overlap_result <- try({
      intersected_extent <- intersect(current_extent, clip_extent)
      if (is.null(intersected_extent) || any(is.na(as.vector(intersected_extent)))) {
        stop("No valid overlap")
      }
      extent_values <- as.vector(intersected_extent)
      if (length(extent_values) != 4 || extent_values[1] >= extent_values[2] || 
          extent_values[3] >= extent_values[4]) {
        stop("Invalid extent dimensions")
      }
      intersected_extent
    }, silent = TRUE)
    
    if (inherits(overlap_result, "try-error")) {
      final_extent <- current_extent
    } else {
      final_extent <- overlap_result
    }
  } else {
    final_extent <- current_extent
  }
  
  # Aligning extent to BC standard grid
  aligned_extent <- align_to_bc_grid_extended(final_extent, target_resolution)
  
  # Creating template and resampling
  grid_template <- rast(extent = aligned_extent, 
                        resolution = target_resolution,
                        crs = bc_albers_crs)
  
  resampled_raster <- resample(input_raster, grid_template, method = "bilinear")
  
  # Final precise clipping to boundary
  if (!is.null(clip_boundary)) {
    boundary_vector <- vect(clip_boundary)
    resampled_raster <- mask(resampled_raster, boundary_vector)
    resampled_raster <- crop(resampled_raster, boundary_vector)
  }
  
  # Writing output with metadata
  output_description <- paste0("BC Extended Grid Standard (500km buffer) - ",
                               "Origin: 159587.5,173787.5 - ",
                               "Resolution: 100m - ",
                               "Method: Clip-first then bilinear resampling to extended BC grid")
  
  writeRaster(resampled_raster, output_path, 
              overwrite = TRUE,
              gdal = c("COMPRESS=LZW", "TILED=YES",
                       "AREA_OR_POINT=AREA",
                       paste0("TIFFTAG_IMAGEDESCRIPTION=", output_description)))
  
  cat("  Completed:", basename(output_path), "\n")
  
  return(resampled_raster)
}

# Determining if input is single file or directory
is_single_file <- !dir.exists(input_path) && file.exists(input_path)

# Getting list of TIFF files to process
if (is_single_file) {
  # Processing single file
  if (!grepl("\\.(tif|tiff)$", input_path, ignore.case = TRUE)) {
    stop("Input file is not a TIFF: ", input_path)
  }
  tiff_files <- input_path
} else {
  # Processing directory
  # Using robust file discovery
  tiff_files <- tryCatch({
    dir_ls(input_path, regexp = "\\.(tif|tiff)$", ignore.case = TRUE)
  }, error = function(e) {
    list.files(input_path, pattern = "\\.(tif|tiff)$", 
               ignore.case = TRUE, full.names = TRUE)
  })
}

# Checking if files were found
if (length(tiff_files) == 0) {
  stop("No TIFF files found in: ", input_path)
}

cat("Found", length(tiff_files), "TIFF file(s) to process\n\n")

# Processing each file
processing_results <- list()

for (i in seq_along(tiff_files)) {
  
  current_input_file <- tiff_files[i]
  output_filename <- paste0(tools::file_path_sans_ext(basename(current_input_file)), 
                            "_bc_extended_100m.tif")
  current_output_path <- file.path(output_folder, output_filename)
  
  # Transforming raster
  processing_result <- tryCatch({
    transform_to_bc_extended(current_input_file, current_output_path, boundary_500km)
    "SUCCESS"
  }, error = function(e) {
    error_message <- paste("ERROR:", e$message)
    cat("ERROR processing", basename(current_input_file), ":", e$message, "\n")
    error_message
  })
  
  processing_results[[basename(current_input_file)]] <- processing_result
}

# Printing processing summary
cat("\n=== PROCESSING COMPLETE ===\n")
for (filename in names(processing_results)) {
  cat("", filename, ":", processing_results[[filename]], "\n")
}

# Verifying BC standard compliance for sample output file
output_files <- list.files(output_folder, pattern = "\\.tif$", full.names = TRUE)

if (length(output_files) > 0) {
  
  cat("\n=== BC GRID COMPLIANCE CHECK ===\n")
  sample_output_file <- output_files[1]
  sample_output_raster <- rast(sample_output_file)
  
  # Basic information
  sample_crs_info <- crs(sample_output_raster, describe = TRUE)
  sample_epsg <- ifelse(is.null(sample_crs_info$code), NA, sample_crs_info$code)
  
  cat("Sample file:", basename(sample_output_file), "\n")
  cat("Resolution:", res(sample_output_raster), "meters\n")
  
  # Checking grid alignment
  sample_extent <- ext(sample_output_raster)
  x_alignment_offset <- (sample_extent[1] - bc_origin_x) %% target_resolution
  y_alignment_offset <- (sample_extent[3] - bc_origin_y) %% target_resolution
  
  # Compliance checks
  grid_compliant <- abs(x_alignment_offset) < 0.01 && abs(y_alignment_offset) < 0.01
  projection_compliant <- !is.na(sample_epsg) && sample_epsg == "3005"
  resolution_compliant <- all(abs(res(sample_output_raster) - target_resolution) < 0.01)
  
  cat("Grid alignment:", ifelse(grid_compliant, "COMPLIANT", "NOT COMPLIANT"), "\n")
  cat("Projection (BC Albers):", ifelse(projection_compliant, "COMPLIANT", "NOT COMPLIANT"), "\n")
  cat("Resolution (100m):", ifelse(resolution_compliant, "COMPLIANT", "NOT COMPLIANT"), "\n")
  
  # Coverage information
  extent_values <- as.vector(sample_extent)
  width_kilometers <- round((extent_values[2] - extent_values[1]) / 1000, 1)
  height_kilometers <- round((extent_values[4] - extent_values[3]) / 1000, 1)
  cat("Coverage area:", width_kilometers, "km x", height_kilometers, "km\n")
}

cat("\nAll files processed and saved to:", output_folder, "\n")