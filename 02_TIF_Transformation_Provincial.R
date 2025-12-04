# ==============================================================================
#' 02 TIF_Transformation_Provincial.R
#' 
#' This script transforms multiple TIFF files to BC Albers projection with 
#' 100m resolution and standard grid alignment. Follows BC Government 
#' Common Raster Cell Origin, Shape and Sizes Standard.
#'  
#' Author: Eclipse Geomatics
#' Date: 2025-07-27
# ==============================================================================

# Loading required libraries --------------------------------------------------
library(terra)
library(sf)
library(fs)
library(dplyr)

# Setting up file paths -------------------------------------------------------
input_path <- "G:/My Drive/01 Provincial Connectivity/Provincial_TIFs/Inputs"
output_folder <- "G:/My Drive/01 Provincial Connectivity/Provincial_TIFs/Outputs"
boundary_path <- "G:/My Drive/01 Provincial Connectivity/AIO/BC_Outline.gpkg"

# Checking if input path exists
if (!file.exists(input_path)) {
  stop("Input path does not exist: ", input_path)
}

# Determining if input is a single file or directory
is_single_file <- file.info(input_path)$isdir == FALSE

# Creating output directory if it doesn't exist
dir_create(output_folder)

# Setting BC Albers projection parameters -------------------------------------
bc_albers_crs <- "EPSG:3005"
target_resolution <- 100  # 100 meter resolution
bc_origin_x <- 159587.5   # BC standard origin easting coordinate
bc_origin_y <- 173787.5   # BC standard origin northing coordinate

# Loading area of interest boundary -------------------------------------------
aoi <- st_read(boundary_path, quiet = TRUE)

# Transforming AOI to BC Albers projection if needed
if (st_crs(aoi)$input != bc_albers_crs) {
  aoi <- st_transform(aoi, bc_albers_crs)
}

# Getting AOI extent for clipping operations
aoi_extent <- ext(vect(aoi))

# ==============================================================================
#' Function to align raster extent to BC standard grid
#' 
#' @param extent_obj Terra extent object to be aligned
#' @param resolution Target resolution in meters (typically 100)
#' @return Aligned terra extent object
#' @examples
#' aligned_ext <- align_to_bc_grid(raster_extent, 100)
# ==============================================================================
align_to_bc_grid <- function(extent_obj, resolution) {
  
  # Extracting extent coordinates
  xmin <- extent_obj[1]
  xmax <- extent_obj[2] 
  ymin <- extent_obj[3]
  ymax <- extent_obj[4]
  
  # Calculating grid-aligned boundaries
  # BC grid lines are at: bc_origin + (n * resolution) where n is any integer
  
  # Finding the grid line at or below xmin
  n_x_min <- floor((xmin - bc_origin_x) / resolution)
  aligned_xmin <- bc_origin_x + (n_x_min * resolution)
  
  # Finding the grid line at or above xmax  
  n_x_max <- ceiling((xmax - bc_origin_x) / resolution)
  aligned_xmax <- bc_origin_x + (n_x_max * resolution)
  
  # Finding the grid line at or below ymin
  n_y_min <- floor((ymin - bc_origin_y) / resolution)
  aligned_ymin <- bc_origin_y + (n_y_min * resolution)
  
  # Finding the grid line at or above ymax
  n_y_max <- ceiling((ymax - bc_origin_y) / resolution)
  aligned_ymax <- bc_origin_y + (n_y_max * resolution)
  
  # Verifying perfect alignment
  x_check <- (aligned_xmin - bc_origin_x) %% resolution
  y_check <- (aligned_ymin - bc_origin_y) %% resolution
  
  # Ensuring having a valid extent
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
#' Function to transform single raster to BC standard
#' 
#' @param input_path Path to input TIFF file
#' @param output_path Path for output transformed file
#' @param clip_boundary Boundary object for clipping (optional)
#' @return Transformed raster object
#' @examples
#' result <- transform_to_bc_standard("input.tif", "output.tif", boundary)
# ==============================================================================
transform_to_bc_standard <- function(input_path, output_path, clip_boundary = NULL) {
  
  cat("Processing:", basename(input_path), "\n")
  
  # Loading input raster
  input_raster <- rast(input_path)
  
  # Checking current coordinate reference system
  current_crs <- crs(input_raster, describe = TRUE)
  current_epsg <- ifelse(is.null(current_crs$code), NA, current_crs$code)
  
  # Reprojecting to BC Albers if needed
  if (is.na(current_epsg) || current_epsg != "3005") {
    input_raster <- project(input_raster, bc_albers_crs, method = "bilinear")
  }
  
  # Getting current raster extent
  current_extent <- ext(input_raster)
  
  # Determining final extent for processing
  if (!is.null(clip_boundary)) {
    clip_extent <- ext(vect(clip_boundary))
    
    # Checking if raster and boundary overlap
    overlap_test <- try({
      intersected_extent <- intersect(current_extent, clip_extent)
      if (is.null(intersected_extent) || any(is.na(as.vector(intersected_extent))) || 
          intersected_extent[1] >= intersected_extent[3] || intersected_extent[2] >= intersected_extent[4]) {
        stop("No valid overlap")
      }
      intersected_extent
    }, silent = TRUE)
    
    # Using intersected extent if overlap exists, otherwise use full extent
    if (inherits(overlap_test, "try-error")) {
      final_extent <- current_extent
    } else {
      final_extent <- overlap_test
    }
  } else {
    final_extent <- current_extent
  }
  
  # Aligning extent to BC standard grid
  aligned_extent <- align_to_bc_grid(final_extent, target_resolution)
  
  # Creating template raster with BC standard specifications
  template_raster <- rast(extent = aligned_extent, 
                          resolution = target_resolution,
                          crs = bc_albers_crs)
  
  # Resampling input raster to match template
  resampled_raster <- resample(input_raster, template_raster, method = "bilinear")
  
  # Clipping to boundary if provided
  if (!is.null(clip_boundary)) {
    boundary_vector <- vect(clip_boundary)
    resampled_raster <- mask(resampled_raster, boundary_vector)
    resampled_raster <- crop(resampled_raster, boundary_vector)
  }
  
  # Writing output file with BC standard metadata
  writeRaster(resampled_raster, output_path, 
              overwrite = TRUE,
              gdal = c("COMPRESS=LZW", "TILED=YES",
                       "AREA_OR_POINT=AREA",
                       paste0("TIFFTAG_IMAGEDESCRIPTION=", 
                              "BC Government Standard Compliant Raster - ",
                              "Origin: 159587.5,173787.5 - ",
                              "Resolution: 100m - ",
                              "Method: Bilinear resampling to BC grid")))
  
  cat("Completed:", basename(output_path), "\n")
  
  return(resampled_raster)
}

# Finding TIFF files to process -----------------------------------------------
if (is_single_file) {
  # Processing single file
  if (!grepl("\\.(tif|tiff)$", input_path, ignore.case = TRUE)) {
    stop("Input file is not a TIFF: ", input_path)
  }
  tiff_files <- input_path
} else {
  # Processing directory of files
  tiff_files <- list.files(input_path, 
                           pattern = "\\.(tif|tiff)$", 
                           ignore.case = TRUE, 
                           full.names = TRUE)
}

# Checking if files were found
if (length(tiff_files) == 0) {
  stop("No TIFF files found in: ", input_path)
}

cat("Found", length(tiff_files), "TIFF file(s) to process\n\n")

# Processing each TIFF file ---------------------------------------------------
processing_results <- list()

for (i in seq_along(tiff_files)) {
  
  # Setting up file paths
  input_file <- tiff_files[i]
  output_filename <- paste0(tools::file_path_sans_ext(basename(input_file)), 
                            "_bc100m.tif")
  output_path <- file.path(output_folder, output_filename)
  
  # Processing current file
  processing_attempt <- try({
    transform_to_bc_standard(input_file, output_path, aoi)
    processing_results[[basename(input_file)]] <- "SUCCESS"
  }, silent = TRUE)
  
  # Recording any errors
  if (inherits(processing_attempt, "try-error")) {
    error_message <- as.character(processing_attempt)
    processing_results[[basename(input_file)]] <- paste("ERROR:", error_message)
    cat("Error processing", basename(input_file), "\n")
  }
}

# Displaying processing summary -----------------------------------------------
cat("\n=== PROCESSING SUMMARY ===\n")
successful_files <- 0
failed_files <- 0

for (file_name in names(processing_results)) {
  result_status <- processing_results[[file_name]]
  cat(file_name, ":", result_status, "\n")
  
  if (result_status == "SUCCESS") {
    successful_files <- successful_files + 1
  } else {
    failed_files <- failed_files + 1
  }
}

cat("\nFiles processed successfully:", successful_files, "\n")
cat("Files with errors:", failed_files, "\n")

# Verifying BC standard compliance ---------------------------------------------
output_files <- list.files(output_folder, pattern = "\\.tif$", full.names = TRUE)

if (length(output_files) > 0) {
  
  cat("\n=== BC STANDARD COMPLIANCE CHECK ===\n")
  sample_file <- output_files[1]
  sample_raster <- rast(sample_file)
  
  # Checking projection compliance
  sample_crs <- crs(sample_raster, describe = TRUE)
  sample_epsg <- ifelse(is.null(sample_crs$code), NA, sample_crs$code)
  
  cat("Sample file:", basename(sample_file), "\n")
  cat("Projection:", ifelse(!is.na(sample_epsg) && sample_epsg == "3005", 
                            "BC Albers (COMPLIANT)", "NOT BC Albers (NON-COMPLIANT)"), "\n")
  
  # Checking resolution compliance
  sample_resolution <- res(sample_raster)
  resolution_compliant <- all(abs(sample_resolution - target_resolution) < 0.01)
  cat("Resolution:", ifelse(resolution_compliant, 
                            "100m x 100m (COMPLIANT)", "NOT 100m (NON-COMPLIANT)"), "\n")
  
  # Checking grid alignment compliance
  sample_extent <- ext(sample_raster)
  x_alignment_offset <- (sample_extent[1] - bc_origin_x) %% target_resolution
  y_alignment_offset <- (sample_extent[3] - bc_origin_y) %% target_resolution
  alignment_compliant <- abs(x_alignment_offset) < 0.01 && abs(y_alignment_offset) < 0.01
  cat("Grid alignment:", ifelse(alignment_compliant, 
                                "BC Standard (COMPLIANT)", "NOT aligned (NON-COMPLIANT)"), "\n")
  
  # Overall compliance status
  overall_compliant <- !is.na(sample_epsg) && sample_epsg == "3005" && 
    resolution_compliant && alignment_compliant
  cat("Overall compliance:", ifelse(overall_compliant, "COMPLIANT", "NON-COMPLIANT"), "\n")
}

cat("\nAll processing complete. Output files saved to:", output_folder, "\n")