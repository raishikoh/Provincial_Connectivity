# Provincial Connectivity

This repository contains scripts for processing and standardizing spatial data required for the Provincial Corridor Project to BC Government standards. The scripts handle conversion of shapefiles to GeoPackage format and transformation of raster data to BC Albers projection with standardized grid alignment. Below is the sequence of scripts with detailed descriptions:

## Sequence of Scripts

1. **01_SHP_to_GPKG_Conversion.R**  
   Converts all shapefiles in a directory to GeoPackage format and transforms them to the EPSG:3005 coordinate system. The script processes multiple shapefiles with error handling to continue processing if individual files fail.

2. **02_TIF_Transformation_Provincial.R**  
   Transforms multiple TIFF files to BC Albers projection with 100m resolution and standard grid alignment. The script follows BC Government Common Raster Cell Origin, Shape, and Sizes Standard, clips data to BC boundaries, and implements resampling. It includes grid alignment functions and compliance verification to ensure outputs meet BC standards.

3. **03_TIF_Transformation_National_Global.R**  

   Transforms national and global TIFF datasets to BC Albers 100m resolution with standard grid alignment. The script extends BC grid system up to 500km beyond provincial boundary and uses a clip-first approach for processing large datasets. It maintains BC grid alignment for regional consistency while handling very large raster files efficiently.
