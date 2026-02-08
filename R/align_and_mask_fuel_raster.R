#' Align and mask a fuel biomass raster to lidar-derived processing mask
#'
#' This function reprojects, crops, and masks a 2D fuel biomass raster so that
#' it is spatially aligned with a lidar-derived voxel processing mask.
#' The output raster matches the lidar extent, resolution, and CRS, with
#' ineligible areas set to `NA`.
#'
#' @param fuel_rast A `SpatRaster` containing fuel biomass values.
#'   Must have a valid coordinate reference system.
#' @param voxel_mask A `SpatRaster` produced by
#'   `\link{build_voxel_mask}`, defining lidar-eligible areas
#'   (values of 1) and ineligible areas (`NA`).
#' @param method Character string specifying the resampling method used
#'   if reprojection is required. Passed to `terra::project()`.
#'   Defaults to `"bilinear"`.
#'
#' @return A `SpatRaster` representing the fuel biomass raster
#'   reprojected (if necessary), cropped to the lidar extent, and masked
#'   to eligible voxels.
#'
#' @details
#' If the coordinate reference system (CRS) of `fuel_rast` does not
#' match that of `voxel_mask`, the fuel raster is reprojected to
#' the voxel mask CRS prior to cropping and masking.
#'
#' @examples
#' \dontrun{
#' library(terra)
#' fuel_file <- system.file("extdata", "fb_sample.tif", package = "voxelFuel")
#' las_file <- system.file("extdata", "lidar_sample.laz", package = "voxelFuel")
#' fuel_rast <- rast(fuel_file)
#' las <- readLAS(las_file)
#' voxel_mask <- build_voxel_mask(las, res = 10, min_density = 10, min_height = 5)
#' aligned_fuel <- align_and_mask_fuel_raster(fuel_rast, voxel_mask)
#' }
#'
#' @export
align_and_mask_fuel_raster <- function(
    fuel_rast,
    voxel_mask,
    method = "bilinear"
) {

  # check input object types
  base::stopifnot(inherits(fuel_rast, "SpatRaster"))
  base::stopifnot(inherits(voxel_mask, "SpatRaster"))

  # get crs
  voxel_mask_crs <- terra::crs(voxel_mask)
  fuel_crs <- terra::crs(fuel_rast)

  # check if either crs is NA
  if (base::is.na(fuel_crs)) {
    base::stop("fuel_rast has no CRS")
  } else if (base::is.na(voxel_mask_crs)){
    base::stop("voxel_mask has no CRS")
  }

  # if fuel_crs doesn't match voxel_mask_crs, reproject it
  if (!terra::same.crs(fuel_rast, voxel_mask_crs)) {
    fuel_rast <- terra::project(fuel_rast, voxel_mask, method = method)
  }

  # crop fuel_rast to the extent of the lidar data
  voxel_mask <- terra::crop(voxel_mask, fuel_rast)
  fuel_rast <- terra::crop(fuel_rast, voxel_mask, mask = TRUE)

  # return aligned fuel_rast
  return(fuel_rast)
}
