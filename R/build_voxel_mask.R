#' Build a voxel processing mask based on lidar point density and canopy height
#'
#' This function creates a raster mask identifying areas eligible for voxel-based
#' fuel downscaling. A cell is considered eligible if it meets both a minimum lidar
#' point density threshold and a minimum canopy height threshold.
#'
#' Lidar point density is computed as points per square meter, and canopy height
#' is estimated using a high quantile (0.999) of point elevations within each cell.
#'
#' @param las A `LAS` object containing height-normalized lidar point cloud data.
#' @param res Numeric. Resolution (in map units) used to compute point density and
#'   canopy height rasters.
#' @param min_density Numeric. Minimum lidar point density (points per square meter)
#'   required for a cell to be considered eligible.
#' @param min_height Numeric. Minimum canopy height (same units as `Z`) required for
#'   a cell to be considered eligible.
#'
#' @return A `SpatRaster` where eligible cells are coded as `1` and ineligible cells
#'   are `NA`. The raster is aligned to the input `las` extent and resolution.
#'
#' @details
#' This mask is typically used to constrain subsequent voxel-based fuel allocation
#' to areas with sufficient lidar support and meaningful canopy structure.
#'
#' @examples
#' \dontrun{
#' library(lidR)
#' las_file <- system.file("extdata", "lidar_sample.laz", package = "voxelFuel")
#' las <- readLAS(las_file)
#' voxel_mask <- build_voxel_mask(las, res = 10, min_density = 10, min_height = 5)
#' }
#'
#' @export
build_voxel_mask <- function(
    las,
    res = 10,
    min_density = 10,
    min_height = 5
) {
  # check input
  base::stopifnot(inherits(las, "LAS"))

  # compute point density raster (points per square meter)
  pd <- lidR::pixel_metrics(
    las = las,
    func = ~length(Z),
    res = res
  )
  pd <- pd / (res^2)

  # compute canopy height model (CHM) raster
  chm <- lidR::pixel_metrics(
    las = las,
    func = ~stats::quantile(Z, probs = 0.999),
    res = res
  )
  names(chm) <- "height"

  # build final mask: 1 = eligible, NA = not eligible
  mask <- terra::ifel(pd >= min_density & chm >= min_height, 1, NA)

  # return the mask
  return(mask)
}
