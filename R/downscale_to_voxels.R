#' Downscale a 2D fuel biomass raster into 3D lidar-derived voxels
#'
#' This function redistributes a 2D fuel biomass raster into a vertical stack of
#' 3D voxels using lidar point density as a proxy for vertical fuel distribution.
#' Biomass within each raster cell is apportioned among vertical voxels in
#' proportion to the number of lidar returns in each voxel.
#'
#' Optionally, a Beer–Lambert-style occlusion correction can be applied to increase
#' biomass allocation to lower canopy voxels based on the fraction of points above
#' each voxel.
#'
#' @param las A `LAS` object containing normalized lidar point cloud data.
#' @param fuel_rast A single-layer `SpatRaster` containing 2D fuel biomass values.
#' @param voxel_res Numeric. Vertical and horizontal resolution of voxels
#'   (in map units).
#' @param can_min_z Numeric. Minimum canopy height (inclusive) to consider when
#'   allocating biomass.
#' @param can_max_z Numeric. Maximum canopy height (inclusive) to consider when
#'   allocating biomass.
#' @param occlusion Character. Occlusion correction method to apply. One of
#'   `"none"` (default) or `"beer"`.
#' @param beer_k Numeric. Beer–Lambert attenuation coefficient controlling the
#'   strength of the occlusion correction when `occlusion = "beer"`.
#' @param max_adj Numeric. Maximum allowed adjustment factor applied to any voxel
#'   during occlusion correction.
#'
#' @return
#' A `data.table` where each row represents a voxel and includes:
#' \itemize{
#'   \item `X`, `Y`, `Z`: voxel center coordinates
#'   \item `n`: number of lidar points in the voxel
#'   \item `cell`: index of the source biomass raster cell
#'   \item `bio`: original biomass value of the raster cell
#'   \item `bio_vox`: final biomass allocated to the voxel
#' }
#'
#' Biomass values are guaranteed (within numerical tolerance) to sum back to the
#' original raster cell values.
#'
#' @details
#' Biomass is first allocated to voxels in proportion to lidar point density.
#' If occlusion correction is enabled, voxel biomass is adjusted using an
#' exponential Beer–Lambert-style function based on the fraction of points above
#' each voxel. A final reproportioning step ensures mass conservation.
#'
#' @examples
#' \dontrun{
#' library(lidR)
#' library(terra)
#'
#' # get example data file paths
#' las_file <- system.file("extdata", "lidar_sample.laz", package = "voxelFuel")
#' fuel_file <- system.file("extdata", "fb_sample.tif", package = "voxelFuel")
#'
#' # read in the example data
#' las <- readLAS(las_file)
#' fuel_rast <- rast(fuel_file)
#'
#' # build voxel mask and align raster
#' voxel_mask <- build_voxel_mask(las)
#' fuel_rast <- align_and_mask_fuel_raster(fuel_rast, voxel_mask)
#'
#' # downscale to voxels
#' fuel_vox <- downscale_to_voxels(las, fuel_rast, voxel_res = 1, can_min_z = 1, can_max_z = 50)
#'
#' # rasterize voxel stack (optional)
#' vs <- rasterize_voxel_stack(fuel_vox, fuel_rast, out_dir = tempdir(), base_name = "example_voxel")
#' }
#'
#' @export
downscale_to_voxels <- function(
    las,
    fuel_rast,
    voxel_res = 1,
    can_min_z = 1,
    can_max_z = 50,
    occlusion = c("none", "beer"),
    beer_k = 0.5,
    max_adj = 5
) {

  # check to make sure occlusion is one of the two options
  occlusion <- base::match.arg(occlusion)

  # check input object types
  base::stopifnot(inherits(las, "LAS"))
  base::stopifnot(inherits(fuel_rast, "SpatRaster"))

  # voxelize lidar data with number of points
  vm <- lidR::voxel_metrics(las, ~length(Z), res = voxel_res)
  data.table::setnames(vm, "V1", "n")

  # get biomass values and cell ids for each voxel
  names(fuel_rast) <- "bio"
  vm <- cbind(
    vm,
    terra::extract(
      fuel_rast,
      base::as.matrix(vm[, list(X, Y)]),
      cells = TRUE
    )
  )

  # drop NA biomass cells
  vm <- vm[!base::is.na(bio)]

  # get proportion of points per voxel (per biomass pixel)
  vm[, prop_all := n / sum(n), by = cell]

  # filter voxels to canopy height range
  vm <- vm[Z >= can_min_z & Z <= can_max_z]

  # renormalize within canopy
  vm[, prop_can := prop_all / sum(prop_all), by = cell]
  vm[base::is.na(prop_can), prop_can := 0]

  # get raw voxelized biomass
  vm[, bio_vox_raw := bio * prop_can]

  # sort for cumulative calculations
  data.table::setorder(vm, cell, Z)

  # optional occlusion correction
  if (occlusion == "beer") {

    # get fraction of points above each voxel
    vm[, prop_above := rev(cumsum(rev(prop_all))) - prop_all, by = cell]

    # calculate Beer-Lambert adjustment factor
    vm[, adj_factor := base::pmin(base::exp(beer_k * prop_above), max_adj)]

    # apply to raw biomass (will no longer sum to original pixel value)
    vm[, bio_vox_temp := bio_vox_raw * adj_factor]

  } else {

    # if no correction, just copy over raw biomass
    vm[, bio_vox_temp := bio_vox_raw]

  }

  # reproportion to make sure it sums up to original pixel value
  vm[, bio_vox := bio_vox_temp * (bio / sum(bio_vox_temp)), by = cell]

  # check to see if it does, indeed, sum back up correctly
  chk <- vm[, list(err = base::abs(sum(bio_vox) - bio[1])), by = cell]
  chk <- chk[!base::is.na(err)]  # drop any cells where bio[1] is NA
  if (base::any(chk$err > 1e-6)) {
    warning("One or more cells does not aggregate back up to the original value")
  }

  # clean up columns
  vm <- vm[, list(X, Y, Z, n, cell, bio, bio_vox)]

  # return final voxel data.table
  return(vm[])
}
