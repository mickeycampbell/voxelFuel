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
#' @import data.table
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
#'   strength of the occlusion correction when `occlusion == "beer"`.
#' @param max_adj Numeric. Maximum allowed adjustment factor applied to any voxel
#'   during occlusion correction.
#' @param export Boolean. Defines whether or not you would like to save the
#'   resulting voxel array to a file.
#' @param file_format Character. If `export == TRUE`, there are three different
#'   file formats you can choose from (`laz`, `csv`, or `parquet`).
#' @param out_dir Character. If `export == TRUE`, the directory where you would
#'   like to save the exported file.
#' @param base_name Character. If `export == TRUE`, the base name (file name
#'   without extension) of the output file you intend to save. The file
#'   extension will be automatically assigned based on `file_format`.
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
    max_adj = 5,
    export = F,
    file_format = c("laz", "csv", "parquet"),
    out_dir,
    base_name
) {

  # check inputs
  occlusion <- match.arg(occlusion)
  file_format <- match.arg(file_format)
  stopifnot(inherits(las, "LAS"))
  stopifnot(inherits(fuel_rast, "SpatRaster"))

  # voxelize lidar
  vm <- lidR::voxel_metrics(las, ~length(Z), res = voxel_res)
  data.table::setnames(vm, "V1", "n")

  # extract biomass values
  names(fuel_rast) <- "bio"
  extracted <- terra::extract(
    fuel_rast,
    as.matrix(vm[, .(X, Y)]),
    cells = TRUE
  ) |> data.table::as.data.table()
  vm <- cbind(vm, extracted)

  # drop NA biomass cells
  vm <- vm[!is.na(bio)]

  # proportion of points per voxel (per raster cell)
  vm[, prop_all := n / sum(n), by = .(cell)]

  # filter canopy height
  vm <- vm[Z >= can_min_z & Z <= can_max_z]

  # renormalize within canopy
  vm[, prop_can := prop_all / sum(prop_all), by = .(cell)]
  vm[is.na(prop_can), prop_can := 0]

  # raw voxel biomass
  vm[, bio_vox_raw := bio * prop_can]

  # sort for cumulative calculations
  data.table::setorder(vm, cell, Z)

  # optional occlusion correction
  if (occlusion == "beer") {

    vm[, prop_above := rev(cumsum(rev(prop_all))) - prop_all, by = .(cell)]
    vm[, adj_factor := pmin(exp(beer_k * prop_above), max_adj)]
    vm[, bio_vox_temp := bio_vox_raw * adj_factor]

  } else {

    vm[, bio_vox_temp := bio_vox_raw]

  }

  # reproportion to match original raster cell value
  vm[, bio_vox := bio_vox_temp * (bio / sum(bio_vox_temp)), by = .(cell)]

  # sanity check
  chk <- vm[, .(err = abs(sum(bio_vox) - bio[1])), by = .(cell)]
  if (any(chk$err > 1e-6, na.rm = TRUE)) {
    warning("One or more cells does not aggregate back up to the original value")
  }

  # retain only desired columns
  vm <- vm[, list(X, Y, Z, n, cell, bio, bio_vox)]

  # optionally, export to file
  if (export){

    # as a laz file
    if (file_format == "laz"){

      # coerce to LAS object
      vm_las <- lidR::LAS(
        data = vm,
        crs = lidR::st_crs(las)
      )

      # add custom attributes
      vm_las <- lidR::add_lasattribute(vm_las, vm$n, "n", "Point count")
      vm_las <- lidR::add_lasattribute(vm_las, vm$cell, "cell", "Raster cell ID")
      vm_las <- lidR::add_lasattribute(vm_las, vm$bio, "bio", "Original biomass")
      vm_las <- lidR::add_lasattribute(vm_las, vm$bio_vox, "bio_vox", "Voxel biomass")

      # define output name and save
      out_file <- file.path(out_dir, paste0(base_name, ".laz"))
      lidR::writeLAS(vm_las, out_file)

    }

    # as a csv file
    else if (file_format == "csv"){

      # define output name and save
      out_file <- file.path(out_dir, paste0(base_name, ".csv"))
      data.table::fwrite(vm, out_file)

    }

    # as a parquet file
    else if (file_format == "parquet"){

      # define output name and save
      out_file <- file.path(out_dir, paste0(base_name, ".parquet"))
      nanoparquet::write_parquet(vm, file = out_file)

    }

    # not one of the three options
    else {
      stop("file_format should be one of 'csv', 'laz', or 'parquet'")
    }

  }

  # return the data.table
  return(vm[])

}
