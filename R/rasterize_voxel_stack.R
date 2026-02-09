#' Rasterize voxel biomass into a vertical stack of 2D rasters
#'
#' Converts voxel-level biomass estimates into a multi-layer raster stack,
#' where each layer represents a horizontal slice of the canopy at a fixed
#' height interval. Additional rasters are produced by summing the stack
#' vertically and aggregating back to the resolution of the input fuel raster.
#'
#' @param vm A `data.table` of voxel-level biomass produced by
#'   \code{\link{downscale_to_voxels}}, containing at least `X`, `Y`, `Z`, and
#'   `bio_vox` columns.
#' @param fuel_rast A `SpatRaster` used as a spatial reference for extent,
#'   resolution, and coordinate reference system.
#' @param out_dir Character. Output directory for raster files.
#' @param base_name Character. Base name used for output raster files.
#' @param voxel_res Numeric. Vertical and horizontal resolution of voxels
#'   (in map units).
#' @param can_min_z Numeric. Minimum canopy height included in the stack.
#' @param can_max_z Numeric. Maximum canopy height included in the stack.
#'
#' @details
#' This function writes three raster files to disk:
#' \itemize{
#'   \item A multi-layer raster stack of voxel biomass (`*_stack.tif`)
#'   \item A vertically summed raster (`*_sum.tif`)
#'   \item An aggregated raster matching the resolution of `fuel_rast` (`*_agg.tif`)
#' }
#'
#' Empty voxels are filled with zero to ensure consistent layer structure.
#'
#' @return
#' A named list of `SpatRasters`:
#' \itemize{
#'   \item `stack`
#'   \item `stack_sum`
#'   \item `stack_agg`
#' }
#'
#' @examples
#' \dontrun{
#' vs <- rasterize_voxel_stack(
#'   vm = fuel_vox,
#'   fuel_rast = aligned_fuel,
#'   out_dir = tempdir(),
#'   base_name = "example_voxel",
#'   voxel_res = 1,
#'   can_min_z = 1,
#'   can_max_z = 50
#' )
#' }
#'
#' @export
rasterize_voxel_stack <- function(
    vm,
    fuel_rast,
    out_dir,
    base_name,
    voxel_res = 1,
    can_min_z = 1,
    can_max_z = 50
) {

  # get number and names of layers
  z_vals <- base::seq(can_min_z, can_max_z, by = voxel_res)
  n_layers <- base::length(z_vals)
  layer_names <- base::paste0("z_", z_vals)

  # create empty raster stack aligned to fuel raster
  stack <- terra::rast(
    extent = terra::ext(fuel_rast),
    resolution = voxel_res,
    crs = terra::crs(fuel_rast),
    nlyrs = n_layers
  )
  names(stack) <- layer_names

  # fill layers
  pb <- progress::progress_bar$new(
    total = n_layers,
    format = ":current/:total [:bar] ETA: :eta"
  )
  for (i in seq_along(z_vals)) {
    pb$tick()
    z_val <- z_vals[i]
    vm_z <- vm[Z == z_val]
    if (base::nrow(vm_z) == 0) next
    xy <- base::as.matrix(vm_z[, .(X,Y)])
    cell_ids <- terra::cellFromXY(stack[[i]], xy)
    stack[[i]][cell_ids] <- vm_z[["bio_vox"]]
    stack[[i]][base::is.na(stack[[i]])] <- 0
  }

  # ensure directories exists
  base::dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  # write output file
  stack_out_file <- base::file.path(out_dir, base::paste0(base_name, "_stack.tif"))
  terra::writeRaster(stack, stack_out_file, overwrite = TRUE)

  # sum it up
  sum_out_file <- base::file.path(out_dir, base::paste0(base_name, "_sum.tif"))
  stack_sum <- terra::app(stack, sum, na.rm = TRUE)
  terra::writeRaster(stack_sum, sum_out_file, overwrite = TRUE)

  # aggregate
  agg_out_file <- base::file.path(out_dir, base::paste0(base_name, "_agg.tif"))
  agg_fact <- terra::res(fuel_rast)[1] / voxel_res
  stack_agg <- terra::aggregate(stack_sum, agg_fact, sum)
  terra::writeRaster(stack_agg, agg_out_file, overwrite = TRUE)

  # returns
  return(list(
    stack = stack,
    stack_sum = stack_sum,
    stack_agg = stack_agg
  ))
}
