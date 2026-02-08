voxelFuel
================

# voxelFuel

`voxelFuel` provides tools for downscaling 2D fuel biomass maps into 3D
voxel representations using high-resolution LiDAR data. This package is
intended for researchers and practitioners working with wildfire fuel
data and LiDAR-derived canopy structure.

## Installation

You can install the development version from GitHub using `remotes`:

``` r
# install.packages("remotes") # if not already installed
remotes::install_github("mickeycampbell/voxelFuel")
```

------------------------------------------------------------------------

### Example Data

Example LiDAR and fuel raster data are included in the package:

``` r
library(voxelFuel)

las_file <- system.file("extdata/lidar.laz", package = "voxelFuel")
fuel_file <- system.file("extdata/fb.tif", package = "voxelFuel")
```

------------------------------------------------------------------------

### Quick Start / Minimal Example

This example demonstrates the core workflow:

``` r
library(lidR)
library(terra)
library(data.table)
library(voxelFuel)

# read example data
las <- readLAS(las_file)
fuel_rast <- rast(fuel_file)

# --- Build voxel eligibility mask ---
voxel_mask <- build_voxel_mask(las)

# --- Align fuel raster to voxel mask ---
fuel_rast <- align_and_mask_fuel_raster(fuel_rast, voxel_mask)

# --- Downscale 2D fuel to 3D voxels ---
fuel_vox <- downscale_to_voxels(
  las = las,
  fuel_rast = fuel_rast,
  voxel_res = 1,
  can_min_z = 1,
  can_max_z = 50,
  occlusion = "none"
)

# --- Rasterize voxel stack ---
vs <- rasterize_voxel_stack(
  vm = fuel_vox,
  fuel_rast = fuel_rast,
  out_dir = tempdir(),
  base_name = "example"
)

# Inspect outputs
vs$stack       # 3D voxel stack (all layers)
vs$stack_sum   # sum across vertical layers
vs$stack_agg   # aggregated back to original raster resolution
```

### Notes

- The package currently only supports *single LAS files*. For multiple
  tiles, process each file separately and mosaic the outputs.
  Alternatively, it could be easily parallelized with `foreach` and
  `dopar` on several lidar tiles at once.

- The biomass raster and LiDAR should be roughly aligned in space and
  projected to the same CRS.

- Voxel sizes, canopy height ranges, and occlusion parameters can be
  adjusted via function arguments.

- Functions return either a SpatRaster (for masks and stacks) or
  data.table (for voxel-level data).
