This repository currently contains two utility programs for working with very large, global DEMs. Originally developed to reprocess the EarthEnv-90 DEM, the programs are the following:
  * **calcslope**: takes a global 3 arc-second DEM that has been partitioned into tiles and calculates the slope at the original resolution
  * **topostats**: decimates the input elevation and slope files to an output resolution specified at runtime, providing median elevation and slope, and fractional coverage of slopes in bins
