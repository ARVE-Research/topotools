#!/bin/bash

source ~/.bash_profile

cd calcslope

./calcslope -180/180/-90/90 ../EarthEnv-DEM90.nc globaldem_wslope.nc

cd ..

cd stats

#./topostats -180/180/-90/90 300 ../calcslope/globaldem_wslope.nc globaldem05m.nc 
./topostats -180/180/-90/90 300 infiles.namelist globaldem05m.nc
