#!/bin/bash
#
#---------------------------------------
# bash script in cyclone
#---------------------------------------
#set -xe
#
# loading modules in cyclone for bash
#. /shared/apps/easybuild/modules/all
module load CDO/1.9.5-foss-2018b

#expname=historical
#expname=rcp85
expname=rcp60
modelname=NorESM1-M

dirin1=/shared/projects/gfi/share/ModData2/CMIP5/atmos/${expname}/va/6hrPlev/${modelname}/r1i1p1
#infile_va1=va_6hrPlev_${modelname}_${expname}_r1i1p1_197*
#infile_va2=va_6hrPlev_${modelname}_${expname}_r1i1p1_198*
#infile_va3=va_6hrPlev_${modelname}_${expname}_r1i1p1_199*
infile_va1=va_6hrPlev_${modelname}_${expname}_r1i1p1_207*
infile_va2=va_6hrPlev_${modelname}_${expname}_r1i1p1_208*
infile_va3=va_6hrPlev_${modelname}_${expname}_r1i1p1_209*
dirin2=/shared/projects/gfi/share/ModData2/CMIP5/atmos/${expname}/ua/6hrPlev/${modelname}/r1i1p1
#infile_ua1=ua_6hrPlev_${modelname}_${expname}_r1i1p1_197*
#infile_ua2=ua_6hrPlev_${modelname}_${expname}_r1i1p1_198*
#infile_ua3=ua_6hrPlev_${modelname}_${expname}_r1i1p1_199*
infile_ua1=ua_6hrPlev_${modelname}_${expname}_r1i1p1_207*
infile_ua2=ua_6hrPlev_${modelname}_${expname}_r1i1p1_208*
infile_ua3=ua_6hrPlev_${modelname}_${expname}_r1i1p1_209*

#ln -sf ${dirin}/{INFILE} .
#ln -sf ${dirin}/${FILEATM} .

# Merge all the files in time from 2070 to 2100 for the RCP8.5 data
cdo mergetime ${dirin1}/${infile_va1} ${dirin1}/${infile_va2} ${dirin1}/${infile_va3} va_6hrPlev_${modelname}_${expname}_r1i1p1_2071_2100.nc 
cdo mergetime ${dirin2}/${infile_ua1} ${dirin2}/${infile_ua2} ${dirin2}/${infile_ua3} ua_6hrPlev_${modelname}_${expname}_r1i1p1_2071_2100.nc

# Select the pressure level of 85000Pa and merge all the files in a single one
cdo sellevel,85000 va_6hrPlev_${modelname}_${expname}_r1i1p1_2071_2100.nc va_6hrPlev_${modelname}_${expname}_r1i1p1_2071_2100_850hPa.nc
cdo sellevel,85000 ua_6hrPlev_${modelname}_${expname}_r1i1p1_2071_2100.nc ua_6hrPlev_${modelname}_${expname}_r1i1p1_2071_2100_850hPa.nc

# remove the files after using them to not exceed the quota
rm va_6hrPlev_${modelname}_${expname}_r1i1p1_2071_2100.nc
rm ua_6hrPlev_${modelname}_${expname}_r1i1p1_2071_2100.nc

# Interpolate to 12h data
cdo inttime,2071-01-01,00:00:00,12hour va_6hrPlev_${modelname}_${expname}_r1i1p1_2071_2100_850hPa.nc va_6hrPlev_${modelname}_${expname}_r1i1p1_2071_2100_850hPa_6h.nc
cdo inttime,2071-01-01,00:00:00,12hour ua_6hrPlev_${modelname}_${expname}_r1i1p1_2071_2100_850hPa.nc ua_6hrPlev_${modelname}_${expname}_r1i1p1_2071_2100_850hPa_6h.nc

rm *850hPa.nc

# Calculate the climatological mean 12hmean
cdo yhourmean va_6hrPlev_${modelname}_${expname}_r1i1p1_2071_2100_850hPa_6h.nc va_6hrPlev_${modelname}_${expname}_r1i1p1_2071_2100_850hPa_6hmean.nc
cdo yhourmean ua_6hrPlev_${modelname}_${expname}_r1i1p1_2071_2100_850hPa_6h.nc ua_6hrPlev_${modelname}_${expname}_r1i1p1_2071_2100_850hPa_6hmean.nc

rm *_6h.nc

# Interpolate to GREB model grid (94x48) 
cdo remapcon,r96x48 va_6hrPlev_${modelname}_${expname}_r1i1p1_2071_2100_850hPa_6hmean.nc va_6hrPlev_${modelname}_${expname}_r1i1p1_2071_2100_850hPa_6hmean_GREBgrid.nc
cdo remapcon,r96x48 ua_6hrPlev_${modelname}_${expname}_r1i1p1_2071_2100_850hPa_6hmean.nc ua_6hrPlev_${modelname}_${expname}_r1i1p1_2071_2100_850hPa_6hmean_GREBgrid.nc

mv va_6hrPlev_${modelname}_${expname}_r1i1p1_2071_2100_850hPa_6hmean_GREBgrid.nc v_${modelname}_${expname}.nc
mv ua_6hrPlev_${modelname}_${expname}_r1i1p1_2071_2100_850hPa_6hmean_GREBgrid.nc u_${modelname}_${expname}.nc

# delete unnecessary files
rm *mean.nc

#EOF


