# ------------------------------------------------------------------------------------------------------------------- #
# Define stuffz (area plotted, size of letters, etc)
gmt set FORMAT_GEO_MAP D
gmt set FORMAT_GEO_MAP D
gmt set FONT_ANNOT_PRIMARY Helvetica
gmt set FONT_ANNOT_PRIMARY 8
gmt set FONT_LABEL Helvetica
gmt set LABEL_FONT_SIZE 7
gmt set MAP_FRAME_TYPE plain

#gmt basemap -JP?+z -R0/100/30/80 -Baf -B+t'-JP10c+f80 -R0/120/0/50 OR -JP10c+z -R0/120/30/80' -c
#gmt basemap -JP?+z -R0/120/3480/6371 -Baf -B+t'-JP10c+fp -R0/120/0/2891 OR -JP10c+z -R0/120/3480/6371' -c
#awk '{print $1, $2, $3}' xyz.txt| gmt xyz2grd -R0/100/0/100 -I1/5 -Gt.nc
awk '{print $1, 6340+$2, $3}' xyz2.txt| gmt xyz2grd -R0/3/6320/6370 -I1/25 -Gt.nc
gmt makecpt -Cmagma -T0/15/1 > seis.cpt
gmt grdimage t.nc -Baf -B+t -R0/3/6320/6370 -Cseis.cpt -JPa5z -Vl -png map
#gmt grdimage t.nc -B -JP15c+z6370 -png map
#gmt grdimage t.nc -B -R0/180/2200/6370 -JP15c+z -png map
#gmt grdimage t.nc -JP6i+f -Bafg -png map
#gmt grdimage t.nc -JP6i+f -Bafg -png map --PROJ_ELLIPSOID=sphere
#gmt subplot end
#gmt end
# ------------------------------------------------------------------------------------------------------------------- #
gmt grdinterpolate S362ANI_kmps.nc -E0/0/180/0+i0.25d+g+p -T25/2890/5 -Gslice.nc