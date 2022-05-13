#gmt makecpt -Cpolar -T-0.03/0.03/0.005 -D+i > pol.cpt

gmt begin test png
gmt set FONT_TITLE 12p,9
gmt set FORMAT_GEO_MAP D
gmt set FORMAT_GEO_MAP D
gmt set FONT_ANNOT_PRIMARY Helvetica
gmt set FONT_ANNOT_PRIMARY 10
gmt set FONT_LABEL 12
gmt set MAP_FRAME_TYPE plain

#gmt basemap -R0/6.3/6280/6370 -Baf  -BWSn -B+t -JPa50z

awk '{print $1, 6370-$2, $3}' xyz_egu_smoothed.txt| gmt xyz2grd -R0/7/6280/6370 -I2.5m/2.5k -Gt_.nc -Vl

gmt grdview t_.nc -JPa30/2.5z -T+s+o0.01p,gray -By10+l"Depth (km)" -Bya5f5 -Bxa1f0.5+l"Distance (km)" -Cpol3.cpt -R0/6.3/6280/6370 -BWsNE
gmt psscale -Dx12.5/0.1+o0/0i+w1.5i/0.1i+h+e -Cpol3.cpt -Baf -Bx+l"Relative amplitude (%)"

gmt end show

#gmt begin polar png
#gmt set FONT_TITLE 12p,9
#
#gmt subplot begin 4x4 -Fs10c/10c -A+o-0.5c/-0.5c
#gmt basemap -JP? -R0/120/0/60 -Baf -B+t'-JP10c -R0/120/0/60' -c
#gmt basemap -JP? -R0/120/30/80 -Baf -B+t'-JP10c -R0/120/30/80' -c
#gmt basemap -JP? -R0/120/0/2891 -Baf -B+t'-JP10c -R0/120/0/2891' -c
#gmt basemap -JP? -R0/120/3480/6371 -Baf -B+t'-JP10c -R0/120/3480/6371 OR -JP10c+fp+z -R0/120/0/2891' -c
#
#gmt basemap -JP?+z -R0/120/0/60 -Baf -B+t'-JP10c+f -R0/120/0/60' -c
#gmt basemap -JP?+z -R0/120/30/80 -Baf -B+t'-JP10c+f80 -R0/120/0/50 OR -JP10c+z -R0/120/30/80' -c
#gmt basemap -JP?+z -R0/120/0/2891 -Baf -B+t'-JP10c+f -R0/120/0/2891' -c
#gmt basemap -JP?+z -R0/120/3480/6371 -Baf -B+t'-JP10c+fp -R0/120/0/2891 OR -JP10c+z -R0/120/3480/6371' -c
#
#gmt basemap -JP?+r -R0/120/0/60 -Baf -B+t'-JP10c+fr -R0/120/0/60' -c
#gmt basemap -JP?+r -R0/120/30/80 -Baf -B+t'-JP10c+fr -R0/120/30/80' -c
#
#gmt subplot end
#gmt end show