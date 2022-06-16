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

awk '{print $1, 6370-$2, $3}' xyz_smoothed_test.txt| gmt xyz2grd -R0/7/6280/6370 -I2.5m/2.5k -Gt_.nc -Vl

gmt grdview t_.nc -JPa30/2.5z -T+s+o0.01p,gray -By10+l"Depth (km)" -Bya5f5 -Bxa1f0.5+l"Distance (km)" -Cpol3.cpt -R0/6.3/6280/6370 -BWsNE
gmt psscale -Dx12.5/0.1+o0/0i+w1.5i/0.1i+h+e -Cpol3.cpt -Baf -Bx+l"Relative amplitude (%)"

gmt end show
