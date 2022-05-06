gmt begin polar png
gmt set FONT_TITLE 12p,9

gmt subplot begin 4x4 -Fs10c/10c -A+o-0.5c/-0.5c

gmt basemap -JP?+z -R0/100/30/80 -Baf -B+t'-JP10c+f80 -R0/120/0/50 OR -JP10c+z -R0/120/30/80' -c
gmt basemap -JP?+z -R0/120/3480/6371 -Baf -B+t'-JP10c+fp -R0/120/0/2891 OR -JP10c+z -R0/120/3480/6371' -c
awk '{print $1, $2, $3}' xyz.txt| gmt xyz2grd -R0/100/0/100 -I1/5 -Gt.nc
gmt grdimage t.nc -Baf -B+t -R0/100/0/50 -JP10c+z -png map
#gmt grdimage t.nc -B -JP15c+z6370 -png map
#gmt grdimage t.nc -B -R0/180/2200/6370 -JP15c+z -png map
#gmt grdimage t.nc -JP6i+f -Bafg -png map
#gmt grdimage t.nc -JP6i+f -Bafg -png map --PROJ_ELLIPSOID=sphere
gmt subplot end
gmt end show