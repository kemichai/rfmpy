gmt begin polar png
gmt set FONT_TITLE 12p,9

gmt subplot begin 4x4 -Fs10c/10c -A+o-0.5c/-0.5c

gmt basemap -JP?+z -R0/120/30/80 -Baf -B+t'-JP10c+f80 -R0/120/0/50 OR -JP10c+z -R0/120/30/80' -c
gmt basemap -JP?+z -R0/120/3480/6371 -Baf -B+t'-JP10c+fp -R0/120/0/2891 OR -JP10c+z -R0/120/3480/6371' -c


gmt grd2xyz slice.nc | awk '{print $1, 6370-$2, $3}' | gmt xyz2grd -R0/180/3480/6345 -I1/5 -Gt.nc
gmt grdimage t.nc -B -R0/180/3480/6370 -JP15c+z -png map
gmt subplot end
gmt end show