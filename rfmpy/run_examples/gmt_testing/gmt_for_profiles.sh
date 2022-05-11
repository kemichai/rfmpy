# ------------------------------------------------------------------------------------------------------------------- #
# Define stuffz (area plotted, size of letters, etc)
#gmt begin polar png

gmt set FONT_TITLE 12p,9

#gmt subplot begin 4x4 -Fs10c/10c -A+o-0.5c/-0.5c

#gmt basemap -JP?+z -R0/100/30/80 -Baf -B+t'-JP10c+f80 -R0/120/0/50 OR -JP10c+z -R0/120/30/80' -c
#gmt basemap -JP?+z -R0/120/3480/6371 -Baf -B+t'-JP10c+fp -R0/120/0/2891 OR -JP10c+z -R0/120/3480/6371' -c
#awk '{print $1, $2, $3}' xyz.txt| gmt xyz2grd -R0/100/0/100 -I1/5 -Gt.nc
#awk '{print $1, 6370-$2, $3}' xyz_example.txt| gmt xyz2grd -R0/3/6320/6370 -I30m/10k -Gt.nc -Vl

#awk '{print $1, 6370-$2, $3}' xyz_example.txt| gmt surface -R0/5/6280/6370 -I5m/10k -Gt.nc
awk '{print $1, 6370-$2, $3}' xyz_example.txt| gmt surface -R0/6/6270/6380 -I1m/1k -Gt.nc
gmt makecpt -Chot -T0/20/1 > seis.cpt
gmt makecpt -Cpolar -T-0.03/0.03/0.005 -D+i > pol.cpt
#gmt grdimage t.nc -Baf -B+t -R0/5/6280/6370  -JPa30z  -png cross_section
#gmt grdimage t.nc -Baf -B+t -R0/5/6280/6370 -Cpol.cpt -JPa30z -png cross_section

awk '{print $1, 6370-$2, $3}' xyz.txt| gmt surface -R0/3/6270/6380 -I0.5m/0.1k -Gt.nc -T0.7
gmt grdimage t.nc -Baf -B+t -R0/3/6280/6370 -Cpol.cpt -JPa30z  -png example



#gmt grdimage t.nc -B -JP15c+z6370 -png map
#gmt grdimage t.nc -B -R0/180/2200/6370 -JP15c+z -png map
#gmt grdimage t.nc -JP6i+f -Bafg -png map
#gmt grdimage t.nc -JP6i+f -Bafg -png map --PROJ_ELLIPSOID=sphere
#gmt subplot end
#gmt end
# ------------------------------------------------------------------------------------------------------------------- #
#gmt grdinterpolate S362ANI_kmps.nc -E0/0/180/0+i0.25d+g+p -T25/2890/5 -Gslice.nc
