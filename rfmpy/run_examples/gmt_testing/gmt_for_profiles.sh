# ------------------------------------------------------------------------------------------------------------------- #
# Define stuffz (area plotted, size of letters, etc)
#gmt begin polar png

gmt set FONT_TITLE 12p,9
gmt makecpt -Cturbo -T0/20/1 > seis.cpt
gmt makecpt -Cpolar -T-0.03/0.03/0.005 -D+i > pol.cpt
#gmt subplot begin 4x4 -Fs10c/10c -A+o-0.5c/-0.5c


# Generic example
#gmt basemap -JP?+z -R0/100/30/80 -Baf -B+t'-JP10c+f80 -R0/120/0/50 OR -JP10c+z -R0/120/30/80' -c
#gmt basemap -JP?+z -R0/120/3480/6371 -Baf -B+t'-JP10c+fp -R0/120/0/2891 OR -JP10c+z -R0/120/3480/6371' -c
#awk '{print $1, $2, $3}' xyz.txt| gmt xyz2grd -R0/100/0/100 -I1/5 -Gt.nc
#awk '{print $1, 6370-$2, $3}' xyz_example.txt| gmt xyz2grd -R0/6/6280/6370 -I10m/20k -Gt.nc -Vl

#awk '{print $1, 6370-$2, $3}' xyz_example.txt| gmt surface -R0/5/6280/6370 -I5m/10k -Gt.nc
#awk '{print $1, 6370-$2, $3}' xyz_example.txt| gmt surface -R0/6/6270/6380 -I10m/10k -Gt.nc
#gmt grdimage t.nc -Baf -B+t -R0/5/6280/6370  -JPa30z  -png cross_section
#gmt grdview t.nc -JPa30z  -Baf -B+t -Cseis.cpt -S5 -T+o+s -R0/6/6280/6370 -Wc -png grd_view
#gmt grdimage t.nc -Baf -B+t -R0/5/6280/6370 -Cpol.cpt -JPa30z -png cross_section


# RF migration example
#awk '{print $1, 6370-$2, $3}' xyz.txt| gmt surface -R0/3/6280/6370 -I5m/5k -Gt_.nc -T0.7
awk '{print $1, 6370-$2, $3}' xyz.txt| gmt xyz2grd -R0/3/6280/6370 -I2.5m/2.5k -Gt_.nc -Vl
#gmt grdimage t.nc -Baf -B+t -R0/3/6280/6370 -Cpol.cpt -JPa30z  -png example
gmt grdview t_.nc -JPa30z -T+o+s -Baf -B+t -Cpol.cpt -R0/3/6280/6370  -png example
Or use -T+o instead of -Qsm


