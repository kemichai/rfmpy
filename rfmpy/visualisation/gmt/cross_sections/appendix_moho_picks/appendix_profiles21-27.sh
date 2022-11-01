#######################################################################################################################
# Description: Cross-section of migrated receiver function amplitudes
#
# KM
# Chavannes
# June 2022
#######################################################################################################################
# Create .cpt file
gmt makecpt -Cpolar -T-0.05/0.05/0.005 -D+i > pol.cpt
gmt makecpt -Cvik.cpt -T-0.11/0.11/0.002 > pol_vik.cpt
gmt begin profiles_long pdf
gmt set FONT_TITLE 12p,9
gmt set FORMAT_GEO_MAP D
gmt set FORMAT_GEO_MAP D
gmt set FONT_ANNOT_PRIMARY Helvetica
gmt set FONT_ANNOT_PRIMARY 10
gmt set FONT_LABEL 12
gmt set MAP_FRAME_TYPE plain


awk '{print $1, 6370-$2, $3}' Cross-section_0.txt| gmt xyz2grd -R0/15.4/6290/6370 -I2.m/2k -Gt_.nc -Vl
gmt grdview t_.nc -JPa30/2.5z -T+s0.01p,gray -By10+l"Depth (km)" -Bya5f5 -Bxa1f0.5+l"Distance (km)" -Cpol_vik.cpt -R0/15.3/6290/6370 -BWsNE+t"Cross section 24"
gmt psscale -Dx12.5/-.4+o0/0i+w1.5i/0.1i+h+e -Cpol_vik.cpt -Baf -Bx+l"Relative amplitude (%)"
awk '{print($2-43, 6370-$3)}' picks_cc0.txt | gmt psxy -Gwhite -Sd.25 -W1.2p,black -l"Moho picks"
awk '{print($2-43, 6370-$3)}' unc_picks_cc0.txt | gmt psxy -Gdimgrey -Sd.25 -W1.2p,black -l"Unc. Moho picks"
gmt legend -DjTR+o1.1c -F+pthick+gwhite --FONT_ANNOT_PRIMARY=11p
#
#
#awk '{print $1, 6370-$2, $3}' Cross-section_0.txt| gmt xyz2grd -R0/15.4/6290/6370 -I2.m/2k -Gt_.nc -Vl
#gmt grdview t_.nc -JPa50/7z -T+s0.01p,gray -By10+l"Depth (km)" -Bya5f5 -Bxa1f0.5+l"Distance (km)" -Cpol_vik.cpt -R0/15.3/6290/6370 -BWsNE+t"Cross section 23" -Yh+0.8c
#awk '{print($2-43, 6370-$3)}' picks_cc0.txt | gmt psxy -Gwhite -Sd.25 -W1.2p,black
#awk '{print($2-43, 6370-$3)}' unc_picks_cc0.txt | gmt psxy -Gdimgrey -Sd.25 -W1.2p,black
#
#awk '{print $1, 6370-$2, $3}' Cross-section_0.txt| gmt xyz2grd -R0/15.4/6290/6370 -I2.m/2k -Gt_.nc -Vl
#gmt grdview t_.nc -JPa50/7z -T+s0.01p,gray -By10+l"Depth (km)" -Bya5f5 -Bxa1f0.5+l"Distance (km)" -Cpol_vik.cpt -R0/15.3/6290/6370 -BWsNE+t"Cross section 22" -Yh+0.8c
#awk '{print($2-43, 6370-$3)}' picks_cc0.txt | gmt psxy -Gwhite -Sd.25 -W1.2p,black
#awk '{print($2-43, 6370-$3)}' unc_picks_cc0.txt | gmt psxy -Gdimgrey -Sd.25 -W1.2p,black
#

awk '{print $1, 6370-$2, $3}' Cross-section_21.txt| gmt xyz2grd -R0/15.4/6290/6370 -I2.m/2k -Gt_.nc -Vl
gmt grdview t_.nc -JPa30/2.5z -T+s0.01p,gray -By10+l"Depth (km)" -Bya5f5 -Bxa1f0.5+l"Distance (km)" -Cpol_vik.cpt -R0/15.3/6290/6370 -BWsNE+t"Cross section 21" -Yh+0.8c
awk '{print($1-2, 6370-$3)}' picks_cc21.txt | gmt psxy -Gwhite -Sd.25 -W1.2p,black
awk '{print($1-2, 6370-$3)}' unc_picks_cc21.txt | gmt psxy -Gdimgrey -Sd.25 -W1.2p,black


gmt end show
