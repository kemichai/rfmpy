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
gmt begin profiles21-24 pdf
gmt set FONT_TITLE 12p,9
gmt set FORMAT_GEO_MAP D
gmt set FORMAT_GEO_MAP D
gmt set FONT_ANNOT_PRIMARY Helvetica
gmt set FONT_ANNOT_PRIMARY 10
gmt set FONT_LABEL 12
gmt set MAP_FRAME_TYPE plain


awk '{print $1, 6370-$2, $3}' cc_files/Cross-section_24.txt| gmt xyz2grd -R0/15.4/6290/6370 -I2.m/2k -Gt_.nc -Vl
gmt grdview t_.nc -JPa50/7z -T+s0.01p,gray -By10+l"Depth (km)" -Bya5f5 -Bxa1f0.5+l"Distance (km)" -Cpol_vik.cpt -R0/15.3/6290/6370 -BWsNE+t"Cross section 24"
gmt psscale -Dx23./-.4+o0/0i+w1.5i/0.1i+h+e -Cpol_vik.cpt -Baf -Bx+l"Relative amplitude (%)"
awk '{print($4, 6370-$3)}' pick_files/picks_cc24.txt | gmt psxy -Gwhite -Sd.25 -W1.2p,black -l"Moho picks"
awk '{print($4, 6370-$3)}' pick_files/unc_picks_cc24.txt | gmt psxy -Gdimgrey -Sd.25 -W1.2p,black -l"Unc. Moho picks"
gmt legend -DjTR+o2.c -F+pthick+gwhite --FONT_ANNOT_PRIMARY=11p
#
#
awk '{print $1, 6370-$2, $3}' cc_files/Cross-section_23.txt| gmt xyz2grd -R0/15.4/6290/6370 -I2.m/2k -Gt_.nc -Vl
gmt grdview t_.nc -JPa50/7z -T+s0.01p,gray -By10+l"Depth (km)" -Bya5f5 -Bxa1f0.5+l"Distance (km)" -Cpol_vik.cpt -R0/15.3/6290/6370 -BWsNE+t"Cross section 23" -Yh+0.2c
awk '{print($4, 6370-$3)}' pick_files/picks_cc23.txt| gmt psxy -Gwhite -Sd.25 -W1.2p,black
awk '{print($4, 6370-$3)}' pick_files/unc_picks_cc23.txt | gmt psxy -Gdimgrey -Sd.25 -W1.2p,black
#
awk '{print $1, 6370-$2, $3}' cc_files/Cross-section_22.txt| gmt xyz2grd -R0/15.4/6290/6370 -I2.m/2k -Gt_.nc -Vl
gmt grdview t_.nc -JPa50/7z -T+s0.01p,gray -By10+l"Depth (km)" -Bya5f5 -Bxa1f0.5+l"Distance (km)" -Cpol_vik.cpt -R0/15.3/6290/6370 -BWsNE+t"Cross section 22" -Yh+0.2c
awk '{print($4, 6370-$3)}' pick_files/picks_cc22.txt| gmt psxy -Gwhite -Sd.25 -W1.2p,black
awk '{print($4, 6370-$3)}' pick_files/unc_picks_cc22.txt | gmt psxy -Gdimgrey -Sd.25 -W1.2p,black


awk '{print $1, 6370-$2, $3}' cc_files/Cross-section_21.txt| gmt xyz2grd -R0/15.4/6290/6370 -I2.m/2k -Gt_.nc -Vl
gmt grdview t_.nc -JPa50/7z -T+s0.01p,gray -By10+l"Depth (km)" -Bya5f5 -Bxa1f0.5+l"Distance (km)" -Cpol_vik.cpt -R0/15.3/6290/6370 -BWsNE+t"Cross section 21" -Yh+0.2c
awk '{print($4, 6370-$3)}' pick_files/picks_cc21.txt | gmt psxy -Gwhite -Sd.25 -W1.2p,black
awk '{print($4, 6370-$3)}' pick_files/unc_picks_cc21.txt | gmt psxy -Gdimgrey -Sd.25 -W1.2p,black


gmt end show
