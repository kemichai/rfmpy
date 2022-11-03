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
gmt begin profile0-6 pdf
gmt set FONT_TITLE 12p,9
gmt set FORMAT_GEO_MAP D
gmt set FORMAT_GEO_MAP D
gmt set FONT_ANNOT_PRIMARY Helvetica
gmt set FONT_ANNOT_PRIMARY 10
gmt set FONT_LABEL 12
gmt set MAP_FRAME_TYPE plain


awk '{print $1, 6370-$2, $3}' cc_files/Cross-section_6.txt| gmt xyz2grd -R0/9.1/6290/6370 -I2.m/2k -Gt_.nc -Vl
gmt grdview t_.nc -JPa30/3.5z -T+s0.01p,gray -By10+l"Depth (km)" -Bya5f5 -Bxa1f0.5+l"Distance (km)" -Cpol_vik.cpt -R0/9/6290/6370 -BWsNE+t"Cross section 6"
gmt psscale -Dx12.5/-.4+o0/0i+w1.5i/0.1i+h+e -Cpol_vik.cpt -Baf -Bx+l"Relative amplitude (%)"
awk '{print($2-43, 6370-$3)}' pick_files/moho_depths_Cross-section_6.txt | gmt psxy -Gwhite -Sd.25 -W1.2p,black -l"Moho picks"
awk '{print($2-43, 6370-$3)}' pick_files/unc_moho_depths_Cross-section_6.txt | gmt psxy -Gdimgrey -Sd.25 -W1.2p,black -l"Unc. Moho picks"
gmt legend -DjTR+o1.1c -F+pthick+gwhite --FONT_ANNOT_PRIMARY=11p


awk '{print $1, 6370-$2, $3}' cc_files/Cross-section_5.txt| gmt xyz2grd -R0/9.1/6290/6370 -I2.m/2k -Gt_.nc -Vl
gmt grdview t_.nc -JPa30/3.5z -T+s0.01p,gray -By10+l"Depth (km)" -Bya5f5 -Bxa1f0.5+l"Distance (km)" -Cpol_vik.cpt -R0/9/6290/6370 -BWsNE+t"Cross section 5" -Yh+0.8c
awk '{print($2-43, 6370-$3)}' pick_files/moho_depths_Cross-section_5.txt | gmt psxy -Gwhite -Sd.25 -W1.2p,black
awk '{print($2-43, 6370-$3)}' pick_files/unc_moho_depths_Cross-section_5.txt | gmt psxy -Gdimgrey -Sd.25 -W1.2p,black

awk '{print $1, 6370-$2, $3}' cc_files/Cross-section_4.txt| gmt xyz2grd -R0/9.1/6290/6370 -I2.m/2k -Gt_.nc -Vl
gmt grdview t_.nc -JPa30/3.5z -T+s0.01p,gray -By10+l"Depth (km)" -Bya5f5 -Bxa1f0.5+l"Distance (km)" -Cpol_vik.cpt -R0/9/6290/6370 -BWsNE+t"Cross section 4" -Yh+0.8c
awk '{print($2-43, 6370-$3)}' pick_files/moho_depths_Cross-section_4.txt | gmt psxy -Gwhite -Sd.25 -W1.2p,black
awk '{print($2-43, 6370-$3)}' pick_files/unc_moho_depths_Cross-section_4.txt | gmt psxy -Gdimgrey -Sd.25 -W1.2p,black


awk '{print $1, 6370-$2, $3}' cc_files/Cross-section_3.txt| gmt xyz2grd -R0/9.1/6290/6370 -I2.m/2k -Gt_.nc -Vl
gmt grdview t_.nc -JPa30/3.5z -T+s0.01p,gray -By10+l"Depth (km)" -Bya5f5 -Bxa1f0.5+l"Distance (km)" -Cpol_vik.cpt -R0/9/6290/6370 -BWsNE+t"Cross section 3" -Yh+0.8c
awk '{print($2-43, 6370-$3)}' pick_files/moho_depths_Cross-section_3.txt | gmt psxy -Gwhite -Sd.25 -W1.2p,black
awk '{print($2-43, 6370-$3)}' pick_files/unc_moho_depths_Cross-section_3.txt | gmt psxy -Gdimgrey -Sd.25 -W1.2p,black

awk '{print $1, 6370-$2, $3}' cc_files/Cross-section_2.txt| gmt xyz2grd -R0/9.1/6290/6370 -I2.m/2k -Gt_.nc -Vl
gmt grdview t_.nc -JPa30/3.5z -T+s0.01p,gray -By10+l"Depth (km)" -Bya5f5 -Bxa1f0.5+l"Distance (km)" -Cpol_vik.cpt -R0/9/6290/6370 -BWsNE+t"Cross section 2" -Yh+0.8c
awk '{print($2-43, 6370-$3)}' pick_files/moho_depths_Cross-section_2.txt | gmt psxy -Gwhite -Sd.25 -W1.2p,black
awk '{print($2-43, 6370-$3)}' pick_files/unc_moho_depths_Cross-section_2.txt | gmt psxy -Gdimgrey -Sd.25 -W1.2p,black

awk '{print $1, 6370-$2, $3}' cc_files/Cross-section_1.txt| gmt xyz2grd -R0/9.1/6290/6370 -I2.m/2k -Gt_.nc -Vl
gmt grdview t_.nc -JPa30/3.5z -T+s0.01p,gray -By10+l"Depth (km)" -Bya5f5 -Bxa1f0.5+l"Distance (km)" -Cpol_vik.cpt -R0/9/6290/6370 -BWsNE+t"Cross section 1" -Yh+0.8c
awk '{print($2-43, 6370-$3)}' pick_files/moho_depths_Cross-section_1.txt | gmt psxy -Gwhite -Sd.25 -W1.2p,black
awk '{print($2-43, 6370-$3)}' pick_files/unc_moho_depths_Cross-section_1.txt | gmt psxy -Gdimgrey -Sd.25 -W1.2p,black


awk '{print $1, 6370-$2, $3}' cc_files/Cross-section_0.txt| gmt xyz2grd -R0/9.1/6290/6370 -I2.m/2k -Gt_.nc -Vl
gmt grdview t_.nc -JPa30/3.5z -T+s0.01p,gray -By10+l"Depth (km)" -Bya5f5 -Bxa1f0.5+l"Distance (km)" -Cpol_vik.cpt -R0/9/6290/6370 -BWsNE+t"Cross section 0" -Yh+0.8c
awk '{print($2-43, 6370-$3)}' pick_files/moho_depths_Cross-section_0.txt | gmt psxy -Gwhite -Sd.25 -W1.2p,black
awk '{print($2-43, 6370-$3)}' pick_files/unc_moho_depths_Cross-section_0.txt | gmt psxy -Gdimgrey -Sd.25 -W1.2p,black


gmt end show
