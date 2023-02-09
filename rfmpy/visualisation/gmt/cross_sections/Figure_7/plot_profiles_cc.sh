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
gmt begin cc pdf
gmt set FONT_TITLE 12p,9
gmt set FORMAT_GEO_MAP D
gmt set FORMAT_GEO_MAP D
gmt set FONT_ANNOT_PRIMARY Helvetica
gmt set FONT_ANNOT_PRIMARY 10
gmt set FONT_LABEL 12
gmt set MAP_FRAME_TYPE plain

#gmt basemap -R0/6.3/6280/6370 -Baf  -BWSn -B+t -JPa50z

# Read xyz file created with python codes...
# x is distance along the profile in degrees
# y is depth in km
# z is the amplitude (smoothed using a Gaussian filter in this case)
awk '{print $1, 6370-$2, $3}' C_C.txt| gmt xyz2grd -R0/8.6/6290/6370 -I2.m/2k -Gt_.nc -Vl

# Plot
gmt grdview t_.nc -JPa30/2.5z -T+s0.01p,gray -By10+l"Depth (km)" -Bya5f5 -Bxa1f0.5+l"Distance (km)" -Cpol_vik.cpt -R0/8.5/6290/6370 -BWsNE
# to remove the frame take out +o
gmt psscale -Dx12.5/-.4+o0/0i+w1.5i/0.1i+h+e -Cpol_vik.cpt -Baf -Bx+l"Relative amplitude (%)"

#start_lon_C='11'
#start_lat_C='45.5'
#end_lon_C='22'
#end_lat_C='50'
# Attempt to compare to Grad 2007
#gmt grd2xyz ../Figure_EASI/Europe_moho_depth_2007.grd > grad.xyz
awk '{print $1, $2, $3}' ../Figure_EASI/grad.xyz | gmt project -C11/45.5 -E22/50 -W-5/5 -Q -Fpz -S > grad.dat
awk '{print($1/111.11, 6370-$2)}' grad.dat | gmt psxy -W1.1,black -t0
# Attempt to compare to SPADA
awk '{print $1, $2, $3}' ../Figure_EASI/Spada_moho.dat | gmt project -C11/45.5 -E22/50 -W-5/5  -Q -Fpz -S > spada_B.dat
awk '{print($1/111, 6370-$2, $3)}' spada_B.dat | gmt psxy -W1.,gray20,-- -t0

awk '{print $1, $2, $3}' profile_cc_certain_picks.xyz | gmt project -C11/45.5 -E22/50 -W-5/5 -Q -Fpz -S > m.dat
awk '{print($1/111.11, 6370-$2)}' m.dat | gmt psxy -Gwhite -Sd.25 -W1.2p,black -l"Moho picks"

awk '{print $1, $2, $3}' profile_cc_uncertain_picks.xyz | gmt project -C11/45.5 -E22/50 -W-10/10 -Q -Fpz -S > m.dat
awk '{print($1/111.11, 6370-$2)}' m.dat | gmt psxy -Gdimgrey -Sd.25 -W1.2p,black -l"Unc. Moho picks"

gmt end show
