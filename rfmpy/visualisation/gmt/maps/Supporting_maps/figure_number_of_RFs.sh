#######################################################################################################################
# Description: Map showing the seismic networks used
#
# Konstantinos Michailos
# Lausanne
# January 2022
#######################################################################################################################
# ------------------------------------------------------------------------------------------------------------------- #
# Output name
out=number_of_RFs.eps
# ------------------------------------------------------------------------------------------------------------------- #
# Define stuffz (area plotted, size of letters, etc)
gmt set FORMAT_GEO_MAP D
gmt set FONT_ANNOT_PRIMARY Helvetica
gmt set FONT_ANNOT_PRIMARY 8
gmt set FONT_LABEL Helvetica
gmt set LABEL_FONT_SIZE 7
gmt set MAP_FRAME_TYPE plain
# Directory containing topo grd file
topodir="/home/kmichailos/Desktop/topo"
# Map boundaries
north=52
south=41
east=25
west=0
proj='-JB10/45/25/45/5i'
# first two / / define the center of the map
#gmt coast -R110/140/20/35 -JB125/20/25/45/5i -Bag -Dl -Ggreen -Wthinnest -A250 -pdf GMT_albers

# ------------------------------------------------------------------------------------------------------------------- #
echo Make basemap...
gmt pscoast -W1/0.05 -Dl $proj -R$west/$east/$south/$north -K -B5WSen -P -X1 -Y10 > $out
# ------------------------------------------------------------------------------------------------------------------- #
echo Plot topo....
#gmt grdimage -R -J /home/kmichall/Desktop/topo/topo.0.20.40.55.3sec.grd -CFrance2.cpt -O -K >> $out
gmt grdimage -R -J $topodir/ETOPO1_Bed_g_gmt4.grd -Cmy_topo.cpt -O -K >> $out
# ------------------------------------------------------------------------------------------------------------------- #
gmt pscoast -W1/0.05 -Df -J -R -K -O -P -Sazure1 -N1/0.05p,black -L3.4/49.7/48/200+l+u >> $out


echo Plotting faults and stuff...
# 250km distance line from a smoothed 800m elevation contour
#gmt psxy -R -J d250km.dat -W1.5p,gray20 -O -K >> $out

# ---------
echo Create cpt...
gmt makecpt -Cviridis -T40/110/10  > seis.cpt
gmt makecpt -Chot -T0/400/50 -D+i -I > seis.cpt

echo Plot scale...
#gmt psscale -Dx2.5/8.9+o0/0i+w1.2i/0.08i+h+e -R -J -Cmy_topo.cpt -Bx500f250 -Bx+l"Topography (m)" -O -K  >> $out
#gmt psscale -Dx7.0/8.9+o0/0i+w1.2i/0.08i+h+e -R -J -Cseis.cpt -Bxa400f200 -Bx+l"Number of waveforms" -O -K  >> $out
gmt psscale -Dx7.0/8.9+o0/0i+w1.2i/0.08i+h+e -R -J -Cseis.cpt -Bxa100f50 -Bx+l"Number of RFs" -O -K  >> $out

#-Bxaf+l"topography" -By+lkm

echo Plot country names...
#gmt pstext -R -J -O -K  -F+f6p,Helvetica,black+jBL+a0 -Gwhite >> $out << END
#15.2 45.8 SLOVENIA
#15.7 45.5 CROATIA
#END
# -=================================================================================================================- #
# ------------------------------------------------------------------------------------------------------------------- #
# -=================================================================================================================- #

echo Plot initial 3D grid...
#awk '{print $1, $2}' files/initial_grid.txt |
#    gmt psxy -R -J -Sx.22 -W1.5p -Gred -O -K -t20 >> $out



echo Plot seismic stations...
awk '{print $3, $2, $4}' ../files/number_of_rf_calculated.txt | gmt psxy -i0,1,2 -Si.2 -R -J \
-O -K -W.5p -Cseis.cpt -t10 >> $out
# RF plot statiions
awk '{print $3, $2, $4}' ../files/rf_plot_stations.txt | gmt psxy -i0,1,2 -Sa.35 -R -J \
-O -K -W.5p -Cseis.cpt  >> $out
awk '{print $3, $2, $5}' ../files/rf_plot_stations.txt | gmt pstext -R -J -O -K -F+f5p,Helvetica,gray10  >> $out

# -=================================================================================================================- #
# ------------------------------------------------------------------------------------------------------------------- #
echo Create legend...
gmt set FONT_ANNOT_PRIMARY 7
gmt pslegend <<END -R -J -Dx3.5i/0.3i+w0i/0.0i/TC -C0.07i/0.1i -F+gwhite+pthin -P -O -K --FONT_ANNOT_PRIMARY=5p >> $out
G -0.05i
H 7 Seismic networks
D0.1i 0.5p
G .04i
S .04i i .11i red 0.3p 0.18i AASN, EASI, CIFALPS, PACASE
G .05i
END

#echo Plot cross section lines
#start_lon='8'
#start_lat='45.5'
#end_lon='15'
#end_lat='50'
#
#gmt psxy << END -R -J -O -W3,dodgerblue -K >> $out
#$start_lon $start_lat
#$end_lon $end_lat
#END
#gmt pstext -R -J -D0/0.23 -O -K -F+f12p,Helvetica,gray10 -TO -Gwhite -W0.1 >> $out << END
#$start_lon $start_lat A
#END
#gmt pstext -R -J -D0/0.23 -O -K -F+f12p,Helvetica,gray10 -TO -Gwhite -W0.1 >> $out << END
#$end_lon $end_lat B
#END
echo Make basemap...
#gmt pscoast -W1/0.05 -Dl $proj -R$west/$east/$south/$north -K -O -Y11 -X-0.01 -B0.5WSen -P >> $out
#gmt psbasemap -R -J -O -K -DjTR+w1.7i+o0.2i/-1.8i+stmp -F+gwhite+p1p >> $out
#gmt pscoast -Rg -JE10/40/3.0i -Dc -A -Bg -Glightgray -Wthinnest -O -K -X0 -Y-8.5 >> $out
#gmt coast -Rg -JE10/40/110/4.4i -Bg  -BWsNE -Dc -A10000 -Glightgray -Wthinnest

#gmt makecpt -Chot -T0/200/50 -H > n.cpt

#awk '{print $1, $2, $3/1000}' ../tele_events_sample.txt | gmt psxy -R -J -i0,1,2s0.2 -Sc.05 -Wthinnest -O -K -t5 -Cn.cpt >> $out
##	# Fancy line
#gmt psxy -R -J -O -K -B -W2.5p,red >> $out << END
#0	45
#20 45
#20	50
#0	  50
#0 45
#END
#gmt psscale -DjLM+w5.0c+jRM+o1.5c+e -Cn.cpt -Bxa50+lDEPTH -By+lkm >> $out
#


#gmt psxy << END -R -J -O -W1.,black,- -K>> $out
#$start_lon $start_lat
#$end_lon $end_lat
#END
#	gmt plot -R0/19/0/25 -Jx1c -B0 -W2p -X-6c -Y-13.5c <<- EOF
#	3	13
#	20 13
#	20 33
#	3  33
#	3  13
#	EOF

# ------------------------------------------------------------------------------------------------------------------- #
gmt psxy -R -J -T -O >> $out
gmt psconvert -Tf -A $out
gmt psconvert -Tg -A+r $out
evince ${out%.*}.pdf
# Delete eps file as it is large
rm ${out%.*}.eps
