#######################################################################################################################
# Description: Map showing the seismic networks used
#
# Konstantinos Michailos
# Lausanne
# January 2022
#######################################################################################################################
# ------------------------------------------------------------------------------------------------------------------- #
# Output name
out=bad_rfs.eps
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
gmt makecpt -C../files/grayC.cpt -T0/4000 -D+i > my_topo.cpt

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
gmt makecpt -Chot -T0/100  -I -M > seis.cpt
gmt makecpt -C../files/bamako.cpt -T0/100 -M -I > seis.cpt

echo Plot scale...
#gmt psscale -Dx2.5/8.9+o0/0i+w1.2i/0.08i+h+e -R -J -Cmy_topo.cpt -Bx500f250 -Bx+l"Topography (m)" -O -K  >> $out
#gmt psscale -Dx7.0/8.9+o0/0i+w1.2i/0.08i+h+e -R -J -Cseis.cpt -Bxa400f200 -Bx+l"Number of waveforms" -O -K  >> $out
gmt psscale -Dx7.0/8.9+o0/0i+w1.2i/0.08i+h+e -R -J -Cseis.cpt -Bxa20f10 -Bx+l"Percentage of bad RFs over all RFs" -O -K  >> $out
# -=================================================================================================================- #
# ------------------------------------------------------------------------------------------------------------------- #
# -=================================================================================================================- #
echo Plot seismic stations...
#awk '{print $3, $2, $4}' ../files/number_of_rf_calculated.txt | gmt psxy -i0,1,2 -Si.15 -R -J \
#-O -K -W.5p -Gred -t5 >> $out
#awk '{print $3, $2, $4}' ../files/number_of_rf_calculated.txt | gmt psxy -i0,1,2 -Si.2 -R -J \
#-O -K -W.5p,firebrick1  -t10 >> $out
#awk '{print $3, $2, $4}' ../files/bad_RFs.txt | gmt psxy -i0,1,2 -Si.2 -R -J \
#-O -K -W.5p -Cseis.cpt -t10 >> $out

awk '{print $3, $2, $4}' ../files/percentage_of_bad_rfs.txt | gmt psxy -i0,1,2 -Si.2 -R -J \
-O -K -W.5p -Cseis.cpt >> $out

# -=================================================================================================================- #
# ------------------------------------------------------------------------------------------------------------------- #
echo Create legend...
gmt set FONT_ANNOT_PRIMARY 7
gmt pslegend <<END -R -J -Dx3.5i/0.3i+w0i/0.0i/TC -C0.07i/0.1i -F+gwhite+pthin -P -O -K --FONT_ANNOT_PRIMARY=5p >> $out
G -0.05i
H 7 Seismic networks
D0.1i 0.5p
G .04i
S .04i i .11i white 0.8p 0.18i AASN, EASI, CIFALPS, PACASE
G .05i
END
# ------------------------------------------------------------------------------------------------------------------- #
gmt psxy -R -J -T -O >> $out
gmt psconvert -Tf -A $out
evince ${out%.*}.pdf