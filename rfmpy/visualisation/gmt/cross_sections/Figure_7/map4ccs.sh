#######################################################################################################################
# Description: Map showing the seismic networks used
#
# Konstantinos Michailos
# Lausanne
# January 2022
#######################################################################################################################
# ------------------------------------------------------------------------------------------------------------------- #
# Output name
out=map.eps
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

gmt makecpt -C../../maps/files/bamako.cpt -T0/400/50 -D+i -I > seis.cpt
gmt makecpt -C../../maps/files/grayC.cpt -T0/4000 -D+i > my_topo.cpt


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
#echo Create cpt...
#gmt makecpt -Cviridis -T40/110/10  > seis.cpt
#gmt makecpt -Chot -T0/400/50 -D+i -I > seis.cpt

echo Plot scale...
#gmt psscale -Dx2.5/8.9+o0/0i+w1.2i/0.08i+h+e -R -J -Cmy_topo.cpt -Bx500f250 -Bx+l"Topography (m)" -O -K  >> $out
#gmt psscale -Dx7.0/8.9+o0/0i+w1.2i/0.08i+h+e -R -J -Cseis.cpt -Bxa400f200 -Bx+l"Number of waveforms" -O -K  >> $out
#gmt psscale -Dx7.0/8.9+o0/0i+w1.2i/0.08i+h+e -R -J -Cseis.cpt -Bxa100f50 -Bx+l"Number of RFs" -O -K  >> $out

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
awk '{print $3, $2, $4}' ../../maps/files/number_of_rf_calculated.txt | gmt psxy -i0,1,2 -Si.2 -R -J \
-O -K -W.7p,dodgerblue  -t10 >> $out
awk '{print $3, $2, $4}' ../../maps/files/number_of_waveforms_EASI.txt | gmt psxy -i0,1,2 -St.2 -R -J \
-O -K -W.7p,dodgerblue  -t10 >> $out
awk '{print $3, $2, $4}' ../../maps/files/number_of_rf_calculated_PACASE.txt | gmt psxy -i0,1,2 -Sd.2 -R -J \
-O -K -W.7p,dodgerblue  -t10 >> $out
awk '{print $3, $2, $4}' ../../maps/files/number_of_rf_calculated_CIFALPS.txt | gmt psxy -i0,1,2 -Ss.2 -R -J \
-O -K -W.7p,dodgerblue  -t10 >> $out


echo Plot country names...
gmt pstext -R -J -O -K  -F+f6p,Helvetica,black+jBL+a0 -Gwhite >> $out << END
13 44 Adria
5 50.8 Europe
10 41.5 Liguria
END



# -=================================================================================================================- #
# ------------------------------------------------------------------------------------------------------------------- #
echo Create legend...
gmt pslegend <<END -R -J -Dx3.7i/0.3i+w0i/0.0i/TC -C0.1i/0.1i -F+gwhite+pthin -P -O -K --FONT_ANNOT_PRIMARY=6p >> $out
G -0.05i
H 7 Seismic networks
D0.1i 0.5p
G .04i
S .04i i .09i white 0.5p 0.18i AASN
G .05i
S .04i s .09i white 0.5p 0.18i CIFALPS
G .05i
S .04i t .09i white 0.5p 0.18i EASI
G .05i
S .04i d .09i white 0.5p 0.18i PACASE
END

# EASI
#-C13.33/50.6 -E13.33/45.6
start_lon='13.33'
start_lat='50.6'
end_lon='13.33'
end_lat='45.6'
gmt psxy << END -R -J -O -W0.8,black -K>> $out
$start_lon $start_lat
$end_lon $end_lat
END
gmt psxy << END -R -J -O -W0.15c,black,1_18 -K>> $out
$start_lon $start_lat
$end_lon $end_lat
END
gmt psxy << END -R -J -O -Sc.13 -Gblack -W0.8,black -K>> $out
$start_lon $start_lat
$end_lon $end_lat
END
gmt pstext -R -J -D0/0.23 -O -K -F+f10p,Helvetica,gray10+jB -TO -Gwhite -W0.1 >> $out << END
$start_lon $start_lat D
END
gmt pstext -R -J -D0/0.23 -O -K -F+f10p,Helvetica,gray10+jB -TO -Gwhite -W0.1 >> $out << END
$end_lon $end_lat D'
END

# Cross-section A-A' start and end points
start_lon_A='3'
start_lat_A='44.1'
end_lon_A='9'
end_lat_A='44.8'
gmt psxy << END -R -J -O -W0.8,black -K>> $out
$start_lon_A $start_lat_A
$end_lon_A $end_lat_A
END
gmt psxy << END -R -J -O -W0.15c,black,1_18 -K>> $out
$start_lon_A $start_lat_A
$end_lon_A $end_lat_A
END
gmt psxy << END -R -J -O -Sc.13 -Gblack -W0.8,black -K>> $out
$start_lon_A $start_lat_A
$end_lon_A $end_lat_A
END
gmt pstext -R -J -D0/0.23 -O -K -F+f10p,Helvetica,gray10+jB -TO -Gwhite -W0.1 >> $out << END
$start_lon_A $start_lat_A A
END
gmt pstext -R -J -D0/0.23 -O -K -F+f10p,Helvetica,gray10+jB -TO -Gwhite -W0.1 >> $out << END
$end_lon_A $end_lat_A A'
END
# Cross-section B-B' start and end points
start_lon_B='6'
start_lat_B='49'
end_lon_B='11.5'
end_lat_B='44'
gmt psxy << END -R -J -O -W0.8,black -K>> $out
$start_lon_B $start_lat_B
$end_lon_B $end_lat_B
END
gmt psxy << END -R -J -O -W0.15c,black,1_18 -K>> $out
$start_lon_B $start_lat_B
$end_lon_B $end_lat_B
END
gmt psxy << END -R -J -O -Sc.13 -Gblack -W0.8,black -K>> $out
$start_lon_B $start_lat_B
$end_lon_B $end_lat_B
END
gmt pstext -R -J -D0/0.23 -O -K -F+f10p,Helvetica,gray10+jB -TO -Gwhite -W0.1 >> $out << END
$start_lon_B $start_lat_B B
END
gmt pstext -R -J -D0/0.23 -O -K -F+f10p,Helvetica,gray10+jB -TO -Gwhite -W0.1 >> $out << END
$end_lon_B $end_lat_B B'
END
# Cross-section C-C' start and end points
start_lon_C='11'
start_lat_C='45.5'
end_lon_C='22'
end_lat_C='50'
echo aaa
gmt psxy << END -R -J -O -W0.8,black,solid -K>> $out
$start_lon_C $start_lat_C
$end_lon_C $end_lat_C
END
gmt psxy << END -R -J -O -W0.15c,black,1_18 -K>> $out
$start_lon_C $start_lat_C
$end_lon_C $end_lat_C
END
gmt psxy << END -R -J -O -Sc.13 -Gblack -W0.8,black -K>> $out
$start_lon_C $start_lat_C
$end_lon_C $end_lat_C
END
gmt pstext -R -J -D0/0.23 -O -K -F+f10p,Helvetica,gray10+jB -TO -Gwhite -W0.1 >> $out << END
$start_lon_C $start_lat_C C
END
gmt pstext -R -J -D0/0.23 -O -K -F+f10p,Helvetica,gray10+jB -TO -Gwhite -W0.1 >> $out << END
$end_lon_C $end_lat_C C'
END

# ------------------------------------------------------------------------------------------------------------------- #
gmt psxy -R -J -T -O >> $out
gmt psconvert -Tf -A $out
gmt psconvert -Tg -A+r $out
evince ${out%.*}.pdf
# Delete eps file as it is large
rm ${out%.*}.eps
