#######################################################################################################################
# Description: Map showing the seismic networks used
#
# Konstantinos Michailos
# Lausanne
# January 2022
#######################################################################################################################
# ------------------------------------------------------------------------------------------------------------------- #
# Output name
out=map_1.eps
# ------------------------------------------------------------------------------------------------------------------- #
# Define stuffz (area plotted, size of letters, etc)
gmt set FORMAT_GEO_MAP D
gmt set FORMAT_GEO_MAP D
gmt set FONT_ANNOT_PRIMARY Helvetica
gmt set FONT_ANNOT_PRIMARY 8
gmt set FONT_LABEL Helvetica
gmt set LABEL_FONT_SIZE 7
gmt set MAP_FRAME_TYPE plain
# Directory containing topo grd file
topodir="/home/kmichall/Desktop/topo"
# Map boundaries
north=52
south=41
east=20
west=0
proj='-JM6i'
# ------------------------------------------------------------------------------------------------------------------- #
echo Make basemap...
gmt pscoast -W1/0.05 -Dl $proj -R$west/$east/$south/$north -K -Y1 -B5WSen -P > $out
# ------------------------------------------------------------------------------------------------------------------- #
echo Plot topo....
#gmt grdimage -R -J /home/kmichall/Desktop/topo/topo.0.20.40.55.3sec.grd -CFrance2.cpt -O -K >> $out
gmt grdimage -R -J /home/kmichall/Desktop/topo/ETOPO1_Bed_g_gmt4.grd -Cmy_topo.cpt -O -K >> $out
# ------------------------------------------------------------------------------------------------------------------- #
gmt pscoast -W1/0.05 -Df -J -R -K -O -P -Sazure1 -N1/0.05p,black -L2.2/47.8/48/200+l+u >> $out


echo Plotting faults and stuff...
# 250km distance line from a smoothed 800m elevation contour
gmt psxy -R -J d250km.dat -W1.5p,gray20 -O -K >> $out

# ---------
echo Create cpt...
gmt makecpt -Cviridis -T40/110/10  > seis.cpt
gmt makecpt -Chot -T0/250/50 -D+i -I > seis.cpt

echo Plot scale...
gmt psscale -Dx0.25/10+o0/0i+w1.2i/0.08i+h+e -R -J -Cmy_topo.cpt -Bx500f250 -Bx+l"Topography (m)" \
 -O -K --FONT_ANNOT_PRIMARY=8p >> $out
#gmt psscale -Dx0.2/8.5+o0/0i+w1.2i/0.08i+h+e -R -J  -Cseis.cpt -Bx50f50 -Bx+l"Number of RFs" \
#-O -K --FONT_ANNOT_PRIMARY=8p >> $out

echo Plot country names...
#gmt pstext -R -J -O -K  -F+f6p,Helvetica,black+jBL+a0 -Gwhite >> $out << END
#15.2 45.8 SLOVENIA
#15.7 45.5 CROATIA
#16.5 45.1 BOSNIA AND HERZEGOVINA
#END
# -=================================================================================================================- #
# ------------------------------------------------------------------------------------------------------------------- #
# -=================================================================================================================- #

echo Plot initial 3D grid...
awk '{print $1, $2}' files/initial_grid.txt |
    gmt psxy -R -J -Sx.22 -W1.5p -Gred -O -K -t20 >> $out

echo Plot seismic stations...
#awk '{print $3, $2}' rfs_calculated.txt |
#    gmt psxy -R -J -Si.22 -W0.5p -Gred -O -K -t20 >> $out
#awk '{print $3, $2, $4}' files/rfs_calculated.txt | gmt psxy -i0,1,2 -Si.25 -R -J \
#-O -K -W.5p -Cseis.cpt -t5 >> $out
#awk '{print $3, $2, $1}' files/rfs_calculated.txt | gmt pstext -R -J -O -K -F+f2p,Helvetica,gray10 -Gwhite >> $out

#awk '{print $3, $2, $4}' number_of_waveforms.txt | gmt psxy -i0,1,2 -Si.25 -R -J \
#-O -K -W.5p -Cseis.cpt -t10 >> $out

# New data
#awk '{print $3, $2}' new_data.txt |
#    gmt psxy -R -J -St.22 -W0.5p,gray -Gdodgerblue -O -K -t0 >> $out
#awk '{print $3, $2}' ZJ.txt |
#    gmt psxy -R -J -St.25 -W0.5p,black -Gdodgerblue -t0 -O -K  >> $out
#awk '{print $3, $2}' missing_FR.txt |
#    gmt psxy -R -J -St.25 -W0.5p -Gorange -t0 -O -K  >> $out
#awk '{print $3, $2}' missing_FR_2.txt |
#    gmt psxy -R -J -St.25 -W0.5p -Gyellow -t0 -O -K  >> $out
#awk '{print $3, $2}' stations_3.txt |
#    gmt psxy -R -J -Si.25 -W0.5p -Gred -t0 -O -K  >> $out
#awk '{print $3, $2}' cifalps.txt |
#    gmt psxy -R -J -St.25 -W0.5p -Gpurple -t0 -O -K  >> $out
#awk '{print $3, $2}' station_5.txt |
#    gmt psxy -R -J -St.25 -W0.5p -Gyellow -t0 -O -K  >> $out
# -=================================================================================================================- #
# ------------------------------------------------------------------------------------------------------------------- #
echo Create legend...
#gmt set FONT_ANNOT_PRIMARY 8
#gmt pslegend <<END -R -J -Dx0.05i/4.15i+w0i/0.0i/TC -C0.07i/0.1i -F+gwhite+pthin -P -O -K >> $out
#G -0.05i
#H 9 Seismic networks
#D0.1i 0.5p
#G .04i
#S .04i i .11i red 0.5p 0.18i AASN, EASI, CIFALPS
#G .05i
##S .04i t .11i orange 0.2p 0.18i ZJ (Data to be processed)
##G .07i
##S .04i t .11i orange 0.2p,red 0.18i PL, HU, CZ, SK (Data to be processed)
##G 0.05i
##S .04i - .14i gray0 thick 0.18i Active faults
#END

echo Plot cross section lines
start_lon='8'
start_lat='45.5'
end_lon='15'
end_lat='50'

gmt psxy << END -R -J -O -W3,dodgerblue -K>> $out
$start_lon $start_lat
$end_lon $end_lat
END
gmt pstext -R -J -D0/0.23 -O -K -F+f12p,Helvetica,gray10 -TO -Gwhite -W0.1 >> $out << END
$start_lon $start_lat A
END
gmt pstext -R -J -D0/0.23 -O -K -F+f12p,Helvetica,gray10 -TO -Gwhite -W0.1 >> $out << END
$end_lon $end_lat B
END




# ------------------------------------------------------------------------------------------------------------------- #
gmt psxy -R -J -T -O >> $out
gmt psconvert -Tf -A $out
evince ${out%.*}.pdf

#awk '{print $1, $2, $3}' xyz.txt| gmt xyz2grd -R0/100/0/100 -I1/5 -Gt.nc
#gmt grdimage t.nc -Baf -B+t -R0/100/0/100 -JP10c+z -png map