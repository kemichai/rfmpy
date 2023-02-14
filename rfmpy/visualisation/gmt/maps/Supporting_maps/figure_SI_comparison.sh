#######################################################################################################################
# Description: Map showing the seismic networks used
#
# Konstantinos Michailos
# Lausanne
# January 2022
#######################################################################################################################
# ------------------------------------------------------------------------------------------------------------------- #
# Output name
out=compare2previous_studies.eps
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

echo Create cpt...
#gmt makecpt -C../files/imola.cpt -T20/80 -D+i -I > seis.cpt
#gmt makecpt -C../files/bamako.cpt -T20/80/5 -D+i -I > seis.cpt
gmt makecpt -Cvik.cpt -T-25/25 -D+i -I > seis.cpt
#gmt makecpt -C../files/moho_picks/batlow.cpt -T0/1200/200 -D+i > seis.cpt
#gmt makecpt -Chot -T0/1200/200 -D+i > seis.cpt
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



echo Plot scale...
gmt psscale -Dx2.5/8.9+o0/0i+w1.2i/0.08i+h+e -R -J -Cmy_topo.cpt -Bx1000f500 -Bx+l"Topography (m)" -O -K  >> $out
gmt psscale -Dx7.0/8.9+o0/0i+w1.2i/0.08i+h+e -R -J -Cseis.cpt -Bxa10f5 -Bx+l"Moho depth residuals (km)" -O -K  >> $out

#-Bxaf+l"topography" -By+lkm

echo Plot country names...
gmt pstext -R -J -O -K  -F+f6p,Helvetica,black+jBL+a0 -Gwhite >> $out << END
13 44 Adria
5 50.8 Europe
10 41.5 Liguria
END
gmt pstext -R -J -O -K  -F+f13p,Helvetica,black+jBL+a0 -Gwhite -W.5p,black >> $out << END
1. 51 a)
END

# -=================================================================================================================- #
# ------------------------------------------------------------------------------------------------------------------- #
# -=================================================================================================================- #

cat <<- END > adria.txt
17.0169     47.1074
15.4473     47.1069
14.9241     46.9806
14.0970     46.8961
13.7595     46.7905
11.2616     46.7694
10.0464     46.3680
8.6624     46.3891
7.4135     45.3116
7.3629     44.4665
7.8186     44.1708
8.2405     44.4454
9.3713     44.3187
10.7722     44.0651
13.4557     42.2694
END

cat <<- END > liguria.txt
8.0106     41.3210
7.7082     41.6486
7.7182     41.8845
7.4158     42.2776
7.4057     42.6708
8.0005     43.6143
8.0307     43.7322
8.0408     44.3874
8.2223     44.3874
8.2324     44.4791
9.4320     44.4922
9.4924     44.3874
9.7545     44.4005
9.7848     44.3088
10.2686    44.3088
10.2787    44.2039
10.9037    44.2039
10.8836    44.0991
11.4464    44.1122
12.5469    43.1949
12.5973    43.2342
13.6154    42.1990
13.6961    42.2252
13.9884    41.9107
END

cat <<- END > europe.txt
8.0363    44.1327
7.6331    44.4079
7.4234    44.4079
7.4073    45.1024
7.9073    45.6265
8.1008    45.6265
8.3105    45.9148
8.4718    45.9148
8.7298    46.2162
9.8105    46.2162
9.8427    46.3210
10.6653    46.3079
11.1815    46.6224
12.7782    46.6224
12.8427    46.7273
14.6976    46.7142
14.7782    46.8190
15.2782    46.8059
15.4556    47.0156
17.0524    47.0156
END

gmt psxy liguria.txt -Wthick,black -O -K -R -J >> $out
gmt psxy europe.txt -Wthick,black -O -K -R -J >> $out
gmt psxy adria.txt -Wthick,black -O -K -R -J >> $out


echo Plot initial 3D grid...
#awk '{print $1, $2}' files/initial_grid.txt |
#    gmt psxy -R -J -Sx.22 -W1.5p -Gred -O -K -t20 >> $out

echo Plot seismic stations...
#awk '{print $3, $2, $4}' files/rfs_calculated.txt | gmt psxy -i0,1,2 -Si.25 -R -J \
#-O -K -W.5p -Cseis.cpt -t5 >> $out
#awk '{print $3, $2, $1}' files/rfs_calculated.txt | gmt pstext -R -J -O -K -F+f2p,Helvetica,gray10 -Gwhite >> $out
#awk '{print $1, $2, $3}' ../files/moho_picks/moho_depths_1-6.txt | gmt psxy -i0,1,2 -Sd.15 -R -J \
#-O -K  -Cseis.cpt -t10 >> $out

awk '{print $1, $2, $3}' diff2spada.txt | gmt psxy -i0,1,2 -Sc.20 -R -J -O -K  -Cseis.cpt -t5 >> $out
awk '{print $1, $2}' ../files/moho_picks/moho_depths_all.dat | gmt psxy -R -J -Sc.03 -Gblack -O -K -t20 >> $out

# -=================================================================================================================- #
# ------------------------------------------------------------------------------------------------------------------- #
echo Create legend...
gmt pslegend <<END -R -J -Dx3.3i/0.3i+w0i/0.0i/TC -C0.1i/0.1i -F+gwhite+pthin -P -O -K --FONT_ANNOT_PRIMARY=6p >> $out
G -0.05i
H 6 Comparison with Spada et al. 2013
D0.1i 0.5p
G .04i
S .05i c .07i dodgerblue2 0.7p,dodgerblue2 0.18i Moho - MohoSP
G .06i
S .05i c .02i black 0.2p,black 0.18i Moho depth pick locations
G .03i
END



gmt pscoast -W1/0.05 -Df $proj -R$west/$east/$south/$north -K -B5WSen -O -P -Y-9 >> $out
gmt grdimage -R -J $topodir/ETOPO1_Bed_g_gmt4.grd -Cmy_topo.cpt -O -K >> $out
gmt pscoast -W1/0.05 -Df -J -R -K -O -P -Sazure1 -N1/0.05p,black -L3.4/49.7/48/200+l+u >> $out


awk '{print $1, $2, $3}' diff2Grad.txt | gmt psxy -i0,1,2 -Sc.20 -R -J -O -K  -Cseis.cpt -t5 >> $out
awk '{print $1, $2}' ../files/moho_picks/moho_depths_all.dat | gmt psxy -R -J -Sc.03 -Gblack -O -K -t20 >> $out




gmt psxy liguria.txt -Wthick,black -O -K -R -J >> $out
gmt psxy europe.txt -Wthick,black -O -K -R -J >> $out
gmt psxy adria.txt -Wthick,black -O -K -R -J >> $out
rm -f adria.txt liguria.txt europe.txt

echo Plot country names...
gmt pstext -R -J -O -K  -F+f6p,Helvetica,black+jBL+a0 -Gwhite >> $out << END
13 44 Adria
5 50.8 Europe
10 41.5 Liguria
END
gmt pstext -R -J -O -K  -F+f13p,Helvetica,black+jBL+a0 -Gwhite -W.5p,black >> $out << END
1. 51 b)
END

echo Create legend...
gmt pslegend <<END -R -J -Dx3.3i/0.3i+w0i/0.0i/TC -C0.1i/0.1i -F+gwhite+pthin -P -O -K --FONT_ANNOT_PRIMARY=6p >> $out
G -0.05i
H 6 Comparison with Grad and Tira 2009
D0.1i 0.5p
G .04i
S .05i c .07i dodgerblue2 0.7p,dodgerblue2 0.18i Moho - MohoGD
G .06i
S .05i c .02i black 0.2p,black 0.18i Moho depth pick locations
G .03i
END

# ------------------------------------------------------------------------------------------------------------------- #
gmt psxy -R -J -T -O >> $out
gmt psconvert -Tf -A $out
evince ${out%.*}.pdf
