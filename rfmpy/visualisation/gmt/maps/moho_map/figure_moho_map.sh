#######################################################################################################################
# Description: Map showing the seismic networks used
#
# Konstantinos Michailos
# Lausanne
# January 2022
#######################################################################################################################
# ------------------------------------------------------------------------------------------------------------------- #
# Output name
out=Moho_map.eps
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
gmt makecpt -C../files/imola.cpt -T20/80 -D+i -I > seis.cpt
gmt makecpt -C../files/bamako.cpt -T20/80 -D+i -I > seis.cpt
#gmt makecpt -C../files/batlow.cpt -T0/1200/200 -D+i > seis.cpt
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
gmt psscale -Dx7.0/8.9+o0/0i+w1.2i/0.08i+h+e -R -J -Cseis.cpt -Bxa20f10 -Bx+l"Moho depth (km)" -O -K  >> $out

#-Bxaf+l"topography" -By+lkm

echo Plot country names...
gmt pstext -R -J -O -K  -F+f6p,Helvetica,black+jBL+a0 -Gwhite >> $out << END
13 44 Adria
5 50.8 Europe
10 41.5 Liguria
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
rm -f adria.txt liguria.txt europe.txt


echo Plot initial 3D grid...
#awk '{print $1, $2}' files/initial_grid.txt |
#    gmt psxy -R -J -Sx.22 -W1.5p -Gred -O -K -t20 >> $out

echo Plot seismic stations...
#awk '{print $3, $2, $4}' files/rfs_calculated.txt | gmt psxy -i0,1,2 -Si.25 -R -J \
#-O -K -W.5p -Cseis.cpt -t5 >> $out
#awk '{print $3, $2, $1}' files/rfs_calculated.txt | gmt pstext -R -J -O -K -F+f2p,Helvetica,gray10 -Gwhite >> $out
#awk '{print $1, $2, $3}' ../files/moho_depths_1-6.txt | gmt psxy -i0,1,2 -Sd.18 -R -J \
#-O -K  -Cseis.cpt -t10 >> $out
awk '{print $1, $2, $3}' ../files/unc_moho_depths_Cross-section_1.txt | gmt psxy -i0,1,2 -Sd.15 -R -J -O -K -W.1p,black -Cseis.cpt -t10 >> $out
awk '{print $1, $2, $3}' ../files/unc_moho_depths_Cross-section_2.txt | gmt psxy -i0,1,2 -Sd.15 -R -J -O -K -W.1p,black -Cseis.cpt -t10 >> $out
awk '{print $1, $2, $3}' ../files/unc_moho_depths_Cross-section_3.txt | gmt psxy -i0,1,2 -Sd.15 -R -J -O -K -W.1p,black -Cseis.cpt -t10 >> $out
awk '{print $1, $2, $3}' ../files/unc_moho_depths_Cross-section_4.txt | gmt psxy -i0,1,2 -Sd.15 -R -J -O -K -W.1p,black -Cseis.cpt -t10 >> $out
awk '{print $1, $2, $3}' ../files/unc_moho_depths_Cross-section_5.txt | gmt psxy -i0,1,2 -Sd.15 -R -J -O -K -W.1p,black -Cseis.cpt -t10 >> $out
awk '{print $1, $2, $3}' ../files/unc_moho_depths_Cross-section_6.txt | gmt psxy -i0,1,2 -Sd.15 -R -J -O -K -W.1p,black -Cseis.cpt -t10 >> $out
awk '{print $1, $2, $3}' ../files/unc_moho_depths_Cross-section_7.txt | gmt psxy -i0,1,2 -Sd.15 -R -J -O -K -W.1p,black -Cseis.cpt -t10 >> $out
awk '{print $1, $2, $3}' ../files/unc_moho_depths_Cross-section_8.txt | gmt psxy -i0,1,2 -Sd.15 -R -J -O -K -W.1p,black -Cseis.cpt -t10 >> $out
awk '{print $1, $2, $3}' ../files/unc_moho_depths_Cross-section_9.txt | gmt psxy -i0,1,2 -Sd.15 -R -J -O -K -W.1p,black -Cseis.cpt -t10 >> $out
awk '{print $1, $2, $3}' ../files/unc_moho_depths_Cross-section_10.txt | gmt psxy -i0,1,2 -Sd.15 -R -J -O -K -W.1p,black -Cseis.cpt -t10 >> $out
awk '{print $1, $2, $3}' ../files/unc_moho_depths_Cross-section_11.txt | gmt psxy -i0,1,2 -Sd.15 -R -J -O -K -W.1p,black -Cseis.cpt -t10 >> $out
awk '{print $1, $2, $3}' ../files/unc_moho_depths_Cross-section_12.txt | gmt psxy -i0,1,2 -Sd.15 -R -J -O -K -W.1p,black -Cseis.cpt -t10 >> $out
awk '{print $1, $2, $3}' ../files/unc_moho_depths_Cross-section_13.txt | gmt psxy -i0,1,2 -Sd.15 -R -J -O -K -W.1p,black -Cseis.cpt -t10 >> $out
awk '{print $1, $2, $3}' ../files/unc_moho_depths_Cross-section_14.txt | gmt psxy -i0,1,2 -Sd.15 -R -J -O -K -W.1p,black -Cseis.cpt -t10 >> $out
awk '{print $1, $2, $3}' ../files/unc_moho_depths_Cross-section_15.txt | gmt psxy -i0,1,2 -Sd.15 -R -J -O -K -W.1p,black -Cseis.cpt -t10 >> $out
awk '{print $1, $2, $3}' ../files/unc_moho_depths_Cross-section_16.txt | gmt psxy -i0,1,2 -Sd.15 -R -J -O -K -W.1p,black -Cseis.cpt -t10 >> $out
awk '{print $1, $2, $3}' ../files/unc_moho_depths_Cross-section_17.txt | gmt psxy -i0,1,2 -Sd.15 -R -J -O -K -W.1p,black -Cseis.cpt -t10 >> $out
awk '{print $1, $2, $3}' ../files/unc_moho_depths_Cross-section_18.txt | gmt psxy -i0,1,2 -Sd.15 -R -J -O -K -W.1p,black -Cseis.cpt -t10 >> $out
awk '{print $1, $2, $3}' ../files/unc_moho_depths_Cross-section_19.txt | gmt psxy -i0,1,2 -Sd.15 -R -J -O -K -W.1p,black -Cseis.cpt -t10 >> $out
awk '{print $1, $2, $3}' ../files/unc_moho_depths_Cross-section_0.txt | gmt psxy -i0,1,2 -Sd.15 -R -J -O -K -W.1p,black -Cseis.cpt -t10 >> $out
awk '{print $1, $2, $3}' ../files/unc_moho_depths_Cross-section_-1.txt | gmt psxy -i0,1,2 -Sd.15 -R -J -O -K -W.1p,black -Cseis.cpt -t10 >> $out
awk '{print $1, $2, $3}' ../files/unc_moho_depths_Cross-section_20a.txt | gmt psxy -i0,1,2 -Sd.15 -R -J -O -K -W.1p,black -Cseis.cpt -t10 >> $out
awk '{print $1, $2, $3}' ../files/unc_moho_depths_Cross-section_20b.txt | gmt psxy -i0,1,2 -Sd.15 -R -J -O -K -W.1p,black -Cseis.cpt -t10 >> $out
awk '{print $1, $2, $3}' ../files/unc_moho_depths_Cross-section_20c.txt | gmt psxy -i0,1,2 -Sd.15 -R -J -O -K -W.1p,black -Cseis.cpt -t10 >> $out
awk '{print $1, $2, $3}' ../files/unc_moho_depths_Cross-section_21a.txt | gmt psxy -i0,1,2 -Sd.15 -R -J -O -K -W.1p,black -Cseis.cpt -t10 >> $out
awk '{print $1, $2, $3}' ../files/unc_moho_depths_Cross-section_21b.txt | gmt psxy -i0,1,2 -Sd.15 -R -J -O -K -W.1p,black -Cseis.cpt -t10 >> $out
awk '{print $1, $2, $3}' ../files/unc_moho_depths_Cross-section_21c.txt | gmt psxy -i0,1,2 -Sd.15 -R -J -O -K -W.1p,black -Cseis.cpt -t10 >> $out
awk '{print $1, $2, $3}' ../files/unc_moho_depths_Cross-section_22a.txt | gmt psxy -i0,1,2 -Sd.15 -R -J -O -K -W.1p,black -Cseis.cpt -t10 >> $out
awk '{print $1, $2, $3}' ../files/unc_moho_depths_Cross-section_22b.txt | gmt psxy -i0,1,2 -Sd.15 -R -J -O -K -W.1p,black -Cseis.cpt -t10 >> $out
#awk '{print $1, $2, $3}' ../files/unc_moho_depths_Cross-section_22c.txt | gmt psxy -i0,1,2 -Sd.15 -R -J -O -K -W.1p,black -Cseis.cpt -t10 >> $out
awk '{print $1, $2, $3}' ../files/unc_moho_depths_Cross-section_23a.txt | gmt psxy -i0,1,2 -Sd.15 -R -J -O -K -W.1p,black -Cseis.cpt -t10 >> $out
awk '{print $1, $2, $3}' ../files/unc_moho_depths_Cross-section_23b.txt | gmt psxy -i0,1,2 -Sd.15 -R -J -O -K -W.1p,black -Cseis.cpt -t10 >> $out
awk '{print $1, $2, $3}' ../files/unc_moho_depths_Cross-section_23c.txt | gmt psxy -i0,1,2 -Sd.15 -R -J -O -K -W.1p,black -Cseis.cpt -t10 >> $out
awk '{print $1, $2, $3}' ../files/unc_moho_depths_Cross-section_24a.txt | gmt psxy -i0,1,2 -Sd.15 -R -J -O -K -W.1p,black -Cseis.cpt -t10 >> $out
awk '{print $1, $2, $3}' ../files/unc_moho_depths_Cross-section_24b.txt | gmt psxy -i0,1,2 -Sd.15 -R -J -O -K -W.1p,black -Cseis.cpt -t10 >> $out
awk '{print $1, $2, $3}' ../files/unc_moho_depths_Cross-section_24c.txt | gmt psxy -i0,1,2 -Sd.15 -R -J -O -K -W.1p,black -Cseis.cpt -t10 >> $out
awk '{print $1, $2, $3}' ../files/unc_moho_depths_Cross-section_25a.txt | gmt psxy -i0,1,2 -Sd.15 -R -J -O -K -W.1p,black -Cseis.cpt -t10 >> $out
awk '{print $1, $2, $3}' ../files/unc_moho_depths_Cross-section_25b.txt | gmt psxy -i0,1,2 -Sd.15 -R -J -O -K -W.1p,black -Cseis.cpt -t10 >> $out
awk '{print $1, $2, $3}' ../files/unc_moho_depths_Cross-section_25c.txt | gmt psxy -i0,1,2 -Sd.15 -R -J -O -K -W.1p,black -Cseis.cpt -t10 >> $out
awk '{print $1, $2, $3}' ../files/unc_moho_depths_Cross-section_26a.txt | gmt psxy -i0,1,2 -Sd.15 -R -J -O -K -W.1p,black -Cseis.cpt -t10 >> $out
awk '{print $1, $2, $3}' ../files/unc_moho_depths_Cross-section_26b.txt | gmt psxy -i0,1,2 -Sd.15 -R -J -O -K -W.1p,black -Cseis.cpt -t10 >> $out
awk '{print $1, $2, $3}' ../files/unc_moho_depths_Cross-section_26c.txt | gmt psxy -i0,1,2 -Sd.15 -R -J -O -K -W.1p,black -Cseis.cpt -t10 >> $out
awk '{print $1, $2, $3}' ../files/unc_moho_depths_Cross-section_27a.txt | gmt psxy -i0,1,2 -Sd.15 -R -J -O -K -W.1p,black -Cseis.cpt -t10 >> $out
awk '{print $1, $2, $3}' ../files/unc_moho_depths_Cross-section_27b.txt | gmt psxy -i0,1,2 -Sd.15 -R -J -O -K -W.1p,black -Cseis.cpt -t10 >> $out
awk '{print $1, $2, $3}' ../files/unc_moho_depths_Cross-section_27c.txt | gmt psxy -i0,1,2 -Sd.15 -R -J -O -K -W.1p,black -Cseis.cpt -t10 >> $out
awk '{print $1, $2, $3}' ../files/unc_moho_depths_Cross-section_28a.txt | gmt psxy -i0,1,2 -Sd.15 -R -J -O -K -W.1p,black -Cseis.cpt -t10 >> $out
awk '{print $1, $2, $3}' ../files/unc_moho_depths_Cross-section_28b.txt | gmt psxy -i0,1,2 -Sd.15 -R -J -O -K -W.1p,black -Cseis.cpt -t10 >> $out
awk '{print $1, $2, $3}' ../files/unc_moho_depths_Cross-section_28c.txt | gmt psxy -i0,1,2 -Sd.15 -R -J -O -K -W.1p,black -Cseis.cpt -t10 >> $out
awk '{print $1, $2, $3}' ../files/unc_moho_depths_Cross-section_29a.txt | gmt psxy -i0,1,2 -Sd.15 -R -J -O -K -W.1p,black -Cseis.cpt -t10 >> $out
awk '{print $1, $2, $3}' ../files/unc_moho_depths_Cross-section_29b.txt | gmt psxy -i0,1,2 -Sd.15 -R -J -O -K -W.1p,black -Cseis.cpt -t10 >> $out
awk '{print $1, $2, $3}' ../files/unc_moho_depths_Cross-section_29c.txt | gmt psxy -i0,1,2 -Sd.15 -R -J -O -K -W.1p,black -Cseis.cpt -t10 >> $out
awk '{print $1, $2, $3}' ../files/unc_moho_depths_Cross-section_30a.txt | gmt psxy -i0,1,2 -Sd.15 -R -J -O -K -W.1p,black -Cseis.cpt -t10 >> $out
awk '{print $1, $2, $3}' ../files/unc_moho_depths_Cross-section_30b.txt | gmt psxy -i0,1,2 -Sd.15 -R -J -O -K -W.1p,black -Cseis.cpt -t10 >> $out
awk '{print $1, $2, $3}' ../files/unc_moho_depths_Cross-section_30c.txt | gmt psxy -i0,1,2 -Sd.15 -R -J -O -K -W.1p,black -Cseis.cpt -t10 >> $out




awk '{print $1, $2, $3}' ../files/moho_depths_Cross-section_1.txt | gmt psxy -i0,1,2 -Sd.18 -R -J -O -K  -Cseis.cpt -t5 >> $out
awk '{print $1, $2, $3}' ../files/moho_depths_Cross-section_2.txt | gmt psxy -i0,1,2 -Sd.18 -R -J -O -K  -Cseis.cpt -t5 >> $out
awk '{print $1, $2, $3}' ../files/moho_depths_Cross-section_3.txt | gmt psxy -i0,1,2 -Sd.18 -R -J -O -K  -Cseis.cpt -t5 >> $out
awk '{print $1, $2, $3}' ../files/moho_depths_Cross-section_4.txt | gmt psxy -i0,1,2 -Sd.18 -R -J -O -K  -Cseis.cpt -t5 >> $out
awk '{print $1, $2, $3}' ../files/moho_depths_Cross-section_5.txt | gmt psxy -i0,1,2 -Sd.18 -R -J -O -K  -Cseis.cpt -t5 >> $out
awk '{print $1, $2, $3}' ../files/moho_depths_Cross-section_6.txt | gmt psxy -i0,1,2 -Sd.18 -R -J -O -K  -Cseis.cpt -t5 >> $out
awk '{print $1, $2, $3}' ../files/moho_depths_Cross-section_7.txt | gmt psxy -i0,1,2 -Sd.18 -R -J -O -K  -Cseis.cpt -t5 >> $out
awk '{print $1, $2, $3}' ../files/moho_depths_Cross-section_8.txt | gmt psxy -i0,1,2 -Sd.18 -R -J -O -K  -Cseis.cpt -t5 >> $out
awk '{print $1, $2, $3}' ../files/moho_depths_Cross-section_9.txt | gmt psxy -i0,1,2 -Sd.18 -R -J -O -K  -Cseis.cpt -t5 >> $out
awk '{print $1, $2, $3}' ../files/moho_depths_Cross-section_10.txt | gmt psxy -i0,1,2 -Sd.18 -R -J -O -K  -Cseis.cpt -t5 >> $out
awk '{print $1, $2, $3}' ../files/moho_depths_Cross-section_11.txt | gmt psxy -i0,1,2 -Sd.18 -R -J -O -K  -Cseis.cpt -t5 >> $out
awk '{print $1, $2, $3}' ../files/moho_depths_Cross-section_12.txt | gmt psxy -i0,1,2 -Sd.18 -R -J -O -K  -Cseis.cpt -t5 >> $out
awk '{print $1, $2, $3}' ../files/moho_depths_Cross-section_13.txt | gmt psxy -i0,1,2 -Sd.18 -R -J -O -K  -Cseis.cpt -t5 >> $out
awk '{print $1, $2, $3}' ../files/moho_depths_Cross-section_14.txt | gmt psxy -i0,1,2 -Sd.18 -R -J -O -K  -Cseis.cpt -t5 >> $out
awk '{print $1, $2, $3}' ../files/moho_depths_Cross-section_15.txt | gmt psxy -i0,1,2 -Sd.18 -R -J -O -K  -Cseis.cpt -t5 >> $out
awk '{print $1, $2, $3}' ../files/moho_depths_Cross-section_16.txt | gmt psxy -i0,1,2 -Sd.18 -R -J -O -K  -Cseis.cpt -t5 >> $out
awk '{print $1, $2, $3}' ../files/moho_depths_Cross-section_17.txt | gmt psxy -i0,1,2 -Sd.18 -R -J -O -K  -Cseis.cpt -t5 >> $out
awk '{print $1, $2, $3}' ../files/moho_depths_Cross-section_18.txt | gmt psxy -i0,1,2 -Sd.18 -R -J -O -K  -Cseis.cpt -t5 >> $out
awk '{print $1, $2, $3}' ../files/moho_depths_Cross-section_19.txt | gmt psxy -i0,1,2 -Sd.18 -R -J -O -K  -Cseis.cpt -t5 >> $out
awk '{print $1, $2, $3}' ../files/moho_depths_Cross-section_0.txt | gmt psxy -i0,1,2 -Sd.18 -R -J -O -K  -Cseis.cpt -t5 >> $out
awk '{print $1, $2, $3}' ../files/moho_depths_Cross-section_-1.txt | gmt psxy -i0,1,2 -Sd.18 -R -J -O -K  -Cseis.cpt -t5 >> $out
awk '{print $1, $2, $3}' ../files/moho_depths_Cross-section_20a.txt | gmt psxy -i0,1,2 -Sd.18 -R -J -O -K  -Cseis.cpt -t5 >> $out
#awk '{print $1, $2, $3}' ../files/moho_depths_Cross-section_20b.txt | gmt psxy -i0,1,2 -Sd.18 -R -J -O -K  -Cseis.cpt -t10 >> $out
#awk '{print $1, $2, $3}' ../files/moho_depths_Cross-section_20c.txt | gmt psxy -i0,1,2 -Sd.18 -R -J -O -K  -Cseis.cpt -t10 >> $out
awk '{print $1, $2, $3}' ../files/moho_depths_Cross-section_21a.txt | gmt psxy -i0,1,2 -Sd.18 -R -J -O -K  -Cseis.cpt -t5 >> $out
awk '{print $1, $2, $3}' ../files/moho_depths_Cross-section_21b.txt | gmt psxy -i0,1,2 -Sd.18 -R -J -O -K  -Cseis.cpt -t5 >> $out
#awk '{print $1, $2, $3}' ../files/moho_depths_Cross-section_21c.txt | gmt psxy -i0,1,2 -Sd.18 -R -J -O -K  -Cseis.cpt -t10 >> $out
awk '{print $1, $2, $3}' ../files/moho_depths_Cross-section_22a.txt | gmt psxy -i0,1,2 -Sd.18 -R -J -O -K  -Cseis.cpt -t5 >> $out
awk '{print $1, $2, $3}' ../files/moho_depths_Cross-section_22b.txt | gmt psxy -i0,1,2 -Sd.18 -R -J -O -K  -Cseis.cpt -t5 >> $out
#awk '{print $1, $2, $3}' ../files/moho_depths_Cross-section_22c.txt | gmt psxy -i0,1,2 -Sd.18 -R -J -O -K  -Cseis.cpt -t10 >> $out
awk '{print $1, $2, $3}' ../files/moho_depths_Cross-section_23a.txt | gmt psxy -i0,1,2 -Sd.18 -R -J -O -K  -Cseis.cpt -t5 >> $out
awk '{print $1, $2, $3}' ../files/moho_depths_Cross-section_23b.txt | gmt psxy -i0,1,2 -Sd.18 -R -J -O -K  -Cseis.cpt -t5 >> $out
#awk '{print $1, $2, $3}' ../files/moho_depths_Cross-section_23c.txt | gmt psxy -i0,1,2 -Sd.18 -R -J -O -K  -Cseis.cpt -t10 >> $out
awk '{print $1, $2, $3}' ../files/moho_depths_Cross-section_24a.txt | gmt psxy -i0,1,2 -Sd.18 -R -J -O -K  -Cseis.cpt -t5 >> $out
awk '{print $1, $2, $3}' ../files/moho_depths_Cross-section_24b.txt | gmt psxy -i0,1,2 -Sd.18 -R -J -O -K  -Cseis.cpt -t5 >> $out
awk '{print $1, $2, $3}' ../files/moho_depths_Cross-section_24c.txt | gmt psxy -i0,1,2 -Sd.18 -R -J -O -K  -Cseis.cpt -t5 >> $out
awk '{print $1, $2, $3}' ../files/moho_depths_Cross-section_25a.txt | gmt psxy -i0,1,2 -Sd.18 -R -J -O -K  -Cseis.cpt -t5 >> $out
awk '{print $1, $2, $3}' ../files/moho_depths_Cross-section_25b.txt | gmt psxy -i0,1,2 -Sd.18 -R -J -O -K  -Cseis.cpt -t5 >> $out
awk '{print $1, $2, $3}' ../files/moho_depths_Cross-section_25c.txt | gmt psxy -i0,1,2 -Sd.18 -R -J -O -K  -Cseis.cpt -t5 >> $out
awk '{print $1, $2, $3}' ../files/moho_depths_Cross-section_26a.txt | gmt psxy -i0,1,2 -Sd.18 -R -J -O -K  -Cseis.cpt -t5 >> $out
awk '{print $1, $2, $3}' ../files/moho_depths_Cross-section_26b.txt | gmt psxy -i0,1,2 -Sd.18 -R -J -O -K  -Cseis.cpt -t5 >> $out
awk '{print $1, $2, $3}' ../files/moho_depths_Cross-section_26c.txt | gmt psxy -i0,1,2 -Sd.18 -R -J -O -K  -Cseis.cpt -t5 >> $out
awk '{print $1, $2, $3}' ../files/moho_depths_Cross-section_27a.txt | gmt psxy -i0,1,2 -Sd.18 -R -J -O -K  -Cseis.cpt -t5 >> $out
awk '{print $1, $2, $3}' ../files/moho_depths_Cross-section_27b.txt | gmt psxy -i0,1,2 -Sd.18 -R -J -O -K  -Cseis.cpt -t5 >> $out
#awk '{print $1, $2, $3}' ../files/moho_depths_Cross-section_27c.txt | gmt psxy -i0,1,2 -Sd.18 -R -J -O -K  -Cseis.cpt -t10 >> $out
awk '{print $1, $2, $3}' ../files/moho_depths_Cross-section_28a.txt | gmt psxy -i0,1,2 -Sd.18 -R -J -O -K  -Cseis.cpt -t5 >> $out
awk '{print $1, $2, $3}' ../files/moho_depths_Cross-section_28b.txt | gmt psxy -i0,1,2 -Sd.18 -R -J -O -K  -Cseis.cpt -t5 >> $out
awk '{print $1, $2, $3}' ../files/moho_depths_Cross-section_28c.txt | gmt psxy -i0,1,2 -Sd.18 -R -J -O -K  -Cseis.cpt -t5 >> $out
#awk '{print $1, $2, $3}' ../files/moho_depths_Cross-section_29a.txt | gmt psxy -i0,1,2 -Sd.18 -R -J -O -K  -Cseis.cpt -t10 >> $out
awk '{print $1, $2, $3}' ../files/moho_depths_Cross-section_29b.txt | gmt psxy -i0,1,2 -Sd.18 -R -J -O -K  -Cseis.cpt -t5 >> $out
awk '{print $1, $2, $3}' ../files/moho_depths_Cross-section_29c.txt | gmt psxy -i0,1,2 -Sd.18 -R -J -O -K  -Cseis.cpt -t5 >> $out
#awk '{print $1, $2, $3}' ../files/moho_depths_Cross-section_30a.txt | gmt psxy -i0,1,2 -Sd.18 -R -J -O -K  -Cseis.cpt -t10 >> $out
awk '{print $1, $2, $3}' ../files/moho_depths_Cross-section_30b.txt | gmt psxy -i0,1,2 -Sd.18 -R -J -O -K  -Cseis.cpt -t5 >> $out
awk '{print $1, $2, $3}' ../files/moho_depths_Cross-section_30c.txt | gmt psxy -i0,1,2 -Sd.18 -R -J -O -K  -Cseis.cpt -t5 >> $out


#awk '{print $1, $2, $3}' ../files/unc_moho_depths_Cross-section_7.txt | gmt psxy -i0,1,2 -Sd.18 -R -J -O -K -W.5p,black -Cseis.cpt -t10 >> $out
#awk '{print $1, $2, $3}' ../files/moho_depths_Cross-section_10.txt | gmt psxy -i0,1,2 -Sd.18 -R -J -O -K  -Cseis.cpt -t10 >> $out
#awk '{print $1, $2, $3}' ../files/unc_moho_depths_Cross-section_10.txt | gmt psxy -i0,1,2 -Sd.18 -R -J -O -K -W.5p,black -Cseis.cpt -t10 >> $out


#awk '{print $3, $2, $4}' ../files/pacase.txt | gmt psxy -i0,1,2 -Ss.2 -R -J \
#-O -K -W.5p -Cseis.cpt -t10 >> $out
# -=================================================================================================================- #
# ------------------------------------------------------------------------------------------------------------------- #
echo Create legend...
gmt set FONT_ANNOT_PRIMARY 7
gmt pslegend <<END -R -J -Dx3.5i/0.3i+w0i/0.0i/TC -C0.07i/0.1i -F+gwhite+pthin -P -O -K --FONT_ANNOT_PRIMARY=5p >> $out
G -0.05i
H 7 Manual picks
D0.1i 0.5p
G .04i
S .05i d .08i darkolivegreen 0.7p,darkolivegreen 0.18i Moho depth
G .06i
S .05i d .06i darkolivegreen 0.7p,black 0.18i Uncertain Moho depth
G .03i
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



gmt pstext -R -J -O -K  -F+f5,Helvetica,black+jBL+a0 -Gwhite >> $out << END
5 42 1
6 42 2
7 42 3
8 42 4
9 42 5
10 42 6
11 42 7
12 42 8
13 42 9
14 42 10
15 42 11
16 42 12
17 42 13
18 45 14
19 45 15
20 45 16
21 45 17
22 45 18
22 45 19
4 42 0
3 42 -1
1.5 43 20
1.5 43.7 21
1.5 44.4 22
1.5 45.1 23
1.5 45.8 24
1.5 46.5 25
1.5 47.2 26
1.5 47.9 27
1.5 48.6 28
1.5 49.3 29
1.5 50.0 30
END

# ------------------------------------------------------------------------------------------------------------------- #
gmt psxy -R -J -T -O >> $out
gmt psconvert -Tf -A $out
evince ${out%.*}.pdf
