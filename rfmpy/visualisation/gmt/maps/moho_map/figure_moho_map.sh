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
picks_dir="/home/kmichailos/Desktop/codes/github/rfmpy/rfmpy/visualisation/gmt/cross_sections/appendix_moho_picks/pick_files"

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
gmt makecpt -Cviridis -T20/60 -D+i -I > seis.cpt
#gmt makecpt -C$picks_dir/batlow.cpt -T0/1200/200 -D+i > seis.cpt
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
gmt psscale -Dx7.0/8.9+o0/0i+w1.2i/0.08i+h+e -R -J -Cseis.cpt -Bxa10f5 -Bx+l"Moho depth (km)" -O -K  >> $out

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
#awk '{print $1, $2, $3}' $picks_dir/moho_depths_1-6.txt | gmt psxy -i0,1,2 -Sd.15 -R -J \
#-O -K  -Cseis.cpt -t10 >> $out
awk '{print $1, $2, $3}' $picks_dir/unc_moho_depths_Cross-section_1.txt | gmt psxy -i0,1,2 -Sd.15 -R -J -O -K -W.1p,black -Cseis.cpt -t10 >> $out
awk '{print $1, $2, $3}' $picks_dir/unc_moho_depths_Cross-section_2.txt | gmt psxy -i0,1,2 -Sd.15 -R -J -O -K -W.1p,black -Cseis.cpt -t10 >> $out
awk '{print $1, $2, $3}' $picks_dir/unc_moho_depths_Cross-section_3.txt | gmt psxy -i0,1,2 -Sd.15 -R -J -O -K -W.1p,black -Cseis.cpt -t10 >> $out
awk '{print $1, $2, $3}' $picks_dir/unc_moho_depths_Cross-section_4.txt | gmt psxy -i0,1,2 -Sd.15 -R -J -O -K -W.1p,black -Cseis.cpt -t10 >> $out
awk '{print $1, $2, $3}' $picks_dir/unc_moho_depths_Cross-section_5.txt | gmt psxy -i0,1,2 -Sd.15 -R -J -O -K -W.1p,black -Cseis.cpt -t10 >> $out
awk '{print $1, $2, $3}' $picks_dir/unc_moho_depths_Cross-section_6.txt | gmt psxy -i0,1,2 -Sd.15 -R -J -O -K -W.1p,black -Cseis.cpt -t10 >> $out
awk '{print $1, $2, $3}' $picks_dir/unc_moho_depths_Cross-section_7.txt | gmt psxy -i0,1,2 -Sd.15 -R -J -O -K -W.1p,black -Cseis.cpt -t10 >> $out
awk '{print $1, $2, $3}' $picks_dir/unc_moho_depths_Cross-section_8.txt | gmt psxy -i0,1,2 -Sd.15 -R -J -O -K -W.1p,black -Cseis.cpt -t10 >> $out
awk '{print $1, $2, $3}' $picks_dir/unc_moho_depths_Cross-section_9.txt | gmt psxy -i0,1,2 -Sd.15 -R -J -O -K -W.1p,black -Cseis.cpt -t10 >> $out
awk '{print $1, $2, $3}' $picks_dir/unc_moho_depths_Cross-section_10.txt | gmt psxy -i0,1,2 -Sd.15 -R -J -O -K -W.1p,black -Cseis.cpt -t10 >> $out
awk '{print $1, $2, $3}' $picks_dir/unc_moho_depths_Cross-section_11.txt | gmt psxy -i0,1,2 -Sd.15 -R -J -O -K -W.1p,black -Cseis.cpt -t10 >> $out
awk '{print $1, $2, $3}' $picks_dir/unc_moho_depths_Cross-section_12.txt | gmt psxy -i0,1,2 -Sd.15 -R -J -O -K -W.1p,black -Cseis.cpt -t10 >> $out
awk '{print $1, $2, $3}' $picks_dir/unc_moho_depths_Cross-section_13.txt | gmt psxy -i0,1,2 -Sd.15 -R -J -O -K -W.1p,black -Cseis.cpt -t10 >> $out
awk '{print $1, $2, $3}' $picks_dir/unc_moho_depths_Cross-section_14.txt | gmt psxy -i0,1,2 -Sd.15 -R -J -O -K -W.1p,black -Cseis.cpt -t10 >> $out
awk '{print $1, $2, $3}' $picks_dir/unc_moho_depths_Cross-section_15.txt | gmt psxy -i0,1,2 -Sd.15 -R -J -O -K -W.1p,black -Cseis.cpt -t10 >> $out
awk '{print $1, $2, $3}' $picks_dir/unc_moho_depths_Cross-section_16.txt | gmt psxy -i0,1,2 -Sd.15 -R -J -O -K -W.1p,black -Cseis.cpt -t10 >> $out
awk '{print $1, $2, $3}' $picks_dir/unc_moho_depths_Cross-section_17.txt | gmt psxy -i0,1,2 -Sd.15 -R -J -O -K -W.1p,black -Cseis.cpt -t10 >> $out
awk '{print $1, $2, $3}' $picks_dir/unc_moho_depths_Cross-section_18.txt | gmt psxy -i0,1,2 -Sd.15 -R -J -O -K -W.1p,black -Cseis.cpt -t10 >> $out
awk '{print $1, $2, $3}' $picks_dir/unc_moho_depths_Cross-section_19.txt | gmt psxy -i0,1,2 -Sd.15 -R -J -O -K -W.1p,black -Cseis.cpt -t10 >> $out
awk '{print $1, $2, $3}' $picks_dir/unc_moho_depths_Cross-section_0.txt | gmt psxy -i0,1,2 -Sd.15 -R -J -O -K -W.1p,black -Cseis.cpt -t10 >> $out
awk '{print $1, $2, $3}' $picks_dir/unc_moho_depths_Cross-section_20.txt | gmt psxy -i0,1,2 -Sd.15 -R -J -O -K -W.1p,black -Cseis.cpt -t10 >> $out
awk '{print $1, $2, $3}' $picks_dir/unc_moho_depths_Cross-section_21.txt | gmt psxy -i0,1,2 -Sd.15 -R -J -O -K -W.1p,black -Cseis.cpt -t10 >> $out
awk '{print $1, $2, $3}' $picks_dir/unc_moho_depths_Cross-section_22.txt | gmt psxy -i0,1,2 -Sd.15 -R -J -O -K -W.1p,black -Cseis.cpt -t10 >> $out
awk '{print $1, $2, $3}' $picks_dir/unc_moho_depths_Cross-section_23.txt | gmt psxy -i0,1,2 -Sd.15 -R -J -O -K -W.1p,black -Cseis.cpt -t10 >> $out
awk '{print $1, $2, $3}' $picks_dir/unc_moho_depths_Cross-section_24.txt | gmt psxy -i0,1,2 -Sd.15 -R -J -O -K -W.1p,black -Cseis.cpt -t10 >> $out
awk '{print $1, $2, $3}' $picks_dir/unc_moho_depths_Cross-section_25.txt | gmt psxy -i0,1,2 -Sd.15 -R -J -O -K -W.1p,black -Cseis.cpt -t10 >> $out
awk '{print $1, $2, $3}' $picks_dir/unc_moho_depths_Cross-section_26.txt | gmt psxy -i0,1,2 -Sd.15 -R -J -O -K -W.1p,black -Cseis.cpt -t10 >> $out
awk '{print $1, $2, $3}' $picks_dir/unc_moho_depths_Cross-section_27.txt | gmt psxy -i0,1,2 -Sd.15 -R -J -O -K -W.1p,black -Cseis.cpt -t10 >> $out
awk '{print $1, $2, $3}' $picks_dir/unc_moho_depths_Cross-section_28.txt | gmt psxy -i0,1,2 -Sd.15 -R -J -O -K -W.1p,black -Cseis.cpt -t10 >> $out
awk '{print $1, $2, $3}' $picks_dir/unc_moho_depths_Cross-section_29.txt | gmt psxy -i0,1,2 -Sd.15 -R -J -O -K -W.1p,black -Cseis.cpt -t10 >> $out
awk '{print $1, $2, $3}' $picks_dir/unc_moho_depths_Cross-section_30.txt | gmt psxy -i0,1,2 -Sd.15 -R -J -O -K -W.1p,black -Cseis.cpt -t10 >> $out
awk '{print $1, $2, $3}' $picks_dir/unc_moho_depths_Cross-section_31.txt | gmt psxy -i0,1,2 -Sd.15 -R -J -O -K -W.1p,black -Cseis.cpt -t10 >> $out



awk '{print $1, $2, $3}' $picks_dir/moho_depths_Cross-section_21.txt | gmt psxy -i0,1,2 -Sc.16 -R -J -O -K  -Cseis.cpt -t5 >> $out
awk '{print $1, $2, $3}' $picks_dir/moho_depths_Cross-section_22.txt | gmt psxy -i0,1,2 -Sc.16 -R -J -O -K  -Cseis.cpt -t5 >> $out
awk '{print $1, $2, $3}' $picks_dir/moho_depths_Cross-section_23.txt | gmt psxy -i0,1,2 -Sc.16 -R -J -O -K  -Cseis.cpt -t5 >> $out
awk '{print $1, $2, $3}' $picks_dir/moho_depths_Cross-section_24.txt | gmt psxy -i0,1,2 -Sc.16 -R -J -O -K  -Cseis.cpt -t5 >> $out
awk '{print $1, $2, $3}' $picks_dir/moho_depths_Cross-section_25.txt | gmt psxy -i0,1,2 -Sc.16 -R -J -O -K  -Cseis.cpt -t5 >> $out
awk '{print $1, $2, $3}' $picks_dir/moho_depths_Cross-section_26.txt | gmt psxy -i0,1,2 -Sc.16 -R -J -O -K  -Cseis.cpt -t5 >> $out
awk '{print $1, $2, $3}' $picks_dir/moho_depths_Cross-section_27.txt | gmt psxy -i0,1,2 -Sc.16 -R -J -O -K  -Cseis.cpt -t5 >> $out
awk '{print $1, $2, $3}' $picks_dir/moho_depths_Cross-section_28.txt | gmt psxy -i0,1,2 -Sc.16 -R -J -O -K  -Cseis.cpt -t5 >> $out
awk '{print $1, $2, $3}' $picks_dir/moho_depths_Cross-section_29.txt | gmt psxy -i0,1,2 -Sc.16 -R -J -O -K  -Cseis.cpt -t5 >> $out
awk '{print $1, $2, $3}' $picks_dir/moho_depths_Cross-section_30.txt | gmt psxy -i0,1,2 -Sc.16 -R -J -O -K  -Cseis.cpt -t5 >> $out
awk '{print $1, $2, $3}' $picks_dir/moho_depths_Cross-section_31.txt | gmt psxy -i0,1,2 -Sc.16 -R -J -O -K  -Cseis.cpt -t5 >> $out
awk '{print $1, $2, $3}' $picks_dir/moho_depths_Cross-section_1.txt | gmt psxy -i0,1,2 -Sc.16 -R -J -O -K  -Cseis.cpt -t5 >> $out
awk '{print $1, $2, $3}' $picks_dir/moho_depths_Cross-section_2.txt | gmt psxy -i0,1,2 -Sc.16 -R -J -O -K  -Cseis.cpt -t5 >> $out
awk '{print $1, $2, $3}' $picks_dir/moho_depths_Cross-section_3.txt | gmt psxy -i0,1,2 -Sc.16 -R -J -O -K  -Cseis.cpt -t5 >> $out
awk '{print $1, $2, $3}' $picks_dir/moho_depths_Cross-section_4.txt | gmt psxy -i0,1,2 -Sc.16 -R -J -O -K  -Cseis.cpt -t5 >> $out
awk '{print $1, $2, $3}' $picks_dir/moho_depths_Cross-section_5.txt | gmt psxy -i0,1,2 -Sc.16 -R -J -O -K  -Cseis.cpt -t5 >> $out
awk '{print $1, $2, $3}' $picks_dir/moho_depths_Cross-section_6.txt | gmt psxy -i0,1,2 -Sc.16 -R -J -O -K  -Cseis.cpt -t5 >> $out
awk '{print $1, $2, $3}' $picks_dir/moho_depths_Cross-section_7.txt | gmt psxy -i0,1,2 -Sc.16 -R -J -O -K  -Cseis.cpt -t5 >> $out
awk '{print $1, $2, $3}' $picks_dir/moho_depths_Cross-section_8.txt | gmt psxy -i0,1,2 -Sc.16 -R -J -O -K  -Cseis.cpt -t5 >> $out
awk '{print $1, $2, $3}' $picks_dir/moho_depths_Cross-section_9.txt | gmt psxy -i0,1,2 -Sc.16 -R -J -O -K  -Cseis.cpt -t5 >> $out
awk '{print $1, $2, $3}' $picks_dir/moho_depths_Cross-section_10.txt | gmt psxy -i0,1,2 -Sc.16 -R -J -O -K  -Cseis.cpt -t5 >> $out
awk '{print $1, $2, $3}' $picks_dir/moho_depths_Cross-section_11.txt | gmt psxy -i0,1,2 -Sc.16 -R -J -O -K  -Cseis.cpt -t5 >> $out
awk '{print $1, $2, $3}' $picks_dir/moho_depths_Cross-section_12.txt | gmt psxy -i0,1,2 -Sc.16 -R -J -O -K  -Cseis.cpt -t5 >> $out
awk '{print $1, $2, $3}' $picks_dir/moho_depths_Cross-section_13.txt | gmt psxy -i0,1,2 -Sc.16 -R -J -O -K  -Cseis.cpt -t5 >> $out
awk '{print $1, $2, $3}' $picks_dir/moho_depths_Cross-section_14.txt | gmt psxy -i0,1,2 -Sc.16 -R -J -O -K  -Cseis.cpt -t5 >> $out
awk '{print $1, $2, $3}' $picks_dir/moho_depths_Cross-section_15.txt | gmt psxy -i0,1,2 -Sc.16 -R -J -O -K  -Cseis.cpt -t5 >> $out
awk '{print $1, $2, $3}' $picks_dir/moho_depths_Cross-section_16.txt | gmt psxy -i0,1,2 -Sc.16 -R -J -O -K  -Cseis.cpt -t5 >> $out
awk '{print $1, $2, $3}' $picks_dir/moho_depths_Cross-section_17.txt | gmt psxy -i0,1,2 -Sc.16 -R -J -O -K  -Cseis.cpt -t5 >> $out
awk '{print $1, $2, $3}' $picks_dir/moho_depths_Cross-section_18.txt | gmt psxy -i0,1,2 -Sc.16 -R -J -O -K  -Cseis.cpt -t5 >> $out
awk '{print $1, $2, $3}' $picks_dir/moho_depths_Cross-section_19.txt | gmt psxy -i0,1,2 -Sc.16 -R -J -O -K  -Cseis.cpt -t5 >> $out
awk '{print $1, $2, $3}' $picks_dir/moho_depths_Cross-section_0.txt | gmt psxy -i0,1,2 -Sc.16 -R -J -O -K  -Cseis.cpt -t5 >> $out
awk '{print $1, $2, $3}' $picks_dir/moho_depths_Cross-section_20.txt | gmt psxy -i0,1,2 -Sc.16 -R -J -O -K  -Cseis.cpt -t5 >> $out

#

#







# -=================================================================================================================- #
# ------------------------------------------------------------------------------------------------------------------- #
echo Create legend...
echo Create legend...
gmt pslegend <<END -R -J -Dx3.7i/0.3i+w0i/0.0i/TC -C0.1i/0.1i -F+gwhite+pthin -P -O -K --FONT_ANNOT_PRIMARY=6p >> $out
G -0.05i
H 7 Manual picks
D0.1i 0.5p
G .04i
S .05i c .08i seagreen2 0.7p,seagreen2 0.18i Moho depth
G .06i
S .05i d .06i seagreen2 0.7p,black 0.18i Unc. Moho depth
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


# Cross-section numbers
gmt pstext -R -J -O -K  -F+f5,Helvetica,black+jBL+a0 -Gwhite >> $out << END
5 42 2
6 42 3
7 42 4
8 42 5
9 42 6
10 42 7
11 42 8
12 42 9
13 42 10
14 42 11
15 42 12
16 42 13
17 42 14
18 45 15
19 45 16
20 45 17
21 45 18
22 45 19
23 45 20
4 42 1
3 42 0
1.5 43 21
1.5 43.7 22
1.5 44.4 23
1.5 45.1 24
1.5 45.8 25
1.5 46.5 26
1.5 47.2 27
1.5 47.9 28
1.5 48.6 29
1.5 49.3 30
1.5 50.0 31
END
gmt pscoast -W1/0.05 -Df $proj -R$west/$east/$south/$north -K -B5WSen -O -P -Y-9 >> $out
gmt grdimage -R -J $topodir/ETOPO1_Bed_g_gmt4.grd -Cmy_topo.cpt -O -K >> $out
gmt pscoast -W1/0.05 -Df -J -R -K -O -P -Sazure1 -N1/0.05p,black -L3.4/49.7/48/200+l+u >> $out
gmt makecpt -Cviridis -T20/60/5 -D+i -I > seis.cpt

#gmt makecpt -C../files/bamako.cpt -T20/80/5 -D+i -I > seis.cpt

#http://gmt.soest.hawaii.edu/doc/5.3.2/grdcontour.html
#gmt surface table_5.11 -R -I0.2 -Graws5.nc -T0.5
#gmt grdview raws5.nc -R -J -B -Cex16.cpt -Qs -O -K -Y-3.75i -X-3.5i >> $ps
#echo "3.25 7 surface (tension = 0.5)" | gmt pstext -R -J -O -K -N -F+f18p,Times-Roman+jCB >> $out
#
awk '{print $1, $2, $3}' $picks_dir/moho_depths_all.dat | gmt psxy -i0,1,2 -Sc.05 -R -J -O -K -Cseis.cpt -t5 >> $out

#grid='-I150+k'
##grid='-I1.5'
##awk '{print $1, $2, $3}' $picks_dir/moho_depths_all.dat | gmt psxy -i0,1,2 -Sc.05 -R -J -O -K -W -Cseis.cpt -t5 >> $out
#gmt blockmean $picks_dir/moho_depths_all.dat -R $grid > mean.xyz
#Block average (x,y,z) data tables by L2 norm
#gmt surface mean.xyz -R $grid -T0.3 -Gdata.nc
##surface reads randomly-spaced (x,y,z) triples from standard input [or table] a
## nd produces a binary grid file of gridded values z(x,y) by solving:
##(1 - T) * L (L (z)) + T * L (z) = 0
## gmt grdcontour data.nc -J -B -C1 -A2  -Gd5c -S0.1 -O -K -L20/80 -Wathin,black -Wcthinner,gray30 >> $out
#gmt grdcontour data.nc -J -B -Cseis.cpt -A10+f7p+o -Gd5c -S15 -O -K -L20/60 -Wathin+cl -Wcthin+c >> $out
# -A
#gmt grdview data.nc -R -J -B -Qs -O -K -Cseis.cpt -S5 -Wc.1 >> $out
#gmt grdcontour data.nc -J -B -Cseis.cpt -N -A10+f7p+o -Q15 -Gd5c -S5 -O -K -Lp -W0.1p >> $out


grid='-I70+k'
gmt blockmean ../files/moho_picks/moho_depths_all.dat -R $grid -Sm -V -C > ship_50k.nc
#gmt nearneighbor $region $grid -S40k -Gship_50k.nc ship.b -bi
#gmt grdcontour ship_50k.nc -J -P -B -Gd5c -O -K >> $out
gmt surface ship_50k.nc -R $grid -T0.4 -Gdata.nc
gmt psmask -R $grid ship_50k.nc -J -K -O -Gwhite >> $out
#gmt grdcontour data.nc -J -B -Cseis.cpt -A10+f7p+o -Gd5c -S15 -O -K -Lp -Wathin+cl -Wcthin+c -R$west/$east/$south/$north >> $out
##alt
gmt grdcontour data.nc -J -B -Cseis.cpt -A10+f6p+o -Gd8c -S5 -Nseis.cpt -O -K -Lp -W0.2p,gray15 -R$west/$east/$south/$north >> $out
#gmt grdcontour data.nc -J -B -Cseis.cpt -A- -S5 -N -O -K -Lp -Wathin+cl -Wcthin+c -R$west/$east/$south/$north >> $out
gmt psmask -C -O -K >> $out



#gmt grdcontour data.nc -J -B -Cseis.cpt -Nseis.cpt -A10+f7p+o -Gd5c -S15 -O -K -L20/60 -Wcthinner,black >> $out
gmt pscoast -W1/0.05 -Df -J -R -K -O -P -N1/0.05p,black -L3.4/49.7/48/200+l+u >> $out
#awk '{print $1, $2, $3}' $picks_dir/moho_depths_all.dat | gmt psxy -i0,1,2 -Sc.05 -R -J -O -K -W -Cseis.cpt -t5 >> $out
#awk '{print $1, $2, $3}' $picks_dir/moho_depths_all.dat | gmt psxy -i0,1,2 -Sc.05 -R -J -O -K -Cseis.cpt -t5 >> $out
gmt psscale -Dx9.2/1.9+o0/0i+w1.0i/0.08i+h+e -R -J -F+gwhite+p1p -Cseis.cpt -Bxa10f5 -Bx+l"Moho depth (km)" -O -K  >> $out




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

# ------------------------------------------------------------------------------------------------------------------- #
gmt psxy -R -J -T -O >> $out
gmt psconvert -Tf -A $out
evince ${out%.*}.pdf
