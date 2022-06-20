gmt begin tele
	gmt makecpt -Chot -T0/300/50 -H > n.cpt

	# first do an overhead of the east coast from 160 km altitude point straight down
#	gmt coast -R-180/-20/0/90 -JPoly/4i -Bx30g10 -By10g10 -Dc -A1000 -Glightgray -Wthinnest
#	gmt coast -Rd -JG10/41/4.5i -Bg -Dc -Gwhite -Slightblue -Wthinnest -A10000
	gmt coast -Rg -JE10/40/110/4.4i -Bg  -BWsNE -Dc -A10000 -Glightgray -Wthinnest

#	gmt grdimage /home/kmichall/Desktop/topo/ETOPO1_Bed_g_gmt4.grd -Cmy_topo.cpt
#
#	gmt legend -DjLM+w3.5c+jRM+o2c/0 -F+p+i <<- EOF
#	H 12 Legend
#	D 0.25c 1p
#	S 0.4c s 0.3c blue       0.25p 0.75c Ocean
#	S 0.4c s 0.3c lightblue  0.25p 0.75c Ice front
#	S 0.4c s 0.3c lightbrown 0.25p 0.75c Grounding line
#	EOF
#	gmt plot tele_events_sample.txt -h1 -Scc -i1,2,3,4+s0.025 -Gred -Wthinnest
  gmt colorbar -DjLM+w5.0c+jRM+o1.5c+e -Cn.cpt -Bxa50+lDEPTH -By+lkm
#  gmt colorbar -DJRM+w6.5c/0.5c+o1c/0+mc -Cn.cpt -Bxa1000+lELEVATION -By+lm



  awk '{print $1, $2, $3}' tele_events_sample.txt | gmt psxy -i0,1,2s0.5 -Sc.15 -Wthinnest -t25 -Cn.cpt
#	# Fancy line
	gmt plot -B -W2.5p,red <<- EOF
	0	45
	20 45
	20	50
	0	  50
	0 45
	EOF

	gmt plot -R0/19/0/25 -Jx1c -B0 -W2p -X-6c -Y-13.5c <<- EOF
	3	13
	20 13
	20 33
	3  33
	3  13
	EOF
#	gmt text -F+f18p+jBL -Dj8p/0 <<- EOF
#	3.3 24 b)
#	EOF

gmt end show


