gmt begin tele
	gmt makecpt -Chot -T0/200/50 -H > n.cpt

	# first do an overhead of the east coast from 160 km altitude point straight down
#	gmt coast -R-180/-20/0/90 -JPoly/4i -Bx30g10 -By10g10 -Dc -A1000 -Glightgray -Wthinnest
#	gmt coast -Rd -JG10/41/4.5i -Bg -Dc -Gwhite -Slightblue -Wthinnest -A10000
	gmt coast -Rg -JE10/40/110/5i -Bg  -BWsNE -Dc -A10000 -Glightgray -Wthinnest

#	gmt grdimage /home/kmichall/Desktop/topo/ETOPO1_Bed_g_gmt4.grd -Cmy_topo.cpt

#	gmt legend -DjLM+w3.5c+jRM+o1c/0 -F+p+i <<- EOF
#	H 12 Legend
#	D 0.25c 1p
#	S 0.4c s 0.3c blue       0.25p 0.75c Ocean
#	S 0.4c s 0.3c lightblue  0.25p 0.75c Ice front
#	S 0.4c s 0.3c lightbrown 0.25p 0.75c Grounding line
#	EOF
#	gmt plot tele_events_sample.txt -h1 -Scc -i1,2,3,4+s0.025 -Gred -Wthinnest

  awk '{print $1, $2, $3/1000}' tele_events_sample.txt | gmt psxy -i0,1,2s0.2 -Sc.25 -Wthinnest -t5 -Cn.cpt
#	# Fancy line
	gmt plot -B -W2p,red <<- EOF
	0	45
	20 45
	20	50
	0	  50
	0 45
	EOF


gmt end show


