from shapely.geometry import Polygon, Point
import numpy as np

lon = []
lat = []
dep = []
with open('moho_depths_all.dat', 'r') as f:
        for line in f:
            if line.startswith('#'):
                print(line)
                continue
            else:
                ln = line.split(',')
                lon.append(float(ln[0]))
                lat.append(float(ln[1]))
                dep.append(float(ln[2]))

coords_adria = [[17.0169, 15.4473, 14.9241,14.0970,13.7595,11.2616,10.0464,
                8.6624, 7.4135, 7.3629,7.8186, 8.2405, 9.3713, 10.7722,13.4557,20],
          [47.1074,47.1069,46.9806,46.8961,46.7905,46.7694,46.3680, 46.3891,
          45.3116, 44.4665,44.1708,44.4454, 44.3187, 44.0651,  42.2694, 43]]
coords_liguria = [
    [8.0106 ,
7.7082  ,
7.7182  ,
7.4158  ,
7.4057  ,
8.0005  ,
8.0307  ,
8.0408  ,
8.2223  ,
8.2324  ,
9.4320  ,
9.4924  ,
9.7545  ,
9.7848  ,
10.2686 ,
10.2787 ,
10.9037 ,
10.8836 ,
11.4464 ,
12.5469 ,
12.5973 ,
13.6154 ,
13.6961 ,
13.9884 ,],
    [
41.3210,
41.6486,
41.8845,
42.2776,
42.6708,
43.6143,
43.7322,
44.3874,
44.3874,
44.4791,
44.4922,
44.3874,
44.4005,
44.3088,
44.3088,
44.2039,
44.2039,
44.0991,
44.1122,
43.1949,
43.2342,
42.1990,
42.2252,
41.9107,
    ]]




adria_latitudes = coords_adria[1]
adria_longitudes = coords_adria[0]
p_adria= Polygon((np.asarray(list(zip(adria_latitudes, adria_longitudes)))))
liguria_latitudes = coords_liguria[1]
liguria_longitudes = coords_liguria[0]
p_liguria= Polygon((np.asarray(list(zip(liguria_latitudes, liguria_longitudes)))))

a_lon = []
a_lat = []
a_dep = []
l_lon = []
l_lat = []
l_dep = []
e_lon = []
e_lat = []
e_dep = []
for i, _ in enumerate(lon):
    if p_adria.contains(Point(lat[i], lon[i])):
        a_lon.append(lon[i])
        a_lat.append(lat[i])
        a_dep.append(dep[i])
    elif p_liguria.contains(Point(lat[i], lon[i])):
        l_lon.append(lon[i])
        l_lat.append(lat[i])
        l_dep.append(dep[i])
    else:
        e_lon.append(lon[i])
        e_lat.append(lat[i])
        e_dep.append(dep[i])



for i, ln in enumerate(a_lon):
    with open('adria.txt', 'a') as of:
        of.write('{}, {}, {}\n'.format(a_lon[i], a_lat[i], a_dep[i]))

for i, ln in enumerate(l_lon):
    with open('liguria.txt', 'a') as of:
        of.write('{}, {}, {}\n'.format(l_lon[i], l_lat[i], l_dep[i]))

for i, ln in enumerate(e_lon):
    with open('europe.txt', 'a') as of:
        of.write('{}, {}, {}\n'.format(e_lon[i], e_lat[i], e_dep[i]))

