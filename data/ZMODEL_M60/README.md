Read .dat file and stored the lists of Lon, Lat, dep, Vp and Vs within the following limits 
min_lon = 0
min_lat = 30
max_lon = 40
max_lat = 60
in a .npz file to avoid 
having larger than 50 MB files in the repository. 


Here is the code snippet used for this analysis.
```
import numpy as np

# read the .dat file
p_velocities = []
s_velocities = []
longitudes = []
latitudes = []
depths = []
with open(path_zmodel_m60 + 'ZMODEL_M60.dat', 'r') as f:
 for line in f:
     if line.startswith('#'):
         print('|Reading ZMODEL_M60 velocity model...              |')
         continue
     else:
         ln = line.split()
         lon_ = float(ln[0])
         lat_ = float(ln[1])
         dep_ = float(ln[2])
         vp = float(ln[3])
         vs = float(ln[5])
         if lon_ < 40 and lon_ > 0 and lat_ > 40 and lat_ < 60:
             longitudes.append(lon_)
             latitudes.append(lat_)
             depths.append(dep_)
             p_velocities.append(vp)
             s_velocities.append(vs)

dictionary = {"longitudes": longitudes, "latitudes": latitudes,
                  "depths": depths, "Vp": p_velocities, "Vs":s_velocities}
np.savez('parameters.npz', **dictionary)
```


September 2023
