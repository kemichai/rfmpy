"""
Little script that tests that the re-orientation thingy
we use works fine.

|--------------------------|
|............... az, dip ..|
|vert  component:  0, -90  |
|north component:  0,   0  |
|east  component: 90,   0  |
|--------------------------|

Location: Chavannes-pres-renens, CH
Date: Feb 2022
Author: Konstantinos Michailos
"""
from obspy.signal.rotate import rotate2zne
from obspy import Stream, read

st = read('/home/kmichall/Desktop/Codes/github/rfmpy/rfmpy/tests/orient_horizontal_components/*SAC')

tr_e = st[1].copy()
tr_n = st[0].copy()
tr_z = st[2].copy()
tr_z.data, tr_n.data, tr_e.data = rotate2zne(st[2].data, 0, -90,
                                             st[0].data, 0, 0,
                                             st[1].data, 90, 0, inverse=False)
rot_stream = Stream()
rot_stream.append(tr_e)
rot_stream.append(tr_n)
rot_stream.append(tr_z)
all = Stream()
all = st + rot_stream
# Plotting both streams we get the same traces... for the configuration above
all.plot()
