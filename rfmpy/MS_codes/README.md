### Matteo's codes 

### See notes: 

- 01_Get_Events.py: provided an earthquake 
catalog, browses a raw seismic database 
(daily files for stations and components), 
to find earthquake arrivals and trim traces
 around the time window of interest.
  Each event is then saved for later 
  RF computation

- 02_Compute_RFs.py: provided a seismic catalog
 with events (with single components traces already 
 cut around a teleseismic earthquake arrival), reads
  the traces and computes RFs. In this code, I read 
  the recordings for the all seismic stations of a given network, 
  so that quality control (QC) is applied both on single
   stations and on single stations compared to median of 
   the network (for QC practices see Subedi et al. 2018 
   and Hetényi 2007 - thesis). RFs are then saved after computation.

- RF_Util.py: this contains a list of functions definition
 which are called by the two codes above. The Time-domain 
 iterative deconvolution method (Ligorría and Ammon 1999) 
 is implemented in function IterativeRF

- 03_SynRF.py: this is just wrapping the fortran code Raysum 
from Frederiksen and Bostock 1999 to produce synthetic R,T,Z 
for a given model geometry. I am using this to check that, for 
a simple model setup, the function IterativeRF is working correctly.
