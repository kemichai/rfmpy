# README #

### Download data
In the archive get_data.py is a pair of functions that
act as useful wrappers for IRIS downloads.
They will download data and store it in day-long miniseed
format in directories of the form:

Yyyyy/Rjjj.01

Where yyyy is the year and jjj is the julian day.

Files will be named network.station.location.channel.year.julianday

You can either write your own scripts
around these (e.g. if you want multiple networks),
or just call it from either the command line
(usage: python get_data.py network starttime endtime,
where starttime and endtime are UTCDateTime parsable strings
such as "2009-12-31T00:00:00.0"), or from a python session.


P.s. You can change the formatting of the output file names if you want.
 This is controlled by the nested get_mseed_storage function.



