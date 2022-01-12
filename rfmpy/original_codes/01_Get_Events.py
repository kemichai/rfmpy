import obspy
import sys
from obspy.clients.fdsn import Client
from obspy.geodetics.base import gps2dist_azimuth as gps2dist
import glob
import RF_Util as rfu
from pathlib import Path

# CREATE A SEISMIC EVENTS DATABASE
# Matteo Scarponi 30.11.2021

# Many of the functions called in this main code, are imported from RF_Util.py which provides a sequence of function definitions.

# 1- This code reads or downloads a catalog of seismic events according to the given parameters
# 2- Goes to the seismic raw database, where data is supposed to be stored in the form of daily files
# 3- Reads the daily file, cuts around the expected earthquake arrival time, plus minus a certain window as indicated in the parameters
# 4- Removes instrument response, filters, resamples and save the Z,E,N components to the EVENTS database folder
# The next code 02_Compute_RFs.py will read into the EVENTS database folder and compute RFs accordingly

#####################
# SET UP PARAMETERS #
#####################

path='/Volumes/Elements/working/DOWNLOAD/waveforms/'	# Path to raw data in form of daily files in the following format:  
														# ../NETWORK/STATION/STATION.NETWORK.AREA.CHANNEL.YEAR.JULDAY
pathout='/Volumes/Elements/working/EVENTS/' 			# Where to store event traces
pathmetadata='/Users/mscarpon/Desktop/Projects/Prague/ambient_noise/working/stations_metadata/'	# Path to .xml metadata files for each station

eps=sys.float_info.epsilon
start=obspy.UTCDateTime('2018-06-01T00:00:00.000')
end=obspy.UTCDateTime('2018-08-31T23:59:59.999')

lat0 =45.8	# Avarage profile lat
lon0 =8.3 	# Avarage profile lon
mw  =6 		# Minimum magnitude for event selection
rmin=28 	# Minimum distance for event selection
rmax=95 	# Maximum distance for event selection

fs=20   	 # Downsampling frequency (Hz)
freqmin=0.1  # minimum frequency for bandpass filtering (Hz)
freqmax=2    # maximum frequency for bandpass filtering (Hz)
corners=2    # Filter order

if fs<=2*freqmax:
	print('HARRY NYQUIST is NOT happy')
	quit()

tbefore = 60   # Time before direct P-wave arrival [s]
tafter  = 120  # Time after direct P-wave arrival [s]

networks=['XK','IV']		# Prepare a list of the networks you want to use

net={				

	'XK':['IA01A','IA02A','IA03A','IA04A',
	'IA05A','IA06A','IA07A','IA08A',
	'IA09A','IA10A'],

	'IV':['SATI','VARE']}

print(start)
print(end)
print('Sampling rate: ',fs)
print('Maximum freq: ',freqmax)
print('Minimum freq: ',freqmin)
print('Read time before P-arrival ',str(tbefore)+'s')
print('Read time after P-arrival ',str(tafter)+'s')
print('XK ',net['XK'])
print('IV ',net['IV'])

##############
# ---------- #
# -- MAIN -- #
# ---------- #
##############

# DOWNLOAD OR READ EARTHQUAKE EVENT CATALOG

# client  = Client('IRIS')
# catalog = client.get_events(
# 	starttime=start,
# 	endtime=end,
# 	minmagnitude=mw,
# 	latitude=lat0,
# 	longitude=lon0,
# 	minradius=rmin,
# 	maxradius=rmax,
# 	magnitudetype='mww')

catalog=obspy.core.event.read_events(
	'/Volumes/Elements/working/IRIS_catalog.xml')

print(catalog.__str__(print_all=True))

# Prepare dictionary for coordinates of all seismic stations, for all networks

lon,ele,lat=[],[],[]
allstations=[]
for network in networks:
	for station in net[network]:
		if station not in allstations:
			inv=obspy.read_inventory(pathmetadata+'*'+station+'*xml')
			lon.append(inv[0][0].longitude)
			lat.append(inv[0][0].latitude)
			ele.append(inv[0][0].elevation)
			allstations.append(station)
londict=dict(zip(allstations,lon))
latdict=dict(zip(allstations,lat))
eledict=dict(zip(allstations,ele))

# READ, PICK AND TRIM EVENT TRACES

for event in catalog:

	print(event)

	eventdate=event.origins[0].time

	# for network in networks:
	for network in ['IV']:

		# for station in net[network]:
		for station in ['VARE']:

			print(station)
			print(londict[station],latdict[station])

			# Reads and trims the raw seismic data around the earthquake arrival time for E.N.Z components
			# See function Pick in RF_Utils.py for more details

			stream,baz,az,dist=rfu.Pick(
				event=event,
				station=station,
				network=network,
				path=path,
				latdict=latdict,
				londict=londict,
				tbefore=tbefore,
				tafter=tafter)

			if not stream:

				continue

			elif len(stream)==3:

				# Remove instrument response
				# Currently station metadata is downloaded from EIDA

				# stream=rfu.RemoveResponse(stream)

				# if not stream:

				# 	continue

				# Bandpass filter and resample

				stream.detrend('demean')
				stream.detrend('linear')

				stream.taper(
					max_percentage=0.1,
					type='cosine',
					max_length=5)

				stream.filter(
					'bandpass',
					freqmin=freqmin,
					freqmax=freqmax,
					corners=corners,
					zerophase=True)

				stream.taper(
					max_percentage=0.1,
					type='cosine',
					max_length=5)

				# stream.resample(
				# 	fs,
				# 	window='hanning',
				# 	no_filter=True)

				stream.plot()

				# Save Event seismic recording E N Z components

				stla=latdict[station]
				stlo=londict[station]
				stel=eledict[station]

				info=stla,stlo,stel,baz,az,dist,tbefore

				# TO BE DONE: Browse catalog to look for events which happened in the same second
				# and create folders accordingly. If not the case, then:

				EventFolder=pathout\
					+str(eventdate.year)+'.'\
					+str(eventdate.julday)+'.'\
					+str(eventdate.hour)+'.'\
					+str(eventdate.minute)+'.'\
					+str(eventdate.second)+'/'

				# Save Event traces to the catalog

				Path(EventFolder).mkdir(exist_ok=True)

				rfu.SaveEventTrace(
					trace=stream[0],
					event=event,
					info=info,
					tbefore=tbefore,
					pathout=EventFolder)

				rfu.SaveEventTrace(
					trace=stream[1],
					event=event,
					info=info,
					tbefore=tbefore,
					pathout=EventFolder)

				rfu.SaveEventTrace(
					trace=stream[2],
					event=event,
					info=info,
					tbefore=tbefore,
					pathout=EventFolder)

			else:

				print('! PROBLEM with - Pick - function')
				quit()