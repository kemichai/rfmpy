import obspy
import glob
from obspy.taup import TauPyModel
import numpy as np
import platform
import core.utils.RF_Util as rf_util
from core.utils.signal_processing import rotate_trace, remove_response, ConvGauss
from core.utils.qc import rms_quality_control, rf_quality_control, RFQuality

# Set up parameters and paths
if platform.node().startswith('kmichailos-laptop'):
    data_root_dir = '/media/kmichailos/SEISMIC_DATA/Data_archive'
    codes_root_dir = '/home/kmichailos/Desktop/codes/bitbucket'
    desktop_dir = '/home/kmichailos/Desktop'
else:
    data_root_dir = '/media/kmichall/SEISMIC_DATA/Data_archive'
    codes_root_dir = '/home/kmichall/Desktop/Codes/bitbucket'
    desktop_dir = '/home/kmichall/Desktop'

# COMPUTE SEISMIC RECEIVER FUNCTIONS
# Matteo Scarponi 30.11.2021

# This code reads a database of seismic EVENTS
# EVENTS are stored in the following format ../EVENTS/YEAR.JULDAY.HOUR.MINUTE.SECOND/YEAR.JULDAY.HOUR.MINUTE.SECOND.STATION.COMPONENT.SAC

##############
# PARAMETERS #
##############

path=desktop_dir + '/RF_test/test_data/'
pathout=desktop_dir + '/RF_test/RF/'

iterations=100							# Iterative time-domain RFs computation

freqmax=2 								# SAME as used in Get_Events.py

C1,C2,C3,C4=10,10,1,1 					# Quality control parameters for Z,E,N traces

model=TauPyModel(model='iasp91')

# Define the list of the networks you want to use: in this case ['XK']
# and a dictionary which lists all the stations belonging to a certain network
# In this case for IvreaArray: stations IA##A belonging to network XK
# This is useful as some quality controls concern the single stations while others compare station with respect to the median of the network

networks=['CH','IV']		# Prepare a list of the networks you want to use

net={				

	'XK':['IA01A','IA02A','IA03A','IA04A',
	'IA05A','IA06A','IA07A','IA08A',
	'IA09A','IA10A'],

	'IV':['SATI','VARE']}

channels=['BH'] # Possible channels for seismic data recordings HH BH ... etc etc

##############
# -- MAIN -- #
##############
# Get list of events - as prepared by 01_Get_Events.py

events = glob.glob(path + '*')

for event in events:

	print(event)

	# for network in networks:
	for network in ['CH']:

		Etraces = obspy.Stream()
		Ntraces = obspy.Stream()
		Ztraces = obspy.Stream()

		# For a given network, collect data from all the stations which recorded the earthquake event

		# for station in net[network]:
		for station in ['DAVOX']:

			Ztrace = obspy.read(event + '/*' + station + '*' + channels[0] + '*Z*')
			Etrace = obspy.read(event + '/*' + station + '*' + channels[0] + '*E*')
			Ntrace = obspy.read(event + '/*' + station + '*' + channels[0] + '*N*')
			print(Ztrace)
			if len(Ztrace) < 1:
				# Component not available for this station - event combination
				continue
			elif len(Ztrace) > 1:
				print('! PROBLEM AMBIGUOUS TRACES !')
				quit()
			if len(Etrace) < 1:
				# Component not available for this station - event combination
				continue
			elif len(Etrace) > 1:
				print('! PROBLEM AMBIGUOUS TRACES !')
				quit()
			if len(Ntrace) < 1:
				# Component not available for this station - event combination
				continue
			elif len(Ntrace) > 1:
				print('! PROBLEM AMBIGUOUS TRACES !')
				quit()

			Ztraces.append(Ztrace[0])
			Etraces.append(Etrace[0])
			Ntraces.append(Ntrace[0])

		fullrmsZ = []
		fullrmsE = []
		fullrmsN = []

		rmsbgZ = []
		maxbgZ = []
		maxpkZ = []

		tbefore = Ztrace[
			0].stats.sac.a  # Time before P-arrivaltime - should be same for all traces! (or move it in the loop)
		fs = Ztraces[
			0].stats.sampling_rate  # Sampling rate [Hz]		- should be same for all traces! (or move it in the loop)
		delta = Ztraces[0].stats.delta  # Delta [s] 				- should be same for all traces! (or move it in the loop)

		# Loop on the stations to check quality control
		# of the single station and of the station compared to
		# the other stations of the network

		# Compute signal and rms background for all the components for each station
		# Some of this is used for single station quality control
		# Some of this is usef for station quality control in comparison with the median of the network
		# If only one station in being used, comparison with the network should be successful automatically
		# For further explanation on this quality control, please check Subedi et al. 2018 and Hetenyi 2007 (PhD thesis)

		for i in range(len(Ztraces)):
			fullrmsZ.append(np.sqrt(np.mean(Ztraces[i].data ** 2)))
			fullrmsE.append(np.sqrt(np.mean(Etraces[i].data ** 2)))
			fullrmsN.append(np.sqrt(np.mean(Ntraces[i].data ** 2)))

			i0 = int((tbefore - 30) * fs)
			i1 = int((tbefore - 5) * fs)
			i2 = int((tbefore + 20) * fs)
			rmsbgZ.append(np.sqrt(np.mean((Ztraces[i].data[i0:i1 + 1] ** 2))))
			maxbgZ.append(np.max(Ztraces[i].data[i0:i1 + 1]))
			maxpkZ.append(np.max(Ztraces[i].data[i1:i2 + 1]))

		medianZ = np.median(fullrmsZ)
		medianE = np.median(fullrmsE)
		medianN = np.median(fullrmsN)

		##############################
		# QUALITY CONTROL CONDITIONS #
		##############################

		for i in range(len(Ztraces)):

			E1 = (medianE * C1 >= fullrmsE[i]) and (fullrmsE[i] >= medianE / C2)
			N1 = (medianN * C1 >= fullrmsN[i]) and (fullrmsN[i] >= medianN / C2)
			Z1 = (medianZ * C1 >= fullrmsZ[i]) and (fullrmsZ[i] >= medianZ / C2)

			Z2 = (maxpkZ[i] >= maxbgZ[i] * C3)
			Z3 = (maxpkZ[i] >= rmsbgZ[i] * C4 * np.sqrt(2))

			# If quality control condtions are met on N Z E traces
			# Roration to TRZ is performed and RF is computed

			if E1 and N1 and Z1 and Z2 and Z3:

				# If quality control is successfull

				Z = Ztraces[i].copy()
				T = Etraces[i].copy()
				R = Ntraces[i].copy()

				# Rotate to Vertical (Z), Radial (R) and Tangential (T) components
				# Consists of a horizontal rotation around Z to align along the baz direction
				# Plus sign reversal of the R component so that main peak looks positive

				T.data, R.data, Z.data = rfu.RoTrace(
					E=Etraces[i].data,
					N=Ntraces[i].data,
					Z=Ztraces[i].data,
					baz=Ztraces[i].stats.sac.baz)

				Z.stats.channel = 'HHZ'
				T.stats.channel = 'HHT'
				R.stats.channel = 'HHR'

				# Ready for processing

				processR = R.copy()
				processZ = Z.copy()

				RF = processR.copy()
				RF.stats.channel = 'RRF'

				RF.data, ds = rfu.IterativeRF(
					traceZ=processZ,
					traceR=processR,
					iterations=iterations,
					FlagStages=True,
					FlagSummary=True)

				RFconvolve = RF.copy()
				RFconvolve = rfu.ConvGauss(
					SpikeTrace=RFconvolve,
					freqmax=freqmax,
					delta=RFconvolve.stats.delta)

				RFconvolve.stats.sac.a = ds

				# RF quality control
				# It checks the amplitude and the location of the direct-P peak:
				# It should be around 0 or "ds" if ds != 0 and between 0.1 and 1 in terms of amplitude

				time1, time2, amp1, amp2 = rfu.RFQuality(RFconvolve)

				if time1 and time2 and amp1 and amp2:

					# If the quality control on the Radial RF is successful
					# Tangential component RF is computed and they are saved

					processZ = Z.copy()
					processT = T.copy()

					RFconvolve.plot()

					TRF = processT.copy()
					TRF.stats.channel = 'TRF'

					TRF.data, ds = rfu.IterativeRF(
						traceZ=processZ,
						traceR=processT,
						iterations=iterations,
						FlagStages=False,
						FlagSummary=False)

					# IterativeRF provides a serie of spikes
					# Here the serie of spikes is convolved by a gaussian bell
					# whose width should match the highest frequency kept in the data

					TRFconvolve = TRF.copy()
					TRFconvolve = rfu.ConvGauss(
						SpikeTrace=TRFconvolve,
						freqmax=freqmax,
						delta=TRFconvolve.stats.delta)

				# And saved

				# rfu.RFSave(
				# 	Trace=RFconvolve,
				# 	pathOUT=pathout)

				# rfu.RFSave(
				# 	Trace=TRFconvolve,
				# 	pathOUT=pathout)

				else:

					# Bad RF trace
					continue














