import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import obspy
from conv_gauss import conv_gauss

import subprocess
import obspy.io.sac.sactrace as sactrace
import matplotlib.pyplot as plt

from RF_Util import IterativeRF,ConvGauss

# PARAMETERS

iterations=30
freqmax=1.

#####################
# Frederiksen SynRF #
#####################

# Uses the code from Frederiksen to produce a synthetic Z,R,T seismogram
# for the given model in sample.mod
# and the given ray geometry in sample.geom

wd='/Users/mscarpon/Desktop/SynRF/Raysum/Sample/'

call='./seis-spread sample.mod sample.geom sample.ph sample.arr sample.tr P'
subprocess.call(call,cwd=wd,shell=True)

tracefile=wd+'/sample.tr'

# Raysum outputs Z,R,T components which are read, used for Radial RF computation and plotted here

########################
# Reading and Plotting #
########################

# Frederiksen

with open(tracefile,'r') as file:
	line1=file.readline()
	line2=file.readline()
	line2=line2.split()
	nTraces=int(line2[0])
	nSamples=int(line2[1])
	delta=float(line2[2])
	shift=float(line2[4])
	print('traces ',nTraces)
	print('samples ',nSamples)
	print('delta ',delta)
	print('shift ',shift)
	fsamp=1./delta
	R=np.zeros((nTraces,nSamples),dtype='float')
	T=np.zeros((nTraces,nSamples),dtype='float')
	Z=np.zeros((nTraces,nSamples),dtype='float')
	RF=np.zeros((nTraces,nSamples),dtype='float')
	for trace in range(nTraces):
		line1=file.readline().split('\n')[0]
		line2=file.readline().split('\n')[0]
		line3=file.readline().split('\n')[0]
		print(line2)
		for sample in range(nSamples):
			R[trace,sample],T[trace,sample],Z[trace,sample]=file.readline().split()

# for i in range(nTraces):
for i in [0]:
	Rtrace = obspy.Trace(data=R[i,:]) 
	Ztrace = obspy.Trace(data=Z[i,:])
	Ttrace = obspy.Trace(data=T[i,:])
	Rtrace.stats.sampling_rate=fsamp
	Ztrace.stats.sampling_rate=fsamp
	Ttrace.stats.sampling_rate=fsamp
	RFtrace,dt=IterativeRF(
		traceZ=Ztrace,
		traceR=Rtrace,
		iterations=iterations,
		FlagStages=False,
		FlagSummary=True)
	RFtrace=ConvGauss(
		SpikeTrace=RFtrace,
		freqmax=freqmax,
		delta=delta)
	RF=RFtrace.data

f=plt.figure()
ax=plt.subplot(111)
ax.grid(True,linestyle='--',lw=0.5,alpha=0.5)
ax.set_xlim([-2,25])
ax.set_xlabel('Time [s]')
ax.set_ylabel('Amplitude')
time=np.arange(len(RF))*delta-dt
ax.plot(time,RF,'b-',lw=0.75,label='Syn RF - Frederiksen')
plt.show()
plt.close()