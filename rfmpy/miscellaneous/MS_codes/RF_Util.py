from obspy.taup import TauPyModel
from obspy.geodetics.base import gps2dist_azimuth as gps2dist
from obspy.geodetics.base import kilometers2degrees as km2deg
import os.path
import obspy
import numpy as np
from obspy.clients.fdsn import RoutingClient
import obspy.io.sac.sactrace as sac
from scipy import signal
import matplotlib.pyplot as plt

######################################
# See each topic indicated like this #
######################################

#############################################################
# 1 - Receiver function Time-Domain iterative deconvolution #
#############################################################

client=RoutingClient('eida-routing')

def RFSave(Trace,pathOUT):

    stla = Trace.stats.sac.stla
    stlo = Trace.stats.sac.stlo
    stel = Trace.stats.sac.stel        

    evla = Trace.stats.sac.evla
    evlo = Trace.stats.sac.evlo
    evdp = Trace.stats.sac.evdp

    az = Trace.stats.sac.az
    baz = Trace.stats.sac.baz
    dist = Trace.stats.sac.dist

    y = Trace.stats.sac.nzyear
    d = Trace.stats.sac.nzjday
    h = Trace.stats.sac.nzhour
    m = Trace.stats.sac.nzmin
    s = Trace.stats.sac.nzsec

    staz = Trace.stats.station

    header = {'knetwk' : 'IvreaArray','kcmpnm': Trace.stats.channel, 'kstnm': Trace.stats.station, 'stla': stla, 'stlo': stlo, 'stel': stel,
              'evla': evla, 'evlo': evlo, 'evdp': evdp, 'mag' : Trace.stats.sac.mag,
              'az': az, 'baz': baz, 'dist': dist, 'nzyear': y, 'a': Trace.stats.sac.a,
              'nzjday': d, 'nzhour': h, 'nzmin': m, 'nzsec': s,
              'delta': Trace.stats.sac.delta}

    RFfilename = pathOUT + str(y) + '.' + str(d) + '.' + str(h) + '.' + str(m) + '.' +str(s) + '.RF.' + str(staz) + '.'+Trace.stats.channel+'.SAC'

    RF_to_file = sac.SACTrace(data=Trace.data,**header) 
    RF_to_file.write(RFfilename)
    return RFfilename

def RFQuality(Trace):

    ds=Trace.stats.sac.a

    iRF_edge=round(Trace.stats.sampling_rate*60)
    iRF = np.argmax(np.abs(Trace.data[:iRF_edge]))

    time1 = not (iRF <= (ds-0.6)*Trace.stats.sampling_rate)
    time2 = not (iRF >= (ds+2.1)*Trace.stats.sampling_rate)
    amp1  = not (Trace.data[iRF] <= 0.01)
    amp2  = not (Trace.data[iRF] >= 1)

    return time1,time2,amp1,amp2

def ConvGauss(SpikeTrace,freqmax,delta):

    ##############################################################################################################################
    # Now we convolve the spike train with a gaussian filter whose width is equal to the maximum frequency content of the signal #
    ##############################################################################################################################

    sigma       = 1./(2*np.pi*freqmax)
    time        = np.arange(-sigma*5,sigma*5+delta,delta)
    gauss       = np.exp(-time*time/(2*sigma*sigma))

    SpikeTrace.data = signal.convolve(SpikeTrace.data,gauss,mode='same')

    return SpikeTrace

def IterativeRF(traceZ,traceR,iterations=30,FlagStages=False,FlagSummary=False):

    # Implementation of the Iterative deconvolution method
    # based on Ligorria & Ammon (BSSA,89-5,1999)

    fs=int(traceR.stats.sampling_rate)

    trZ=traceZ.data
    trR=traceR.data
    tbefore=traceR.stats.sac.a

    # Cutting first ds seconds from the Z trace to create delay between R and Z traces 
    # This means that the direct-P arrival (or main reference peak) in the RFs should be at exactly "ds" second from zero.

    ds = 5              # Cut ds seconds from Z so that direct P-arrival appears at t==ds in the final RF trace
    delay = ds*fs
    trZ = trZ[delay:]
    
    nz = len(trZ)
    nr = len(trR)

    if nz > nr:
        trZ = trZ[:nr]
    
    dirac_sum = np.zeros(nr)       # Prepare empty trace where to store cross-correlation maxima iteratively
    xcz = signal.correlate(trZ,trZ,mode='full')
    mxcz = np.max(np.abs(xcz))
    trH = trR
        
    rms = []

    k=(2*nr-1)//2    
    iteration=0
    while iteration < iterations:
        iteration+=1

        dirac  = np.zeros(nr)
        xcr    = signal.correlate(trZ,trH,mode='full',method='fft')
        ixcr   = np.argmax(abs(xcr[nr-nz:nr]))
        shift  = nz - ixcr - 1
        dirac[shift] = xcr[nr-nz+ixcr]/mxcz
             
        newH = signal.convolve(trZ,dirac)
        newH = newH[:nr]
        dirac_sum = dirac_sum + dirac

        conv = signal.convolve(trZ,dirac_sum)
        conv = conv[:len(trR)]
        diff = trR[:] - conv[:]
        diff = np.linalg.norm(diff)
        normConv = np.linalg.norm(conv)
        normR = np.linalg.norm(trR)
        diff = diff/(np.sqrt(normR*normConv))*100
        rms.append(diff)
        
        if FlagStages:
        
            f=plt.figure(1)

            ax = plt.subplot(311)
            ttH=np.arange(len(trH))/fs
            ttZ=np.arange(len(trZ))/fs+ds
            ax.vlines(shift/fs,ymin=0,ymax=np.max(trZ),color='g',ls='--',label='peak')
            ax.plot(ttZ,trZ,'k',label='trZ')
            ax.plot(ttH,trH,'r',label='trH')
            ax.set_ylabel('amplitude')
            ax.set_xlabel('time [s]')
            ax.set_title('iteration: '+str(iteration)+'\nshift: '+str(shift))
            ax.legend(loc='best')
            ax.grid(True)

            ax = plt.subplot(312)
            tt=np.arange(len(xcr))/fs
            ax.plot(tt,xcr,'k',label='xcorr(Z,H)')
            ax.vlines(x=(nr-nz+ixcr)/fs,ymin=0,ymax=xcr[nr-nz+ixcr],color='r',label='max')
            ax.vlines(x=k/fs,ymin=0,ymax=xcr[nr-nz+ixcr],color='g',label='0 shift point')
            ax.set_ylabel('amplitude')
            ax.set_xlabel('time [s]')
            ax.set_title('shift: '+str(shift))
            ax.legend(loc='best')
            ax.grid(True)

            ax = plt.subplot(313)
            tt=np.arange(len(dirac_sum))/fs
            ax.plot(tt,dirac_sum,'k',label='dirac_sum')
            ax.plot(tt,dirac,'r',label='new spike')
            ax.set_ylabel('amplitude')
            ax.set_xlabel('time [s]')
            ax.set_title('Updating RF')
            ax.legend(loc='best')
            ax.grid(True)
    
            plt.tight_layout()
            plt.show()
        
        trH = trH - newH

    if FlagSummary:
        
        f=plt.figure(2)

        ax = plt.subplot(311)
        tt=np.arange(len(dirac_sum))/fs
        ax.plot(tt,dirac_sum,'k',lw=0.5,label='computed RF')
        ax.fill_between(tt,y1=dirac_sum,y2=0,where=dirac_sum>0,color='r')
        ax.fill_between(tt,y1=dirac_sum,y2=0,where=dirac_sum<0,color='b')
        ax.set_title('Radial receiver function '+traceZ.stats.station)
        ax.set_ylabel('amplitude')
        ax.set_xlabel('time (s)')
        ax.grid(True,alpha=0.5,lw=0.2)
        ax.legend(loc='best')
        
        ax = plt.subplot(312)
        tt=np.arange(len(trR))/fs
        ax.plot(tt,trR,'k',label='R component')
        ax.plot(tt,conv,'r',label='R approx')
        ax.set_title(str(iteration)+'th iteration')
        ax.set_ylabel('amplitude [counts]')
        ax.set_xlabel('time (s)')
        ax.legend(loc='best')
        ax.grid(True)
        
        ax = plt.subplot(313)
        ax.plot(range(len(rms)),rms,'ro-',label='rms%')
        ax.set_title('Final rms: '+ str("%.2f" % rms[-1]) +'% error')
        ax.set_ylabel('error [%]')
        ax.set_xlabel('iterations')
        ax.grid(True)
        ax.legend(loc='best')
        
        plt.tight_layout()
        plt.show()

    #     plt.close()

    Time=(ds*10)
    dirac_sum=dirac_sum[:Time*fs]

    return dirac_sum,ds

def IteraTraceRF(traceZ,traceR,iterations=30,flag_stages=False,flag_summary=False):

    RF=traceR.copy()

    RF.data,ds = IterativeRF(
        traceZ,
        traceR,
        iterations=iterations,
        flag_stages=flag_stages,
        flag_summary=flag_summary)

    RF.stats.channel = 'RRF'
    RF.stats.npts = len(traceR.data)

    return RF,ds

################################
# 2 - Original traces rotation #
################################

def RoTrace(E,N,Z,baz):

    # Applies a baz angle rotation
    # Input: E,N,Z,baz
    # Output: T,R,Z
    # E, N, Z = 1-D arrays
    # baz = back-azimuth angle from station to earthquake
    # Takes E, N, Z arrays and rotates them clockwise using the back-azimuth angle
    # E -> T
    # N -> R
    # After the clockwise rotation, both R and T are reversed
    # R reversal is necessary so that direct P-arrival always look positive on the RF

    # MS on 28.11.2021 I don't remember why I reversed T sign
    # For details check the 3D rotation matrix M 

    x=np.zeros((3,len(E)))
    x[0,:]=E
    x[1,:]=N
    x[2,:]=Z

    angle=np.deg2rad(baz)
    c = np.cos(angle)
    s = np.sin(angle)
    M = np.array([[-c,+s,0],[-s,-c,0],[0,0,1]]) 

    y=np.matmul(M,x)

    return y[0,:],y[1,:],y[2,:]

#########################
# 3 - Save Event traces #
#########################

def SaveEventTrace(trace,event,info,tbefore,pathout):

    # Saves Event trace to database for later RFs computation

    stla,stlo,stel,baz,az,dist,tbefore=info

    mag=event.magnitudes[0].mag
    date=event.origins[0].time
    network=trace.stats.network

    header = {
        'knetwk':trace.stats.network,
        'kcmpnm':trace.stats.channel,
        'kstnm':trace.stats.station,
        'stla':stla,'stlo':stlo,'stel':stel,
        'evla':event.origins[0].latitude,
        'evlo':event.origins[0].longitude,
        'evdp':event.origins[0].depth,
        'mag':mag,
        'a':tbefore,
        'az':az,
        'baz':baz,
        'dist':dist,
        'nzyear':date.year,
        'nzjday':date.julday,
        'nzhour':date.hour,
        'nzmin':date.minute,
        'nzsec':date.second,
        'delta':trace.stats.delta}

    filename=pathout\
        +str(date.year)+'.'\
        +str(date.julday)+'.'\
        +str(date.hour)+'.'\
        +str(date.minute)+'.'\
        +str(date.second)+'.'\
        +trace.stats.station+'.'\
        +trace.stats.channel+'.SAC'

    tracetofile=sac.SACTrace(data=trace.data,**header) 
    tracetofile.write(filename)

    return

###########################################
# 4 - Get Event from catalog and raw data #
###########################################

def FileName(station,network,path,date,channel='HH'):

    # Assuming seismic raw database is divided into daily recordings in the following format:
    # ../NETWORK/STATION/STATION.NETWORK.AREA.CHANNEL.YEAR.JULDAY
    # returns the filenames to be read to look for the earthquake event arrival

    area={
        'XK':'IVREA',
        'IV':'INGV'}

    filenameZ=path+network+'/'+station+'/'+station+'.'+network+'.'+area[network]+'.'+channel+'Z'+str(date.year)+'.'+str(date.julday)
    filenameE=path+network+'/'+station+'/'+station+'.'+network+'.'+area[network]+'.'+channel+'E'+str(date.year)+'.'+str(date.julday)
    filenameN=path+network+'/'+station+'/'+station+'.'+network+'.'+area[network]+'.'+channel+'N'+str(date.year)+'.'+str(date.julday)

    return filenameZ,filenameE,filenameN

def Pick(event,station,network,path,latdict,londict,tbefore,tafter):

    model = TauPyModel(model='iasp91') # Select model for Direct-P travel time computation

    stla=latdict[station]
    stlo=londict[station]

    print(stlo,stla)

    # Seismic event location and occurance time

    evla=event.origins[0].latitude
    evlo=event.origins[0].longitude
    evdp=event.origins[0].depth/1000 # [km]
    date=event.origins[0].time

    print(evlo,evla)
    print(evdp,date)

    dist,az,baz=gps2dist(evla,evlo,stla,stlo)

    print('Distance in degrees: ',km2deg(dist/1000))
    print('az: ',az)
    print('baz: ',baz)

    # P-wave traveltime

    Ptraveltime = model.get_travel_times(
        source_depth_in_km=evdp,
        distance_in_degree=km2deg(dist/1000),
        phase_list=["P"])[0].time

    Parrivaltime=date+Ptraveltime

    # Check for recordings

    rawZ,rawE,rawN=FileName(station,network,path,date)

    print('rawZ ',rawZ)
    print('rawE ',rawE)
    print('rawN ',rawN)

    # Read raw stream components 

    if os.path.isfile(rawZ) and os.path.isfile(rawN) and os.path.isfile(rawE): 

        rawZ=obspy.read(rawZ)
        rawE=obspy.read(rawE)
        rawN=obspy.read(rawN)

    else:
        # Missing component
        print('Event not available')
        return False,-9999,-9999,-9999

    # Check if earthquake arrival is included in the data

    if (Parrivaltime+tafter).julday>date.julday:

        print('Event travelling during midnight')

        # Earthquake signal travelling across midnight
        # Following day required

        extraday=obspy.UTCDateTime(
            year=date.year,
            julday=date.julday+1)

        extrarawZ,extrarawE,extrarawN=FileName(station,network,path,extraday)

        rawZ.append(extrarawZ).merge()
        rawE.append(extrarawE).merge()
        rawN.append(extrarawN).merge()

        # Convert to obspy traces

        rawZ=rawZ[0]
        rawE=rawE[0]
        rawN=rawN[0]

    else:

        print('Daily trace available')

        # Convert to obspy traces

        rawZ=rawZ[0]
        rawE=rawE[0]
        rawN=rawN[0]

    # Check raw Z,E,N files contain the target time window

    if rawZ.stats.starttime>Parrivaltime-tbefore or rawZ.stats.endtime<Parrivaltime+tafter:
        print('rawZ.stats.starttime ',rawZ.stats.starttime)
        print('rawZ.stats.endtime ',rawZ.stats.endtime)
        print('Parrivaltime ',Parrivaltime)
        print('tbefore ',tbefore)
        print('tafter ',tafter)
        print('Not enough data')
        # Not enough data
        return False,-9999,-9999,-9999
    elif rawE.stats.starttime>Parrivaltime-tbefore or rawE.stats.endtime<Parrivaltime+tafter:
        print('rawE.stats.starttime ',rawE.stats.starttime)
        print('rawE.stats.endtime ',rawE.stats.endtime)
        print('Parrivaltime ',Parrivaltime)
        print('tbefore ',tbefore)
        print('tafter ',tafter)
        print('Not enough data')
        # Not enough data
        return False,-9999,-9999,-9999
    elif rawN.stats.starttime>Parrivaltime-tbefore or rawN.stats.endtime<Parrivaltime+tafter:
        print('rawN.stats.starttime ',rawN.stats.starttime)
        print('rawN.stats.endtime ',rawN.stats.endtime)
        print('Parrivaltime ',Parrivaltime)
        print('tbefore ',tbefore)
        print('tafter ',tafter)
        print('Not enough data')
        print('Not enough data')
        # Not enough data
        return False,-9999,-9999,-9999

    # Cut data around the window of interest

    threshold=5 # [s]

    dt=tafter+tbefore
    t0=Parrivaltime-tbefore
    t1=Parrivaltime+tafter

    rawZ.trim(
        starttime=t0,
        endtime=t1,
        nearest_sample=False)
    rawE.trim(
        starttime=t0,
        endtime=t1,
        nearest_sample=False)
    rawN.trim(
        starttime=t0,
        endtime=t1,
        nearest_sample=False)

    # Check there is no hole in the window of interest

    if np.abs(rawZ.stats.npts*rawZ.stats.delta-dt)>threshold:
        # More than 5s gap
        return False
    if np.abs(rawN.stats.npts*rawN.stats.delta-dt)>threshold:
        # More than 5s gap
        return False
    if np.abs(rawE.stats.npts*rawE.stats.delta-dt)>threshold:
        # More than 5s gap
        return False

    # Store precious information

    rawZ.stats.station=station
    rawE.stats.station=station
    rawN.stats.station=station

    rawZ.stats.network=network
    rawE.stats.network=network
    rawN.stats.network=network

    # Merge into one stream

    stream=obspy.Stream()
    stream.append(rawZ)
    stream.append(rawE)
    stream.append(rawN)

    stream.plot()

    print(event.__str__())

    return stream,baz,az,dist

def RemoveResponse(stream):

    # Pre-filtering

    newstream=obspy.Stream()

    freqmax=20
    freqmin=0.005

    for trace in stream:

        trace.detrend('demean')
        trace.detrend('linear')
        trace.taper(
            type='cosine',
            max_percentage=0.1,
            max_length=300)
        trace.filter(
            'bandpass',
            freqmin=freqmin,
            freqmax=freqmax,
            corners=4,
            zerophase=True)
        trace.taper(
            type='cosine',
            max_percentage=0.1,
            max_length=300)

        # Remove response

        inv=client.get_stations(
            network=trace.stats.network,
            station=trace.stats.station,
            channel=trace.stats.channel,
            level='response',
            starttime=trace.stats.starttime,
            endtime=trace.stats.endtime)
        try:
            trace.remove_response(
                inventory=inv,
                output='VEL',
                water_level=60,
                pre_filt=None,
                zero_mean=True,
                taper=True,
                taper_fraction=0.05,
                plot=False)
        except ValueError:
            return False

        newstream.append(trace)

    return newstream