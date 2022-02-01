"""
Bunch of functions...

Original codes by Matteo Scarponi on 30.11.2021
Location: Chavannes-pres-renens, CH
Date: Jan 2022
Author: Konstantinos Michailos
"""

from obspy.taup import TauPyModel
from obspy.geodetics.base import gps2dist_azimuth as gps2dist
from obspy.geodetics.base import kilometers2degrees as km2deg
import os.path
import obspy
import numpy as np
import obspy.io.sac.sactrace as sac
import matplotlib.pyplot as plt
from obspy import read
import matplotlib
from pathlib import Path
from scipy import signal
import glob


def store_receiver_functions(trace, path_to_store_rf):
    """
    Function for storing calculated receiver functions in SAC files.

    NOTE from GH:
    -earthquake date, time, latitude, longitude, depth, magnitude
    -data: time series of amplitude and sampling rate
    -back-azimuth, epicentral distance, P wave ray parameter
    -location of P arrival in the RF time series

    :type trace: obspy.core.trace.Trace
    :param trace: Waveform trace to store receiver function (either RRF or TRF).
    :type path_to_store_rf: str
    :param path_to_store_rf: Directory where traces will be stored.

    :returns: Final receiver function waveforms.
    """

    # Station info
    station_lat = trace.stats.sac.stla
    station_lon = trace.stats.sac.stlo
    station_ele = trace.stats.sac.stel
    station_name = trace.stats.station
    # Earthquake info
    event_lat = trace.stats.sac.evla
    event_lon = trace.stats.sac.evlo
    event_dep = trace.stats.sac.evdp
    event_mag = trace.stats.sac.mag
    dist_, az, baz = gps2dist(event_lat, event_lon, station_lat, station_lon)
    dist = km2deg(dist_/1000)
    y = trace.stats.sac.nzyear
    d = trace.stats.sac.nzjday
    h = trace.stats.sac.nzhour
    m = trace.stats.sac.nzmin
    s = trace.stats.sac.nzsec

    header = {'knetwk' : trace.stats.network, 'kcmpnm': trace.stats.channel,
              'kstnm': trace.stats.station, 'stla': station_lat, 'stlo': station_lon, 'stel': station_ele,
              'evla': event_lat, 'evlo': event_lon, 'evdp': event_dep, 'mag' : event_mag,
              'az': az, 'baz': baz, 'dist': dist, 'nzyear': y, 'a': trace.stats.sac.a,
              'nzjday': d, 'nzhour': h, 'nzmin': m, 'nzsec': s,
              'delta': trace.stats.sac.delta}

    julian_day = str(d)
    ev_h = str(h)
    ev_m = str(m)
    ev_s = str(s)
    if len(julian_day) == 1:
        julian_day = '00' + julian_day
    elif len(julian_day) == 2:
        julian_day = '0' + julian_day
    if len(ev_h) == 1:
        ev_h = '0' + ev_h
    if len(ev_m) == 1:
        ev_m = '0' + ev_m
    if len(ev_s) == 1:
        ev_s = '0' + ev_s

    rf_filename = (str(y) + '.' + julian_day + '.' + ev_h + '.' + ev_m + '.' + ev_s + '.' +
                   trace.stats.network + '.' + str(station_name) + '.' + trace.stats.channel + '.SAC')
    RF_to_file = sac.SACTrace(data=trace.data, **header)
    RF_to_file.write(path_to_store_rf + rf_filename)

    return


def IterativeRF(trace_z, trace_r, iterations=100, ds=30, iteration_plots=False, summary_plot=False):
    """
    Implementation of the Iterative deconvolution method.

    Reference: Ligorria, J. P., & Ammon, C. J. (1999). Iterative Deconvolution
    and Receiver-Function Estimation. Bulletin of the Seismological Society of
    America, 89, 5.

    :type trace_z: obspy.core.stream.Stream
    :param trace_z: Vertical component traces.
    :type trace_r: obspy.core.stream.Stream
    :param trace_r: Radial component traces.
    :type iterations: int
    :param iterations: Number of iterations for the deconvolution (default is 120).
    :type ds: int
    :param ds: Seconds from zero to align P-wave arrival (default is 30 seconds).
    :type iteration_plots: bool
    :param iteration_plots: Plot each iteration's plot.
    :type summary_plot: bool
    :param summary_plot: Plot a summary plot in the end.


    :returns:

    """

    fs = int(trace_r.stats.sampling_rate)
    trZ = trace_z.data
    trR = trace_r.data
    # Cutting first ds seconds from the Z trace to create delay between R and Z traces
    # This means that the direct-P arrival (or main reference peak)
    # in the RFs should be at exactly "ds" second from zero.
    # Cut ds seconds from Z so that direct P-arrival appears at t==ds in the final RF trace
    # NOTE we use 30 for this work... which should be the t_before

    delay = ds * fs
    trZ = trZ[delay:]
    nz = len(trZ)
    nr = len(trR)

    if nz > nr:
        trZ = trZ[:nr]
    # Prepare empty trace where to store cross-correlation maxima iteratively
    dirac_sum = np.zeros(nr)
    xcz = signal.correlate(trZ, trZ, mode='full')
    mxcz = np.max(np.abs(xcz))
    trH = trR
    rms = []

    k = (2*nr-1)//2
    iteration = 0
    while iteration < iterations:
        iteration += 1
        dirac = np.zeros(nr)
        xcr = signal.correlate(trZ, trH, mode='full', method='fft')
        ixcr = np.argmax(abs(xcr[nr-nz:nr]))
        shift = nz - ixcr - 1
        dirac[shift] = xcr[nr-nz+ixcr]/mxcz
             
        newH = signal.convolve(trZ, dirac)
        newH = newH[:nr]
        dirac_sum = dirac_sum + dirac

        conv = signal.convolve(trZ, dirac_sum)
        conv = conv[:len(trR)]
        diff = trR[:] - conv[:]
        diff = np.linalg.norm(diff)
        normConv = np.linalg.norm(conv)
        normR = np.linalg.norm(trR)
        diff = diff/(np.sqrt(normR*normConv))*100
        rms.append(diff)
        
        if iteration_plots:
            f = plt.figure(1)
            ax = plt.subplot(311)
            ttH = np.arange(len(trH))/fs
            ttZ = np.arange(len(trZ))/fs+ds
            ax.vlines(shift/fs, ymin=0, ymax=np.max(trZ), color='g', linestyle='--', label='peak')
            ax.plot(ttZ, trZ, 'k', label='trZ')
            ax.plot(ttH, trH, 'r', label='trH')
            ax.set_ylabel('Amplitude')
            ax.set_xlabel('Time [s]')
            ax.set_title('Iteration: ' + str(iteration) + '\nshift: ' + str(shift))
            ax.legend(loc='best')
            ax.grid(True)

            ax = plt.subplot(312)
            tt = np.arange(len(xcr))/fs
            ax.plot(tt, xcr, 'k', label='xcorr(Z,H)')
            ax.vlines(x=(nr-nz+ixcr)/fs, ymin=0, ymax=xcr[nr-nz+ixcr], color='r', label='max')
            ax.vlines(x=k/fs, ymin=0, ymax=xcr[nr-nz+ixcr], color='g', label='0 shift point')
            ax.set_ylabel('Amplitude')
            ax.set_xlabel('Time [s]')
            ax.set_title('Shift: ' + str(shift))
            ax.legend(loc='best')
            ax.grid(True)

            ax = plt.subplot(313)
            tt = np.arange(len(dirac_sum))/fs
            ax.plot(tt, dirac_sum, 'k', label='dirac_sum')
            ax.plot(tt, dirac, 'r', label='new spike')
            ax.set_ylabel('Amplitude')
            ax.set_xlabel('Time [s]')
            ax.set_title('Updating RF')
            ax.legend(loc='best')
            ax.grid(True)
            plt.tight_layout()
            plt.show()
        
        trH = trH - newH

    if summary_plot:
        f = plt.figure(2)
        ax = plt.subplot(311)
        tt = np.arange(len(dirac_sum))/fs
        ax.plot(tt, dirac_sum, 'k', lw=0.5, label='Computed RF')
        ax.fill_between(tt, y1=dirac_sum, y2=0, where=dirac_sum > 0, color='r')
        ax.fill_between(tt, y1=dirac_sum, y2=0, where=dirac_sum < 0, color='b')
        ax.set_title('Radial receiver function ' + trace_z.stats.station)
        ax.set_ylabel('Amplitude')
        ax.set_xlabel('Time (s)')
        ax.grid(True, alpha=0.5, lw=0.2)
        ax.legend(loc='best')
        
        ax = plt.subplot(312)
        tt = np.arange(len(trR))/fs
        ax.plot(tt, trR, 'k', label='R component')
        ax.plot(tt, conv, 'r', label='R approx')
        ax.set_title(str(iteration) + 'th iteration')
        ax.set_ylabel('Amplitude [counts]')
        ax.set_xlabel('Time (s)')
        ax.legend(loc='best')
        ax.grid(True)
        
        ax = plt.subplot(313)
        ax.plot(range(len(rms)), rms, 'ro-', label='rms%')
        ax.set_title('Final rms: ' + str("%.2f" % rms[-1]) + '% error')
        ax.set_ylabel('error [%]')
        ax.set_xlabel('iterations')
        ax.grid(True)
        ax.legend(loc='best')
        
        plt.tight_layout()
        plt.show()

    time = (ds*10)
    dirac_sum = dirac_sum[:time*fs]

    return dirac_sum, ds


def IteraTraceRF(traceZ, traceR, iterations=100, flag_stages=False, flag_summary=False):
    # THIS IS NOT USED NEED TO CHECK WITH MATTEO...
    RF=traceR.copy()
    RF.data,ds = IterativeRF(traceZ, traceR, iterations=iterations, flag_stages=flag_stages, flag_summary=flag_summary)
    RF.stats.channel = 'RRF'
    RF.stats.npts = len(traceR.data)

    return RF, ds


def save_event_traces(trace, event, info, path_out):
    """
    Saves event traces to database for later RFs computation

    :type trace: obspy.core.trace.Trace
    :param trace: Trace to store as a SAC file.
    :type event: obspy.core.event.event.Event
    :param event: Teleseismic event object.
    :type info: list
    :param info: List containing various information about the teleseismic event and the seismic station.
    :type path_out: str
    :param path_out: Path to the directory where the files will be stored.

    :returns: A SAC file that will be later used for RF calculations.
    """

    # Fetch information about the event and the station
    stla, stlo, stel, baz, az, dist, time_before = info
    magnitude = event.magnitudes[0].mag
    date = event.origins[0].time
    origin_time = (event.preferred_origin() or event.origins[0])['time']

    header = {'knetwk': trace.stats.network, 'kcmpnm': trace.stats.channel, 'kstnm': trace.stats.station,
              'stla': stla, 'stlo': stlo, 'stel': stel,
              'evla': event.origins[0].latitude, 'evlo': event.origins[0].longitude, 'evdp': event.origins[0].depth,
              'mag': magnitude, 'a': time_before, 'az': az, 'baz': baz, 'dist': dist,
              'nzyear': date.year, 'nzjday': date.julday, 'nzhour': date.hour, 'nzmin': date.minute,
              'nzsec': date.second, 'delta': trace.stats.delta}

    # Making sure that the length of the file names are the same
    julian_day = str(origin_time.julday)
    ev_hour = str(origin_time.hour)
    ev_min = str(origin_time.minute)
    ev_sec = str(origin_time.second)

    if len(julian_day) == 1:
        julian_day = '00' + julian_day
    elif len(julian_day) == 2:
        julian_day = '0' + julian_day
    if len(ev_hour) == 1:
        ev_hour = '0' + ev_hour
    if len(ev_min) == 1:
        ev_min = '0' + ev_min
    if len(ev_sec) == 1:
        ev_sec = '0' + ev_sec

    # P_year.julian_day.hour.minute.seconds
    folder_name = (path_out + 'P_' + str(origin_time.year) + '.' + julian_day + '.'
                   + str(ev_hour) + '.' + str(ev_min) + '.'
                   + str(ev_sec) + '/')
    Path(folder_name).mkdir(exist_ok=True)
    # Year.julian_day.hour.minute.seconds
    filename = (folder_name + str(origin_time.year) + '.' + julian_day + '.' + ev_hour + '.' + ev_min +
                '.' + ev_sec + '.' + trace.stats.network + '.' + trace.stats.station + '.' + trace.stats.channel
                + '.SAC')

    trace2file = sac.SACTrace(data=trace.data, **header)
    trace2file.write(filename)

    return


def define_filenames(station, network, path, date, channel):
    """
    Defines the paths to the continuous waveform data to read the data for a given station at a given time.

    :type station: str
    :param station: Name of seismic site.
    :type network: str
    :param network: Name of the seismic network.
    :type path: str
    :param path: Path to attach to the filename.
    :type date: obspy.core.utcdatetime.UTCDateTime
    :param date: Origin time of the teleseismic earthquake.
    :type channel: str
    :param channel: Name of the specific channel.

    :returns: String with the path and file name.
    """

    julian_day = str(date.julday)
    if len(julian_day) == 1:
        julian_day = '00' + julian_day
    elif len(julian_day) == 2:
        julian_day = '0' + julian_day

    filename = (path + network + '/Y' + str(date.year) + '/R' + julian_day + '.01' +
                '/' + network + '.' + station + '..' + channel + '.' + str(date.year) + '.' + julian_day)
    return filename


def cut_waveforms(event, station, network, path, sta_latitude, sta_longitude, t_before, t_after, channel_first_letters):
    """
    Cuts waveforms around the P-wave arrival of a teleseismic event.

    :type event: obspy.core.event.event.Event
    :param event: Teleseismic event.
    :type station: str
    :param station: Station name.
    :type network: str
    :param network: Network name.
    :type path: str
    :param path: Directory path.
    :type sta_latitude: float
    :param sta_latitude: Station's latitude.
    :type sta_longitude: float
    :param sta_longitude: Station's longitude.
    :type t_before: float
    :param t_before: Time before the P-wave arrival to start the cut window.
    :type t_after: float
    :param t_after: Time after the P-wave arrival to end the cut window.
    :type channel_first_letters: str
    :param channel_first_letters: Two first letters of the channel.

    :returns: The stream cut around the P arrival along with some info (i.e., baz, az, dist)
    """

    # Select model for Direct-P travel time computation
    model = TauPyModel(model='iasp91')

    # Tele-seismic event location and occurrence time
    evla = event.origins[0].latitude
    evlo = event.origins[0].longitude
    evdp = float(event.origins[0].depth)/1000  # [km]
    date = event.origins[0].time
    print('>>> Origin time:', date)
    print('>>> Longitude:', evlo)
    print('>>> Latitude:', evla)
    print('>>> Depth:', evdp)

    # Calculate distance azimuth and back-azimuth
    dist, az, baz = gps2dist(evla, evlo, sta_latitude, sta_longitude)
    print('>>> Distance in degrees: ', round(km2deg(dist/1000), 2))
    print('>>> Azimuth            : ', round(az, 2))
    print('>>> Backazimuth        : ', round(baz, 2))

    # Calculate the P-wave travel time
    p_travel_time = model.get_travel_times(source_depth_in_km=evdp, distance_in_degree=km2deg(dist/1000),
                                           phase_list=["P"])[0].time
    p_arrival_time = date + p_travel_time

    # Check for recordings
    wav_file_Z = define_filenames(station, network, path, date, channel=channel_first_letters + 'Z')
    wav_file_E = define_filenames(station, network, path, date, channel=channel_first_letters + 'E')
    wav_file_N = define_filenames(station, network, path, date, channel=channel_first_letters + 'N')
    print(wav_file_Z)
    # Read raw stream components
    if os.path.isfile(wav_file_Z) and os.path.isfile(wav_file_N) and os.path.isfile(wav_file_E):
        rawZ = read(wav_file_Z)
        rawE = read(wav_file_E)
        rawN = read(wav_file_N)
    else:
        # Missing component
        print('One of the components in the waveform data missing')
        return [], -9999, -9999, -9999

    # Check if earthquake arrival is included in the data

    if (p_arrival_time + t_after).julday > date.julday:
        print('Event travelling during midnight')
        # Earthquake signal travelling across midnight
        # Following day required
        extraday = obspy.UTCDateTime(year=date.year, julday=date.julday+1)
        extrarawZ = define_filenames(station, network, path, extraday, channel=channel_first_letters + 'Z')
        extrarawE = define_filenames(station, network, path, extraday, channel=channel_first_letters + 'E')
        extrarawN = define_filenames(station, network, path, extraday, channel=channel_first_letters + 'N')
        rawZ.append(extrarawZ).merge()
        rawE.append(extrarawE).merge()
        rawN.append(extrarawN).merge()

        # Convert to obspy traces
        rawZ = rawZ[0]
        rawE = rawE[0]
        rawN = rawN[0]
    else:
        print('>>> Daily trace available!')
        # Convert to obspy traces
        rawZ = rawZ[0]
        rawE = rawE[0]
        rawN = rawN[0]
    # Check raw Z,E,N files contain the target time window
    if rawZ.stats.starttime > p_arrival_time - t_before or rawZ.stats.endtime < p_arrival_time + t_after:
        print('NOT ENOUGH DATA')
        print('rawZ.stats.starttime ', rawZ.stats.starttime)
        print('rawZ.stats.endtime ', rawZ.stats.endtime)
        print('Parrivaltime ', p_arrival_time)
        print('tbefore ', t_before)
        print('tafter ', t_after)
        return [], -9999, -9999, -9999
    elif rawE.stats.starttime > p_arrival_time - t_before or rawE.stats.endtime < p_arrival_time + t_after:
        print('NOT ENOUGH DATA')
        print('rawE.stats.starttime ', rawE.stats.starttime)
        print('rawE.stats.endtime ', rawE.stats.endtime)
        print('Parrivaltime ', p_arrival_time)
        print('tbefore ', t_before)
        print('tafter ', t_after)
        return [], -9999, -9999, -9999
    elif rawN.stats.starttime > p_arrival_time - t_before or rawN.stats.endtime < p_arrival_time + t_after:
        print('NOT ENOUGH DATA')
        print('rawN.stats.starttime ', rawN.stats.starttime)
        print('rawN.stats.endtime ', rawN.stats.endtime)
        print('Parrivaltime ', p_arrival_time)
        print('tbefore ', t_before)
        print('tafter ', t_after)
        print('Not enough data')
        return [], -9999, -9999, -9999

    # Cut data around the window of interest in seconds
    threshold = 5
    dt = t_after + t_before
    t0 = p_arrival_time - t_before
    t1 = p_arrival_time + t_after

    rawZ.trim(starttime=t0, endtime=t1, nearest_sample=False)
    rawE.trim(starttime=t0, endtime=t1, nearest_sample=False)
    rawN.trim(starttime=t0, endtime=t1, nearest_sample=False)
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

    # Store information
    rawZ.stats.station = station
    rawE.stats.station = station
    rawN.stats.station = station
    rawZ.stats.network = network
    rawE.stats.network = network
    rawN.stats.network = network
    # Merge into one stream
    stream = obspy.Stream()
    stream.append(rawZ)
    stream.append(rawE)
    stream.append(rawN)

    return stream, baz, az, dist


def get_station_details(inventory):
    """Gives the details of the seismic sites from a stationxml file."""
    channels = inventory.get_contents()['channels']
    stations = {ch[:-1] + '?': ch[-1] for ch in channels}
    return stations


def get_unique_stations(event_dir):
    """
    Creates a list of all the stations for which we have waveform data.

    NOTE: Currently written so if the waveforms for a given function does not have
          a channel ending in E it will not add the other two channes (N, Z).
          This WILL NOT INCLUDE STATIONS WITH 1, 2 FOR THE HORIZONTAL COMPONENTS.
    """
    # TODO: sort out this function

    import glob
    import os.path
    wav_files = glob.glob(event_dir + '/*SAC')
    unique_station_list = []
    for wav_file in wav_files:
        station = wav_file.split('/')[-1].split('.')[-3]
        channel = wav_file.split('/')[-1].split('.')[-2]
        if station not in unique_station_list:
            unique_station_list.append(station)

    station_list = []
    for station_name in unique_station_list:
        for wav_file in wav_files:
            station_ = wav_file.split('/')[-1].split('.')[-3]
            channel_ = wav_file.split('/')[-1].split('.')[-2]
            if station_name == station_ and channel_[-1] == 'Z' and os.path.isfile(wav_file):
                N_comp = '.'.join(wav_file.split('.')[:11]) + '.' + channel_[0:2] + 'N.SAC'
                E_comp = '.'.join(wav_file.split('.')[:11]) + '.' + channel_[0:2] + 'E.SAC'
                if os.path.isfile(N_comp) and os.path.isfile(E_comp):
                    # Check that we have one stream for each component before we proceed
                    station_list.append(station_ + '.' + channel_)
                    N_channel = channel_[0:2] + 'N'
                    station_list.append(station_ + '.' + N_channel)
                    E_channel = channel_[0:2] + 'E'
                    station_list.append(station_ + '.' + E_channel)
                else:
                    continue
                    print('We do not have data from all three components for station', station_)
                # TODO: at the moment the stations with channels that are not N or E (e.g., 1,2) are not used...
    return station_list


def printing_station_name(station_name, station_network):
    """"""
    if len(station_name) == 3:
        station_name = station_name + '  '
    elif len(station_name) == 4:
        station_name = station_name + ' '

    station_name2print = station_network + '.' + station_name

    return station_name2print


def _dip_azimuth2ZSE_base_vector(dip, azimuth):
    """
    Helper function converting a vector described with azimuth and dip of unit
    length to a vector in the ZSE (Vertical, South, East) base.
    The definition of azimuth and dip is according to the SEED reference
    manual, as are the following examples (they use rounding for small
    numerical inaccuracies - also positive and negative zero are treated as
    equal):
    >>> r = lambda x: np.array([_i if _i != -0.0 else 0.0\
        for _i in np.round(x, 10)])
    >>> r(_dip_azimuth2ZSE_base_vector(-90, 0)) #doctest: +NORMALIZE_WHITESPACE
    array([ 1., 0., 0.])
    >>> r(_dip_azimuth2ZSE_base_vector(90, 0)) #doctest: +NORMALIZE_WHITESPACE
    array([-1., 0., 0.])
    >>> r(_dip_azimuth2ZSE_base_vector(0, 0)) #doctest: +NORMALIZE_WHITESPACE
    array([ 0., -1., 0.])
    >>> r(_dip_azimuth2ZSE_base_vector(0, 180)) #doctest: +NORMALIZE_WHITESPACE
    array([ 0., 1., 0.])
    >>> r(_dip_azimuth2ZSE_base_vector(0, 90)) #doctest: +NORMALIZE_WHITESPACE
    array([ 0., 0., 1.])
    >>> r(_dip_azimuth2ZSE_base_vector(0, 270)) #doctest: +NORMALIZE_WHITESPACE
    array([ 0., 0., -1.])
    """
    # Convert both to radian.
    dip = np.deg2rad(dip)
    azimuth = np.deg2rad(azimuth)

    # Define the rotation axis for the dip.
    c1 = 0.0
    c2 = 0.0
    c3 = -1.0
    # Now the dip rotation matrix.
    dip_rotation_matrix = np.cos(dip) * \
        np.matrix(((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0))) + \
        (1 - np.cos(dip)) * np.matrix(((c1 * c1, c1 * c2, c1 * c3),
                                       (c2 * c1, c2 * c2, c2 * c3),
                                       (c3 * c1, c3 * c2, c3 * c3))) + \
        np.sin(dip) * np.matrix(((0, -c3, c2), (c3, 0, -c1), (-c2, c1, 0)))
    # Do the same for the azimuth.
    c1 = -1.0
    c2 = 0.0
    c3 = 0.0
    azimuth_rotation_matrix = np.cos(azimuth) * \
        np.matrix(((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0))) + \
        (1 - np.cos(azimuth)) * np.matrix(((c1 * c1, c1 * c2, c1 * c3),
                                           (c2 * c1, c2 * c2, c2 * c3),
                                           (c3 * c1, c3 * c2, c3 * c3))) + \
        np.sin(azimuth) * np.matrix(((0, -c3, c2), (c3, 0, -c1), (-c2, c1, 0)))

    # Now simply rotate a north pointing unit vector with both matrixes.
    temp = np.dot(azimuth_rotation_matrix, [[0.0], [-1.0], [0.0]])
    return np.array(np.dot(dip_rotation_matrix, temp)).ravel()


def rotate2ZNE(data_1, azimuth_1, dip_1, data_2, azimuth_2, dip_2, data_3,
               azimuth_3, dip_3):
    """
    Rotates an arbitrarily oriented three-component vector to ZNE.
    Each components orientation is described with a azimuth and a dip. The
    azimuth is defined as the degrees from North, clockwise and the dip is the
    defined as the number of degrees, down from horizontal. Both definitions
    are according to the SEED standard.
    The three components need not be orthogonal to each other but the
    components have to be linearly independent. The function performs a full
    base change to orthogonal Vertical, North, and East orientations.
    :param data_1: Data component 1.
    :param azimuth_1: The azimuth of component 1.
    :param dip_1: The dip of component 1.
    :param data_2: Data component 2.
    :param azimuth_2: The azimuth of component 2.
    :param dip_2: The dip of component 2.
    :param data_3: Data component 3.
    :param azimuth_3: The azimuth of component 3.
    :param dip_3: The dip of component 3.
    :rtype: Tuple of three NumPy arrays.
    :returns: The three rotated components, oriented in Z, N, and E.
    >>> # An input of ZNE yields an output of ZNE
    >>> rotate2ZNE(np.arange(3), 0, -90, np.arange(3) * 2, 0, 0, \
            np.arange(3) * 3, 90, 0) # doctest: +NORMALIZE_WHITESPACE
    (array([ 0., 1., 2.]), array([ 0., 2., 4.]), array([ 0., 3., 6.]))
    >>> # An input of ZSE yields an output of ZNE
    >>> rotate2ZNE(np.arange(3), 0, -90, np.arange(3) * 2, 180, 0, \
            np.arange(3) * 3, 90, 0) # doctest: +NORMALIZE_WHITESPACE
    (array([ 0., 1., 2.]), array([ 0., -2., -4.]), array([ 0., 3., 6.]))
    >>> # Mixed up components should get rotated to ZNE.
    >>> rotate2ZNE(np.arange(3), 0, 0, np.arange(3) * 2, 90, 0, \
            np.arange(3) * 3, 0, -90) # doctest: +NORMALIZE_WHITESPACE
    (array([ 0., 3., 6.]), array([ 0., 1., 2.]), array([ 0., 2., 4.]))
    """
    # Internally works in Vertical, South, and East components; a right handed
    # coordinate system.

    # Define the base vectors of the old base in terms of the new base vectors.
    base_vector_1 = _dip_azimuth2ZSE_base_vector(dip_1, azimuth_1)
    base_vector_2 = _dip_azimuth2ZSE_base_vector(dip_2, azimuth_2)
    base_vector_3 = _dip_azimuth2ZSE_base_vector(dip_3, azimuth_3)

    # Build transformation matrix.
    T = np.matrix([base_vector_1, base_vector_2, base_vector_3]).transpose()

    # Apply it.
    z, s, e = np.dot(T, [data_1, data_2, data_3])
    # Replace all negative zeros. These might confuse some further processing
    # programs.
    z = np.array(z).ravel()
    z[z == -0.0] = 0
    # Return a North pointing array.
    n = -1.0 * np.array(s).ravel()
    n[n == -0.0] = 0
    e = np.array(e).ravel()
    e[e == -0.0] = 0

    return z, n, e


def get_station_info(path_wavs_list):
    """
    Reads all the available waveform files and creates a list
    of the unique seismic stations including latitude, longitude
    and station elevation.

    :type path_wavs_list: list
    :param path_wavs_list: List of directory paths

    :returns: A list of seismic site details.
    """
    stations = []
    for path in path_wavs_list:
        print(path)
        all_event_dir = glob.glob(path + '*')
        for event_dir in all_event_dir:
            wav_files = glob.glob(event_dir + '/*Z.SAC')
            for wav_file in wav_files:
                # print(wav_file)
                tr = obspy.read(wav_file)
                try:
                    lat = str(tr[0].stats.sac.stla)
                    lon = str(tr[0].stats.sac.stlo)
                    # ele = str(tr[0].stats.sac.stel)
                except Exception as e:
                    print(f"Could not read {wav_file} due to {e}")
                    continue
                station = wav_file.split('/')[-1].split('.')[-3]
                channel = wav_file.split('/')[-1].split('.')[-2]
                network = wav_file.split('/')[-1].split('.')[-4]
                # station_name = network + '.' + station + '.' + channel + ' ' + lat + ' ' + lon + ' ' + ele
                station_name = network + '.' + station + ' ' + lat + ' ' + lon # + ' ' + ele
                if station_name not in stations:
                    stations.append(station_name)

    return stations
