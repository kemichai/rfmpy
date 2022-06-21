"""
Functions for signal processing.

Original codes by Matteo Scarponi on 30.11.2021
Location: Chavannes-pres-renens, CH
Date: Jan 2022
Author: Konstantinos Michailos
"""
import numpy as np
import obspy
import logging
import sys

# Create a log file
# a_logger = logging.getLogger()
# a_logger.setLevel(logging.INFO)
# output_file_handler = logging.FileHandler("logfile_SWISS.txt")
# stdout_handler = logging.StreamHandler(sys.stdout)
# a_logger.addHandler(output_file_handler)
# a_logger.addHandler(stdout_handler)



def rotate_trace(east, north, vertical, baz):
    """
    Applies a back-azimuth angle rotation to the traces. After the clockwise rotation, both R and T are
    reversed R reversal is necessary so that direct P-arrival always look positive on the RF.

    Consists of a horizontal rotation around Z to align along the baz direction
    Plus sign reversal of the R component so that main peak looks positive

    :type east: 1D array
    :param east: east component
    :type north: 1D array
    :param north: north component
    :type vertical: 1D array
    :param vertical: vertical component
    :type baz: float
    :param baz: back azimuth angle

    :returns: T, R, Z components

    .. note:: MS on 28.11.2021 I don't remember why I reversed T sign
              For details check the 3D rotation matrix M
    """

    x = np.zeros((3, len(east)))
    x[0, :] = east
    x[1, :] = north
    x[2, :] = vertical

    angle = np.deg2rad(baz)
    M = np.array([[-np.cos(angle), +np.sin(angle), 0],
                  [-np.sin(angle), -np.cos(angle), 0],
                  [0,               0,             1]])
    y = np.matmul(M, x)

    return y[0, :], y[1, :], y[2, :]



def remove_response(stream, high_cut=20, low_cut=0.005, filt_order=4):
    """
    Remove response from waveforms.

    :type stream: obspy.core.stream.Stream
    :param stream: Stream of traces.
    :type filt_order: int
    :param filt_order: Number of corners for filter.
    :type low_cut: float
    :param low_cut: Low cut for bandpass in Hz.
    :type high_cut: float
    :param high_cut: High cut for bandpass in Hz.

    :returns: Detrended traces.
    :rtype: obspy.core.stream.Stream
    """

    from obspy.clients.fdsn import RoutingClient

    client = RoutingClient('iris-federator')

    new_stream = obspy.Stream()
    for trace in stream:
        trace.detrend('demean')
        trace.detrend('linear')
        trace.taper(type='cosine', max_percentage=0.1, max_length=300)
        trace.filter('bandpass', freqmin=low_cut, freqmax=high_cut, corners=filt_order, zerophase=True)
        trace.taper(type='cosine', max_percentage=0.1, max_length=300)
        # Remove response
        inv = client.get_stations(network=trace.stats.network, station=trace.stats.station, channel=trace.stats.channel,
                                  level='response', starttime=trace.stats.starttime, endtime=trace.stats.endtime)
        try:
            trace.remove_response(inventory=inv, output='VEL', water_level=60, pre_filt=None, zero_mean=True,
                                  taper=True, taper_fraction=0.05, plot=False)
        except ValueError:
            return False
        new_stream.append(trace)

    return new_stream


def ConvGauss(spike_trace, high_cut, delta):
    """
    Convolve the spike train with a gaussian filter whose width is equal to
    the maximum frequency content of the signal.

    :type spike_trace: -
    :param spike_trace:
    :type high_cut: float
    :param high_cut: High cut for bandpass in Hz.
    :type delta: float
    :param delta: Sample distance in seconds.

    :returns: Convolved spikes.

    """
    from scipy import signal

    sigma = 1./(2 * np.pi * high_cut)
    # As MS used it.
    # time = np.arange(-sigma * 5, sigma * 5 + delta, delta)
    # As it is found in the Matlab codes.
    time = np.arange(-sigma * 5, sigma * 5, delta)
    gauss = np.exp(-time * time/(2 * sigma * sigma))
    spike_trace.data = signal.convolve(spike_trace.data, gauss, mode='same')

    return spike_trace


def correct_orientations(st_east, st_north, st_vertical, inventory, logfile, comparison_plot=False):
    """
    Corrects misaligned horizontal components of the
    AlpArray seismic sites using station metadata (stationxml file)
    information available from the different seismic
    network operators that contributed to the AlpArray Seismic Network.

    This takes into consideration the epoch of the traces to
    assign the correct azimuth and dip for rotations.

    :type st_east: obspy.core.stream.Stream
    :param st_east: Horizontal component waveform traces (East-West).
    :type st_north: obspy.core.stream.Stream
    :param st_north: Horizontal component waveform traces (North-South).
    :type st_vertical: obspy.core.stream.Stream
    :param st_vertical: Vertical component waveform traces.
    :type inventory: obspy.core.inventory.Inventory
    :param inventory: Inventory containing response information for the stations in st.
    :type comparison_plot: boolean
    :param comparison_plot: If True it will plot the original and rotated traces on top
                            of each other (default value is False).

    :returns: Streams of traces for north, east and vertical components
              with corrected orientations (according to the azimuth and dip values
              of each channel stored in the inventory file.
    """
    from obspy.signal.rotate import rotate2zne
    from obspy import Stream
    from obspy.core import UTCDateTime

    v_corr = []
    e_corr = []
    n_corr = []
    for i, trace in enumerate(st_vertical):
        trace_n = st_north[i]
        trace_e = st_east[i]
        trace_z = st_vertical[i]

        def get_trace_name(tr):
            return tr.stats.network + '.' + tr.stats.station + '.' + tr.stats.location + '.' + tr.stats.channel

        e_trace_name = get_trace_name(trace_e)
        n_trace_name = get_trace_name(trace_n)
        z_trace_name = get_trace_name(trace_z)
        # Time of trace to choose the right epoch
        trace_time = trace.stats.starttime
        for net in inventory:
            for sta in net:
                for cha in sta:
                    if cha.end_date == None:
                        cha.end_date = UTCDateTime('2500-12-31T23:59:59.000000Z')

                    cha_name = net.code + '.' + sta.code + '.' + cha.location_code + '.' + cha.code
                    # both name and time match
                    if e_trace_name == cha_name and trace_time > cha.start_date and trace_time < cha.end_date:
                        e_trace_az = cha.azimuth
                        e_trace_dip = cha.dip

                    if n_trace_name == cha_name and trace_time > cha.start_date and trace_time < cha.end_date:
                        n_trace_az = cha.azimuth
                        n_trace_dip = cha.dip

                    if z_trace_name == cha_name and trace_time > cha.start_date and trace_time < cha.end_date:
                        z_trace_az = cha.azimuth
                        z_trace_dip = cha.dip

        tr_e = trace_e.copy()
        tr_n = trace_n.copy()
        tr_z = trace_z.copy()
        try:
            logfile.info("|-----------------------------------------------|")
            tr_z.data, tr_n.data, tr_e.data = rotate2zne(trace_z.data, z_trace_az, z_trace_dip,
                                                         trace_n.data, n_trace_az, n_trace_dip,
                                                         trace_e.data, e_trace_az, e_trace_dip,
                                                         inverse=False)
            logfile.info(f"| Rotation applied to N trace: {trace_n.stats.station}; AZ: {n_trace_az}; DIP: {n_trace_dip}")
            logfile.info(f"| Rotation applied to E trace: {trace_e.stats.station}; AZ: {e_trace_az}; DIP: {e_trace_dip}")
        except Exception as e:
            # print(f"|No information found for trace: {trace_z.stats.station}")
            logfile.info(f"|No information found for trace: {trace_z.stats.station}")
        logfile.info("|-----------------------------------------------|")
        rot_stream = Stream()
        rot_stream.append(tr_e)
        rot_stream.append(tr_n)
        rot_stream.append(tr_z)

        if comparison_plot:
            orig_stream = Stream()
            orig_stream.append(trace_e)
            orig_stream.append(trace_n)
            orig_stream.append(trace_z)
            all = Stream()
            all = orig_stream + rot_stream
            print(n_trace_az, e_trace_az)
            all.plot()
        # Change channel names for borehole sites
        if tr_n.stats.channel[-1] != 'N':
            tr_n.stats.channel = tr_n.stats.channel[0:2] + 'N'
            tr_n.stats.sac.kcmpnm = tr_n.stats.sac.kcmpnm[0:2] + 'N'
        if tr_e.stats.channel[-1] != 'E':
            tr_e.stats.channel = tr_n.stats.channel[0:2] + 'E'
            tr_e.stats.sac.kcmpnm = tr_n.stats.sac.kcmpnm[0:2] + 'E'
        # Change channel names to BH*
        if tr_n.stats.channel[0] != 'B' or tr_n.stats.channel[0] != 'H':
            tr_n.stats.channel = 'B' + tr_n.stats.channel[1:3]
            tr_n.stats.sac.kcmpnm = 'B' + tr_n.stats.sac.kcmpnm[1:3]
        if tr_e.stats.channel[0] != 'B' or tr_e.stats.channel[0] != 'H':
            tr_e.stats.channel = 'B' + tr_e.stats.channel[1:3]
            tr_e.stats.sac.kcmpnm = 'B' + tr_e.stats.sac.kcmpnm[1:3]
        if tr_z.stats.channel[0] != 'B' or tr_z.stats.channel[0] != 'H':
            tr_z.stats.channel = 'B' + tr_z.stats.channel[1:3]
            tr_z.stats.sac.kcmpnm = 'B' + tr_z.stats.sac.kcmpnm[1:3]
        # Append to the list
        v_corr.append(tr_z)
        e_corr.append(tr_e)
        n_corr.append(tr_n)

    return e_corr, n_corr, v_corr


def pre_processing(R, T, Z, low_cut, high_cut, order, t_bef, t_aft):
    """
    Process waveforms before calculating RFs.
    1) bandpass filter
    2) demean
    3) taper

    :type R: obspy.core.trace.Trace
    :param R: Waveform trace of radial component.
    :type T: obspy.core.trace.Trace
    :param T: Waveform trace of transverse component.
    :type Z: obspy.core.trace.Trace
    :param Z: Waveform trace of vertical component.
    :type high_cut: float
    :param high_cut: High cut for bandpass in Hz.
    :type order: int
    :param order: Filter order to use.
    :type t_bef: int
    :param t_bef: Time before the P-wave arrival to start the cut window.
    :type t_aft: int
    :param t_aft: Time after the P-wave arrival to start the cut window.

    :returns: Processed traces.
    """

    from obspy.signal import filter
    import obspy.io.sac.sactrace as sac

    def process_trace(tr):
        tr.detrend('demean')
        tr.taper(max_percentage=0.5, type='hann', max_length=15, side='both')
        tr.filter(type='bandpass', freqmin=low_cut, freqmax=high_cut, corners=order, zerophase=True)

        # Trim 40 sec before and after the P wave arrival (total 240 seconds long traces)
        # trace_length = int(tr.stats.npts / tr.stats.sampling_rate)
        # P wave arrival should be in the middle of the trace
        # t0 = tr.stats.starttime + trace_length/2 - time_window
        # t1 = tr.stats.starttime + trace_length/2 + time_window

        time_before = tr.stats.sac.a
        fs = tr.stats.sampling_rate
        t0 = int((time_before - t_bef) * fs)
        t1 = int((time_before + t_aft) * fs)

        tr_sliced = tr.copy()
        tr_sliced.data = tr_sliced.data[t0:t1]
        # Updates endtime too
        tr_sliced.stats.starttime = tr_sliced.stats.starttime + t0/fs

        # update SAC header
        tr_sliced.stats.sac.nzyear = tr_sliced.stats.starttime.year
        tr_sliced.stats.sac.nzjday = tr_sliced.stats.starttime.julday
        tr_sliced.stats.sac.nzhour = tr_sliced.stats.starttime.hour
        tr_sliced.stats.sac.nzmin = tr_sliced.stats.starttime.minute
        tr_sliced.stats.sac.nzsec = tr_sliced.stats.starttime.second
        tr_sliced.stats.sac.nzmsec = tr_sliced.stats.starttime.microsecond

        return tr_sliced

    R_ = process_trace(R)
    T_ = process_trace(T)
    Z_ = process_trace(Z)

    return R_, T_, Z_
