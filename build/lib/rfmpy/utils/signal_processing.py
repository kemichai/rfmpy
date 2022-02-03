"""
Functions for signal processing.

Original codes by Matteo Scarponi on 30.11.2021
Location: Chavannes-pres-renens, CH
Date: Jan 2022
Author: Konstantinos Michailos
"""
import numpy as np
import obspy


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
    time = np.arange(-sigma * 5, sigma * 5 + delta, delta)
    gauss = np.exp(-time * time/(2 * sigma * sigma))
    spike_trace.data = signal.convolve(spike_trace.data, gauss, mode='same')

    return spike_trace


def correct_orientations(east, north, vertical):
    """

    """
    from obspy.signal.rotate import rotate2zne
    from obspy import Stream
    from obspy import read_inventory

    inv = read_inventory('/home/kmichall/Desktop/Codes/github/rfmpy/rfmpy/metadata/*.xml')

    v_corr = []
    e_corr = []
    n_corr = []
    for i, trace in enumerate(vertical):
        # print(z_trace, north_comp_traces[i])
        trace_n = north[i]
        trace_e = east[i]
        trace_z = vertical[i]

        orig_stream = Stream()
        orig_stream.append(trace_e)
        orig_stream.append(trace_n)
        orig_stream.append(trace_z)

        # orig_stream.plot()
        ####################
        def get_trace_name(tr):
            return tr.stats.network + '.' + tr.stats.station + '.' + tr.stats.location + '.' + tr.stats.channel

        e_trace_name = get_trace_name(trace_e)
        n_trace_name = get_trace_name(trace_n)
        z_trace_name = get_trace_name(trace_z)
        # Time of trace to choose the right epoch
        trace_time = trace.stats.starttime
        for net in inv:
            for sta in net:
                for cha in sta:
                    cha_name = net.code + '.' + sta.code + '.' + cha.location_code + '.' + cha.code
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
        tr_z.data, tr_n.data, tr_e.data = rotate2zne(trace_z.data, z_trace_az, z_trace_dip,
                                                     trace_n.data, n_trace_az, n_trace_dip,
                                                     trace_e.data, e_trace_az, e_trace_dip,
                                                     inverse=False)
        rot_stream = Stream()
        rot_stream.append(tr_e)
        rot_stream.append(tr_n)
        rot_stream.append(tr_z)
        all = Stream()
        all = orig_stream + rot_stream
        if n_trace_az != 0.0 or e_trace_az != 90.0:
            print(n_trace_az, e_trace_az)
            all.plot()
        # Change channel names for borehole sites
        if tr_n.stats.channel[-1] != 'N':
            print(tr_n.stats.channel)
            tr_n.stats.channel = tr_n.stats.channel[0:2] + 'N'
            tr_n.stats.sac.kcmpnm = tr_n.stats.sac.kcmpnm[0:2] + 'N'
        if tr_e.stats.channel[-1] != 'E':
            print(tr_e.stats.channel)
            tr_e.stats.channel = tr_n.stats.channel[0:2] + 'E'
            tr_e.stats.sac.kcmpnm = tr_n.stats.sac.kcmpnm[0:2] + 'E'

        v_corr.append(tr_z)
        e_corr.append(tr_e)
        n_corr.append(tr_n)

    return e_corr, n_corr, v_corr
