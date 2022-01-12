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
