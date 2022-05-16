"""
Function to create waveforms in the format for calculating RFs

Note: Functions here call the ones from utils.

Original codes by Matteo Scarponi on 30.11.2021
Location: Chavannes-pres-renens, CH
Date: Jan 2022
Author: Konstantinos Michailos
"""

import rfmpy.utils.RF_Util as rf_util
from rfmpy.utils.signal_processing import rotate_trace, remove_response, ConvGauss
import itertools
from pathlib import Path


def prep_wavs4rf(catalog, inventory, wav_directory, output_dir, time_before=60, time_after=120,
                 filt_order=2, highcut=2, low_cut=0.1, samp_rate=20, downsampling=False):
    """
    Function that cuts waveforms around the P wave arrival times of teleseismic events of interest.

    :type catalog: obspy.core.event.Catalog
    :param catalog: Obspy QUAKEML catalog containing info about the teleseismic events.
    :type inventory: obspy.core.inventory.inventory.Inventory
    :param inventory: Obspy inventory class containing station metadata.
    :type wav_directory: str
    :param wav_directory: Directory where the continuous waveform data are stored.
    :type output_dir: str
    :param output_dir: Directory where the cut waveform data will be stored.
    :type time_after: float
    :param time_after: Time after the P-wave arrival to end the cut window (default is 120 s)
    :type time_before: float
    :param time_before: Time before the P-wave arrival to start the cut window (default is 60 s).
    :type filt_order: int
    :param filt_order: Number of corners for filter.
    :type lowcut: float
    :param lowcut: Low cut for bandpass in Hz.
    :type highcut: float
    :param highcut: High cut for bandpass in Hz.
    :type samp_rate: float
    :param samp_rate: Sampling rate desired in Hz.
    :type downsampling: bool
    :param downsampling: Whether to downsample the data or not.

    :returns: SAC waveform files ready for RF calculations
    :rtype: obspy.core.stream.Trace

    .. note::

    """
    if highcut and highcut >= 0.5 * samp_rate:
        raise IOError('Highcut must be lower than the nyquist')

    stations = rf_util.get_station_details(inventory)
    for event, seedid in itertools.product(catalog, stations):
        # Event details
        origin_time = (event.preferred_origin() or event.origins[0])['time']
        # Station details
        net, sta, loc, cha = seedid.split('.')
        try:
            args = (seedid[:-1] + stations[seedid], origin_time)
            coords = inventory.get_coordinates(*args)
        except Exception as e:
            print(e)
            print(">>> Error while trying to use ", args[0])
        stream, baz, az, dist = rf_util.cut_waveforms(event=event, station=sta, network=net, path=wav_directory,
                                                      sta_latitude=coords['latitude'],
                                                      sta_longitude=coords['longitude'],
                                                      t_before=time_before, t_after=time_after,
                                                      channel_first_letters=cha[0:2])
        print(">>> Read data from station", sta, "for earthquake:", origin_time)
        if len(stream) != 3:
            from warnings import warn
            warn('Need 3 component seismograms. %d components '
                 'detected for event %s, station %s.' % (len(stream), origin_time, seedid))
            continue

        # Remove instrument response
        stream = remove_response(stream, high_cut=20, low_cut=0.005, filt_order=4)
        if not stream:
            continue
        # Bandpass filter and resample
        stream.detrend('demean')
        stream.detrend('linear')
        stream.taper(max_percentage=0.1, type='cosine', max_length=5)
        stream.filter('bandpass', freqmin=low_cut, freqmax=highcut, corners=filt_order, zerophase=True)
        stream.taper(max_percentage=0.1, type='cosine', max_length=5)

        if downsampling:
            stream.resample(samp_rate, window='hanning', no_filter=True)

        # Save Event seismic recording E N Z components
        stla = coords['latitude']
        stlo = coords['longitude']
        stel = coords['elevation']

        # Store information to a variable to pass it to the function
        info = stla, stlo, stel, baz, az, dist, time_before

        Path(output_dir).mkdir(exist_ok=True)

        rf_util.save_event_traces(trace=stream[0], event=event, info=info, path_out=output_dir)
        rf_util.save_event_traces(trace=stream[1], event=event, info=info, path_out=output_dir)
        rf_util.save_event_traces(trace=stream[2], event=event, info=info, path_out=output_dir)
    return
