"""
Main functions 1) to create waveforms in the format for calculating RFs
               2) for calculating RFs.

Note: Functions here call the ones from utils.

Original codes by Matteo Scarponi on 30.11.2021
Location: Chavannes-pres-renens, CH
Date: Jan 2022
Author: Konstantinos Michailos
"""

import core.utils.RF_Util as rf_util
from utils.signal_processing import rotate_trace, remove_response, ConvGauss
from core.utils.qc import rms_quality_control, rf_quality_control
from obspy import read_inventory, read_events, UTCDateTime as UTC
import itertools
from pathlib import Path
import glob
import obspy
import numpy as np


def prep_wavs4rf(catalog, inventory, wav_directory, output_dir, time_before=60, time_after=120,
                             filt_order=2, highcut=2, lowcut=0.1, samp_rate=20, downsampling=False):
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
        stream.filter('bandpass', freqmin=lowcut, freqmax=highcut, corners=filt_order, zerophase=True)
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


def calculate_rf(path_ev, path_out, iterations=30, c1=10, c2=10, c3=1, c4=1, max_frequency=2.0, save=True):
    """
    Calculate receiver functions for waveforms trimmed around teleseismic arrivals.

    For choosing the c1-4 parameters values for the quality control please refer
    to the figure stored in the example directory or GH's PhD thesis.

    :type path_ev: str
    :param path_ev: Path where the waveforms are stored.
    :type path_out: str
    :param path_out: Path pointing to the directory for storing the RFs.
    :type iterations: int
    :param iterations: Number of iterations for the iterative time domain deconvolution.
    :type c1: float
    :param c1: Control parameters for quality criteria.
    :type c2: float
    :param c2: Control parameters for quality criteria.
    :type c3: float
    :param c3: Control parameters for quality criteria.
    :type c4: float
    :param c4: Control parameters for quality criteria.
    :type max_frequency: float
    :param max_frequency: High cut for bandpass filter in Hz (default is to 2.0 Hz).
    :type save: bool
    :param save:  Whether to save the figure or not (defaults to True)

    :returns: Receiver functions stored in SAC files.
    """

    # Get list of events - as they were prepared in 01_Get_Events.py
    all_event_dir = glob.glob(path_ev + '*')
    for event_dir in all_event_dir:
        print('Calculating RF for event in: ', event_dir)
        station_list = rf_util.get_list_of_unique_stations(event_dir)
        east_comp_traces = obspy.Stream()
        north_comp_traces = obspy.Stream()
        vert_comp_traces = obspy.Stream()
        for station in station_list:
            single_station_trace = obspy.read(event_dir + '/*' + station + '*')
            if len(single_station_trace) != 1 or len(single_station_trace) != 1 or len(single_station_trace) != 1:
                raise IOError('Either more or less traces than needed!')
            if single_station_trace[0].stats.channel[-1] == 'Z':
                vert_comp_traces.append(single_station_trace[0])
            if single_station_trace[0].stats.channel[-1] == 'E':
                east_comp_traces.append(single_station_trace[0])
            if single_station_trace[0].stats.channel[-1] == 'N':
                north_comp_traces.append(single_station_trace[0])
        # Quality control
        # Time before P-arrival time should be same for all traces!
        tbefore = vert_comp_traces[0].stats.sac.a
        # Sampling rate [Hz]
        fs = vert_comp_traces[0].stats.sampling_rate
        # Delta [s]
        delta = vert_comp_traces[0].stats.delta
        # List of booleans (if True do the calculations)
        quality_control_1 = rms_quality_control(vert_comp_traces, east_comp_traces, north_comp_traces,
                                                c1=c1, c2=c2, c3=c3, c4=c4)
        for i, vertical_trace in enumerate(vert_comp_traces):
            # Station name
            station_name = rf_util.printing_station_name(vertical_trace.stats.station, vertical_trace.stats.network)
            if quality_control_1[i]:
                # If quality control is successful (i.e., qc_control[i] == True)
                Z = vert_comp_traces[i].copy()
                T = east_comp_traces[i].copy()
                R = north_comp_traces[i].copy()
                # TODO: NEED TO ROTATE TO N, E, Z BEFORE WE ROTATE TO RT
                # Rotate traces to Vertical (Z), Radial (R) and Tangential (T) components
                T.data, R.data, Z.data = rotate_trace(east=east_comp_traces[i].data,
                                                      north=north_comp_traces[i].data,
                                                      vertical=vert_comp_traces[i].data,
                                                      baz=vert_comp_traces[i].stats.sac.baz)
                Z.stats.channel = 'HHZ'
                T.stats.channel = 'HHT'
                R.stats.channel = 'HHR'
                # Ready for processing
                processR = R.copy()
                processZ = Z.copy()
                RF = processR.copy()
                RF.stats.channel = 'RRF'
                RF.data, ds = rf_util.IterativeRF(trace_z=processZ, trace_r=processR, iterations=iterations,
                                                  iteration_plots=False, summary_plot=False)
                RFconvolve = RF.copy()
                RFconvolve = ConvGauss(spike_trace=RFconvolve, high_cut=max_frequency,
                                       delta=RFconvolve.stats.delta)
                RFconvolve.stats.sac.a = ds
                # RF quality control
                quality_control_2 = rf_quality_control(RFconvolve)
                # If qc_2 is True
                if quality_control_2:
                    processZ = Z.copy()
                    processT = T.copy()
                    # RFconvolve.plot()
                    TRF = processT.copy()
                    TRF.stats.channel = 'TRF'
                    TRF.data, ds = rf_util.IterativeRF(trace_z=processZ, trace_r=processT, iterations=iterations,
                                                       iteration_plots=False, summary_plot=False)
                    # TODO: add the option to plot stuff in the main function...
                    # IterativeRF provides a serie of spikes
                    # Here the serie of spikes is convolved by a gaussian bell
                    # whose width should match the highest frequency kept in the data
                    # TRFconvolve.plot()
                    TRFconvolve = TRF.copy()
                    TRFconvolve = ConvGauss(spike_trace=TRFconvolve, high_cut=max_frequency,
                                            delta=TRFconvolve.stats.delta)
                    print('>>> Station: ', station_name, ' -- Passed QC!')
                    # And save
                    if save:
                        rf_util.RFSave(Trace=RFconvolve, pathOUT=path_out)
                        rf_util.RFSave(Trace=TRFconvolve, pathOUT=path_out)
                else:
                    print('>>> Station: ', station_name, ' -- Failed on QC 2.')
                    continue
            else:
                print('>>> Station: ', station_name, ' -- Failed on QC 1.')
                continue
    return
