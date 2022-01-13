"""
Function for calculating RFs.

Note: Functions here call the ones from utils.

Original codes by Matteo Scarponi on 30.11.2021
Location: Chavannes-pres-renens, CH
Date: Jan 2022
Author: Konstantinos Michailos
"""

import rfmpy.utils.RF_Util as rf_util
from rfmpy.utils.signal_processing import rotate_trace, remove_response, ConvGauss
from rfmpy.utils.qc import rms_quality_control, rf_quality_control
from obspy import read_inventory, read_events, UTCDateTime as UTC
import itertools
from pathlib import Path
import glob
import obspy
import numpy as np


def calculate_rf(path_ev, path_out, iterations=200, c1=10, c2=10, c3=1, c4=1, max_frequency=1.0, save=True):
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
                # TODO: NEED TO ROTATE TO ZNE BEFORE WE ROTATE TO RT
                # from obspy import read, read_inventory
                # st = read("/path/to/ffbx_unrotated_gaps.mseed")
                # inv = read_inventory("/path/to/ffbx.stationxml")
                # st.rotate(method="->ZNE", inventory=inv)
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
