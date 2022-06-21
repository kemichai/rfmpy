"""
Function for calculating RFs.

TODO: Functions here call the ones from utils. CHANGE THIS...

Location: Chavannes-pres-renens, CH
Date: Jan 2022
Author: Konstantinos Michailos
"""

import rfmpy.utils.RF_Util as rf_util
from rfmpy.utils import signal_processing
from rfmpy.utils import qc
from obspy import read_inventory, read_events, UTCDateTime as UTC
import itertools
from pathlib import Path
import glob
import obspy
import numpy as np
import os
import logging
import sys

# Create a log file
a_logger = logging.getLogger()
a_logger.setLevel(logging.INFO)
output_file_handler = logging.FileHandler("logfile_SWISS.txt")
stdout_handler = logging.StreamHandler(sys.stdout)
a_logger.addHandler(output_file_handler)
a_logger.addHandler(stdout_handler)


def calculate_rf(path_ev, path_out, inventory, iterations=200, ds=30,
                 c1=10, c2=10,
                 max_frequency=1.0,
                 sta_lta_qc=None,
                 pre_processing=None,
                 save=True, plot=True):
    """
    Calculate receiver functions for waveforms cut around teleseismic arrivals.

    For choosing the c1-c4 parameters values for the quality control please refer
    to the figure stored in the example directory or GH's PhD thesis.

    :type path_ev: str
    :param path_ev: Path where the waveforms are stored.
    :type path_out: str
    :param path_out: Path pointing to the directory for storing the RFs.
    :type inventory: obspy.core.inventory.Inventory
    :param inventory: Inventory containing response information for the stations.
    :type iterations: int
    :param iterations: Number of iterations for the iterative time domain deconvolution (default is 200).
    :type ds: int
    :param ds: Seconds from zero to align P-wave arrival during deconvolution (default is 30 seconds).
    :type c1: float
    :param c1: Control parameters for quality criteria.
    :type c2: float
    :param c2: Control parameters for quality criteria.
    :type max_frequency: float
    :param max_frequency: High cut for bandpass filter in Hz (default is to 1.0 Hz).
    :type sta_lta_qc: tuple
    :param sta_lta_qc: Tuple defining the sta/lta parameters for the qc step (NEEDS TO BE DEFINED; default is None).
    :type pre_processing: tuple
    :param pre_processing: Tuple defining some pre-processing parameters (NEEDS TO BE DEFINED; default is None).
    :type save: bool
    :param save: Whether to save the traces or not (defaults to True).
    :type plot: bool
    :param plot: Whether to plot the RFs or not (default is True).


    :returns: Receiver functions stored in SAC files.

    .. Note::
        Tuples for sta_lta_qc and pre_processing NEED TO BE DEFINED while calling the function!
        See example below.

    .. rubric:: Example

    >>> sta_lta_qc_parameters = {'sta': 3, 'lta': 50, 'high_cut': 1.0, 'threshold': 2.5}
    >>> pre_processing_parameters = {'low_cut': 0.05, 'high_cut': 1.0, 'order': 2, 't_before': 40, 't_after': 60}
    >>> calculate_rf(path_ev, path_out, inventory, iterations=200,
    ...              ds=30, c1=10, c2=10, c3=1, c4=1,
    ...              sta_lta_qc=sta_lta_qc_parameters,
    ...              pre_processing=pre_processing_parameters,
    ...              max_frequency=1, save=True, plot=False)
    """

    all_event_dir = glob.glob(path_ev + '*')
    for event_dir in all_event_dir:
        # print('Calculating RF for event in: ', event_dir)
        a_logger.info(f'Calculating RF for event in: {event_dir}')
        # Read waveform triplets (same network, station, channel, location, sps, ntps,)
        vert_comp_traces, north_comp_traces, east_comp_traces = rf_util.fetch_waveforms(event_dir)
        if len(vert_comp_traces) == 0 or len(north_comp_traces) == 0 or len(east_comp_traces) == 0:
            print('>>> No data skip this event.')
            continue
        # Corrects misaligned components using info from stationxml files (i.e., azimuth and dip of each component)
        # for the correct epoch! This include borehole seismometers (e.g., BH2, BH3, HH1, HH2, etc).
        east_comp_traces_corr, \
        north_comp_traces_corr, \
        vert_comp_traces_corr = signal_processing.correct_orientations(st_east=east_comp_traces,
                                                                       st_north=north_comp_traces,
                                                                       st_vertical=vert_comp_traces,
                                                                       inventory=inventory,
                                                                       logfile=a_logger)
        # Quality control -- List of booleans (if True do the calculations)
        quality_control_1 = qc.rms_quality_control(vert_comp_traces_corr,
                                                   east_comp_traces_corr,
                                                   north_comp_traces_corr,
                                                   c1=c1, c2=c2)
        for i, vertical_trace in enumerate(vert_comp_traces_corr):
            station_name = rf_util.printing_station_name(vertical_trace.stats.station, vertical_trace.stats.network)
            if quality_control_1[i]:
                # If quality control is successful (i.e., qc_control[i] == True)
                Z = vert_comp_traces_corr[i].copy()
                T = east_comp_traces_corr[i].copy()
                R = north_comp_traces_corr[i].copy()
                # Rotate traces to Vertical (Z), Radial (R) and Tangential (T) components
                T.data, R.data, Z.data = signal_processing.rotate_trace(east=east_comp_traces_corr[i].data,
                                                                        north=north_comp_traces_corr[i].data,
                                                                        vertical=vert_comp_traces_corr[i].data,
                                                                        baz=vert_comp_traces_corr[i].stats.sac.baz)
                Z.stats.channel = 'HHZ'
                T.stats.channel = 'HHT'
                R.stats.channel = 'HHR'
                Z_sta_lta = Z.copy()
                R_sta_lta = R.copy()
                sta_lta_Z = qc.sta_lta_quality_control(Z_sta_lta,
                                                       sta=sta_lta_qc['sta'],
                                                       lta=sta_lta_qc['lta'],
                                                       high_cut=sta_lta_qc['high_cut'],
                                                       threshold=sta_lta_qc['threshold'])
                sta_lta_R = qc.sta_lta_quality_control(R_sta_lta,
                                                       sta=sta_lta_qc['sta'],
                                                       lta=sta_lta_qc['lta'],
                                                       high_cut=sta_lta_qc['high_cut'],
                                                       threshold=sta_lta_qc['threshold'])
                if sta_lta_R and sta_lta_Z:
                    R_filtered = R.copy()
                    Z_filtered = Z.copy()
                    T_filtered = T.copy()
                    # Processing (i.e., demean, taper, bandpass, cut waveforms)
                    R_filtered,\
                    T_filtered,\
                    Z_filtered = signal_processing.pre_processing(R_filtered, T_filtered, Z_filtered,
                                                                  low_cut=pre_processing['low_cut'],
                                                                  high_cut=pre_processing['high_cut'],
                                                                  order=pre_processing['order'],
                                                                  t_bef=pre_processing['t_before'],
                                                                  t_aft=pre_processing['t_after'])
                    processR = R_filtered.copy()
                    processZ = Z_filtered.copy()
                    RF = processR.copy()
                    RF.stats.channel = 'RRF'
                    RF.data, RF_cc = rf_util.IterativeRF(trace_z=processZ, trace_r=processR, iterations=iterations,
                                                         tshift=ds, iteration_plots=False, summary_plot=plot)
                    # Store cc value in the SAC header (CC between observed R component and the predicted one).
                    RFconvolve = RF.copy()
                    RFconvolve = signal_processing.ConvGauss(spike_trace=RFconvolve, high_cut=max_frequency,
                                                             delta=RFconvolve.stats.delta)
                    RFconvolve.stats.sac.a = ds
                    RFconvolve.stats.sac.cc_value = RF_cc
                    # RF quality control
                    quality_control_2 = qc.rf_quality_control(RFconvolve)
                    # If qc_2 is True
                    if quality_control_2:
                        processZ = Z_filtered.copy()
                        processT = T_filtered.copy()
                        # TRF = processT.copy()
                        # TRF.stats.channel = 'TRF'
                        # TRF.data, TR_cc = rf_util.IterativeRF(trace_z=processZ, trace_r=processT, iterations=iterations,
                        #                                       tshift=ds, iteration_plots=False, summary_plot=False)
                        # TRF.stats.sac.cc_value = TR_cc
                        # TRFconvolve = TRF.copy()
                        # TRFconvolve = signal_processing.ConvGauss(spike_trace=TRFconvolve, high_cut=max_frequency,
                        #                                           delta=TRFconvolve.stats.delta)
                        # print('>>> Station: ', station_name, ' -- Passed QC 1!', ' -- Passed STA/LTA QC!',
                        #       ' -- Passed QC 2!')
                        a_logger.info(f'>>> Station: {station_name} - Passed QC 1! - Passed STA/LTA QC! - Passed QC 2!')
                        # Save receiver functions
                        if save:
                            rf_util.store_receiver_functions(RFconvolve, path_out + 'RF/')
                            # rf_util.store_receiver_functions(TRFconvolve, path_out + 'TRF/')
                    else:
                        # print('>>> Station: ', station_name, ' -- Failed on QC 2.')
                        a_logger.info(f'>>> Station: {station_name} - Failed on QC 2.')
                        continue
                else:
                    # print('>>> Station: ', station_name, ' -- Failed on STA/LTA.')
                    a_logger.info(f'>>> Station: {station_name} - Failed on STA/LTA.')
                    continue
            else:
                # print('>>> Station: ', station_name, ' -- Failed on QC 1.')
                a_logger.info(f'>>> Station: {station_name} - Failed on QC 1.')
                continue
    return
