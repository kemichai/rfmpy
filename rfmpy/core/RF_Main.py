"""
Function for calculating RFs.

Note: Functions here call the ones from utils.

Original codes by Matteo Scarponi on 30.11.2021
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

def calculate_rf(path_ev, path_out, inventory, iterations=200, ds=30, c1=10, c2=10, c3=1, c4=1,
                 max_frequency=1.0, save=True, plot=True):
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
    :type ds: int
    :param ds: Seconds from zero to align P-wave arrival (default is 30 seconds).
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
    #TODO: update docstrings
    :returns: Receiver functions stored in SAC files.
    """

    all_event_dir = glob.glob(path_ev + '*')
    #############################################
    # TEST TO READ ONLY 2 months of data...
    # all_event_dir = []
    # all_event_dir_ = glob.glob(path_ev + 'P*')
    # for ev in all_event_dir_:
    #     ev_name = ev.split('/')[-1]
    #     yr = ev_name.split('.')[0]
    #     dd = int(ev_name.split('.')[1])
    #     if yr == 'P_2016' and dd <= 30:
    #         all_event_dir.append(ev)
    #############################################




    for event_dir in all_event_dir:
        print('Calculating RF for event in: ', event_dir)
        # Read waveform triplets (same network, station, channel, location, sps, ntps,)
        vert_comp_traces, north_comp_traces, east_comp_traces = rf_util.fetch_waveforms(event_dir)
        # Corrects misaligned components using info from stationxml files (i.e., azimuth and dip of each component)
        # for the correct epoch! This include borehole seismometers (e.g., BH2, BH3, HH1, HH2, etc).
        east_comp_traces_corr, north_comp_traces_corr, vert_comp_traces_corr = signal_processing.correct_orientations(
            st_east=east_comp_traces,
            st_north=north_comp_traces,
            st_vertical=vert_comp_traces, inventory=inventory)
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
                sta_lta_Z = qc.sta_lta_quality_control(Z_sta_lta, sta=3, lta=50, high_cut=1.0, threshold=2.5)
                sta_lta_R = qc.sta_lta_quality_control(R_sta_lta, sta=3, lta=50, high_cut=1.0, threshold=2.5)
                if sta_lta_R and sta_lta_Z:
                    R_filtered = R.copy()
                    Z_filtered = Z.copy()
                    T_filtered = T.copy()
                    # Processing (i.e., demean, taper, bandpass, cut waveforms)
                    R_filtered, T_filtered, Z_filtered = signal_processing.rf_processing(R_filtered, T_filtered,
                                                                                         Z_filtered, low_cut=0.05,
                                                                                         high_cut=1.0, order=2,
                                                                                         t_bef=40, t_aft=60)
                    # TODO: add 40 before and 40 s after... option in the main function
                    processR = R_filtered.copy()
                    processZ = Z_filtered.copy()
                    RF = processR.copy()
                    RF.stats.channel = 'RRF'
                    RF.data, RF_cc = rf_util.IterativeRF(trace_z=processZ, trace_r=processR, iterations=iterations,
                                                         tshift=ds, iteration_plots=False, summary_plot=plot)
                    # Store cc value in the SAC header (CC between R component and approximated R component).
                    RF.stats.sac.cc_value = RF_cc
                    RFconvolve = RF.copy()
                    RFconvolve = signal_processing.ConvGauss(spike_trace=RFconvolve, high_cut=max_frequency,
                                                             delta=RFconvolve.stats.delta)
                    RFconvolve.stats.sac.a = ds
                    # RF quality control
                    quality_control_2 = qc.rf_quality_control(RFconvolve, c3=c3, c4=c4)
                    # If qc_2 is True
                    if quality_control_2:
                        processZ = Z_filtered.copy()
                        processT = T_filtered.copy()
                        TRF = processT.copy()
                        TRF.stats.channel = 'TRF'
                        TRF.data, TR_cc = rf_util.IterativeRF(trace_z=processZ, trace_r=processT, iterations=iterations,
                                                              tshift=ds, iteration_plots=False, summary_plot=False)
                        TRF.stats.sac.cc_value = TR_cc
                        TRFconvolve = TRF.copy()
                        TRFconvolve = signal_processing.ConvGauss(spike_trace=TRFconvolve, high_cut=max_frequency,
                                                                  delta=TRFconvolve.stats.delta)
                        print('>>> Station: ', station_name, ' -- Passed QC 1!', ' -- Passed STA/LTA QC!',
                              ' -- Passed QC 2!')
                        # Save receiver functions
                        if save:
                            rf_util.store_receiver_functions(RFconvolve, path_out + 'RF/')
                            rf_util.store_receiver_functions(TRFconvolve, path_out + 'TRF/')
                    else:
                        print('>>> Station: ', station_name, ' -- Failed on QC 2.')
                        continue
                else:
                    print('>>> Station: ', station_name, ' -- Failed on STA/LTA.')
                    continue
            else:
                print('>>> Station: ', station_name, ' -- Failed on QC 1.')
                continue
    return
