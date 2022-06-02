"""
Set of functions for quality control of the data.

Original codes by Matteo Scarponi on 30.11.2021
Location: Chavannes-pres-renens, CH
Date: Jan 2022
Author: Konstantinos Michailos
"""
import numpy as np


def rms_quality_control(ztraces, etraces, ntraces, c1=10, c2=10):
    """
    Quality control applied on the original Z, E, N component seismic traces.

    Loops on the stations to check quality control of the single station and of the station compared to
    the other stations of the network. Computes signal and rms background for all the components for each station
    For more info we refer you to GH's PhD thesis (chapter 3.4).

    :type ztraces: obspy.core.stream.Stream
    :param ztraces: Vertical component waveform traces.
    :type etraces: obspy.core.stream.Stream
    :param etraces: Horizontal component waveform traces (East-West).
    :type ntraces: obspy.core.stream.Stream
    :param ntraces: Horizontal component waveform traces (North-South).
    :type c1: float
    :param c1: Control parameters for quality criteria.
    :type c2: float
    :param c2: Control parameters for quality criteria.

    :returns: Lists of booleans. If they are True means they have passed the
              quality control number 1.
    """

    # Time before P-arrival time - should be same for all traces! (or move it in the loop)
    time_before = ztraces[0].stats.sac.a
    # Sampling rate [Hz]		- should be same for all traces! (or move it in the loop)
    fs = ztraces[0].stats.sampling_rate
    # Delta [s] 				- should be same for all traces! (or move it in the loop)
    delta = ztraces[0].stats.delta
    i0 = int((time_before - 30) * fs)
    i1 = int((time_before - 5) * fs)
    i2 = int((time_before + 20) * fs)

    # Calculate rms values
    full_rms_z_list = []
    full_rms_e_list = []
    full_rms_n_list = []
    for i in range(len(ztraces)):
        try:
            np.max(ztraces[i].data[i0:i1 + 1])
            np.max(ztraces[i].data[i1:i2 + 1])
        except Exception as e:
            print(e)
            print(">>> Error while trying to use trace: ", ztraces[i], ", skipping this station...")
            # In case there is a gap in the data we append the following values that will fail the QC control
            full_rms_z_list.append(0.1)
            full_rms_n_list.append(0.1)
            full_rms_e_list.append(0.1)
            continue

        # Root mean square -rms- of the entire trace.
        full_rms_z = np.sqrt(np.mean(ztraces[i].data ** 2))
        # Root mean square of the entire trace.
        full_rms_e = np.sqrt(np.mean(etraces[i].data ** 2))
        # Root mean square of the entire trace.
        full_rms_n = np.sqrt(np.mean(ntraces[i].data ** 2))

        full_rms_z_list.append(full_rms_z)
        full_rms_e_list.append(full_rms_e)
        full_rms_n_list.append(full_rms_n)

    median_rms_z = np.median(full_rms_z)
    median_rms_e = np.median(full_rms_e)
    median_rms_n = np.median(full_rms_n)

    conditionals = []
    for i in range(len(full_rms_z_list)):
        e1 = (median_rms_e * c1 >= full_rms_e_list[i]) and (full_rms_e_list[i] >= median_rms_e / c2)
        n1 = (median_rms_n * c1 >= full_rms_n_list[i]) and (full_rms_n_list[i] >= median_rms_n / c2)
        z1 = (median_rms_z * c1 >= full_rms_z_list[i]) and (full_rms_z_list[i] >= median_rms_z / c2)
        conditionals.append(e1 and n1 and z1)
        # All three conditions need to be True to keep the waveforms for RF calculations
    return conditionals


def rf_quality_control(trace):
    """
    Second quality control step. Eliminate any RFs with weak signals. We apply this
    to the calculated receiver function (single trace).

    The time criteria (i.e., time_1 and time_2) help eliminate traces with common problems such as
    1) the P wave does not emerge from the background
    2) there are timing problems between the different components of the seismometer.
    The amplitude limits (i.e., amplitude_1 and amplitude_2) determined empirically throw away traces
    that have ringing issues.

    :type trace: obspy.core.trace.Trace
    :param trace: Receiver function trace RF.
    :returns: Boolean. If true the trace passed the quality control.
    """

    # First phase of qc RF RMS
    # Time before P-arrival time - should be same for all traces! (or move it in the loop)
    time_before = trace.stats.sac.a
    # Sampling rate [Hz]		- should be same for all traces! (or move it in the loop)
    fs = trace.stats.sampling_rate
    # Delta [s] 				- should be same for all traces! (or move it in the loop)
    delta = trace.stats.delta
    i0 = int((time_before - 30) * fs)
    i1 = int((time_before - 10) * fs)
    i2 = int((time_before + 2) * fs)
    i3 = int((time_before + 30) * fs)

    # Calculate rms values
    try:
        np.max(trace.data[i0:i1 + 1])
        np.max(trace.data[i1:i2 + 1])
    except Exception as e:
        print(e)
        print(">>> Error while trying to use trace: ", trace, ", skipping this station...")
        # In case there is a gap in the data we append the following values that will fail the QC control below
        rms_signal_z = 1.0
        max_background_z = 10.0

    # Root mean square of the background (between -30 and -10 seconds before the P arrival).
    rms_background_z = np.sqrt(np.mean((trace.data[i0:i1 + 1] ** 2)))
    # Root mean square of the signal (between 2 and 30 seconds after the P arrival).
    rms_signal_z = np.sqrt(np.mean((trace.data[i2:i3 + 1] ** 2)))
    #
    print(rms_signal_z, rms_background_z)
    z4 = (rms_signal_z / rms_background_z > 1.0)
    print(z4)
    # second phase of qc RF peak
    ds = trace.stats.sac.a
    iRF_edge = round(trace.stats.sampling_rate * 60)
    iRF = np.argmax(np.abs(trace.data[:iRF_edge]))

    time_1 = (iRF >= (ds - 0.0) * trace.stats.sampling_rate)
    time_2 = (iRF <= (ds + 2.0) * trace.stats.sampling_rate)
    amplitude_1 = (trace.data[iRF] >= 0.05)
    amplitude_2 = (trace.data[iRF] <= 0.8)

    # All the above need to be true
    qc_test = time_1 and time_2 and amplitude_1 and amplitude_2 and z4

    return qc_test


def sta_lta_quality_control(trace, sta, lta, high_cut, threshold):
    """
    Quality control step applying STA/LTA algorithm to throw away
    traces with weak signals.

    :type trace: obspy.core.trace.Trace
    :param trace: Waveform trace to run STA/LTA.
    :type sta: int
    :param sta: Short-time average (STA) window (default is 3 s).
    :type lta: int
    :param lta: Long-time average (STA) window (default is 50 s).
    :type threshold: float
    :param threshold: Threshold which when passed we keep the traces.
    :type high_cut: float
    :param high_cut: High cut for bandpass in Hz.

    :returns: Boolean. If true the trace passed the quality control.
    """

    from obspy.signal.trigger import classic_sta_lta

    df = trace.stats.sampling_rate
    tr = trace.copy()
    tr.filter("highpass", freq=high_cut)
    a = classic_sta_lta(tr, nsta=int(sta * df), nlta=int(lta * df))
    if max(a) < threshold:
        print('Low STA/LTA...')
        qc = False
    else:
        # print('Passed STA/LTA!')
        qc = True

    return qc


def RFQuality(Trace):
    """Original function. To be used for comparison of the modified codes to
    the original ones MS sent.
    NOTE: CHECKED THIS and we get the same results
    """

    ds=Trace.stats.sac.a

    iRF_edge=round(Trace.stats.sampling_rate*60)
    iRF = np.argmax(np.abs(Trace.data[:iRF_edge]))

    time1 = not (iRF <= (ds-0.6)*Trace.stats.sampling_rate)
    time2 = not (iRF >= (ds+2.1)*Trace.stats.sampling_rate)
    amp1 = (Trace.data[iRF] >= 0.01)
    amp2 = (Trace.data[iRF] <= 1.0)

    return time1,time2,amp1,amp2