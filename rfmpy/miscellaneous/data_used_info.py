"""
Number of traces...
"""
import platform
import glob
import obspy

# Set up parameters and paths
if platform.node().startswith('kmichailos-laptop'):
    data_root_dir = '/media/kmichailos/SEISMIC_DATA/Data_archive'
    codes_root_dir = '/home/kmichailos/Desktop/codes/bitbucket'
    desktop_dir = '/home/kmichailos/Desktop'
    hard_drive_dir = '/media/kmichailos/SEISMIC_DATA/'
else:
    data_root_dir = '/media/kmichall/SEISMIC_DATA/Data_archive'
    codes_root_dir = '/home/kmichall/Desktop/Codes/bitbucket'
    desktop_dir = '/home/kmichall/Desktop'
    hard_drive_dir = '/media/kmichall/SEISMIC_DATA/'

path_wavs_list_part1 = [
    hard_drive_dir + 'RF_data/DATA_RFAA_part_1/SWISS/data/',
    hard_drive_dir + 'RF_data/DATA_RFAA_part_1/EASI/data/',
hard_drive_dir + 'RF_data/DATA_RFAA_part_1/SLOVENIA/data/',
hard_drive_dir + 'RF_data/DATA_RFAA_part_2/OBS/data/',
hard_drive_dir + 'RF_data/DATA_RFAA_part_1/FRANCE/south_Fr_unsort/',
hard_drive_dir + 'RF_data/DATA_RFAA_part_1/FRANCE/data/',
hard_drive_dir + 'RF_data/DATA_RFAA_part_1/North_Italy/events_fri_ven/',
hard_drive_dir + 'RF_data/DATA_RFAA_part_2/Austria/data_AAA_corrected/',
hard_drive_dir + 'RF_data/DATA_RFAA_part_2/CIFAlps/data_YP2012/',
hard_drive_dir + 'RF_data/DATA_RFAA_part_2/data_DINAR/',
hard_drive_dir + 'RF_data/DATA_RFAA_part_2/HU_SK/data/',
hard_drive_dir + 'RF_data/DATA_RFAA_part_3/AARF/DATA_MOBST/data/',
hard_drive_dir + 'RF_data/DATA_RFAA_part_3/AARF/DATA_PERMST/data/',
hard_drive_dir + 'RF_data/DATA_RFAA_part_3/GERMANY/DE_AA_RF/DATA/data/',
hard_drive_dir + 'RF_data/CIFALPS/cifalps_unsort/',
hard_drive_dir + 'RF_data/INGV-Permanent-data/',
hard_drive_dir + 'RF_data/INGV-Temporary-data/data/',
hard_drive_dir + 'RF_data/AAPA/data/',]

unigue_events = []
for path in path_wavs_list_part1:
    all_event_dir = glob.glob(path + 'P*')
    for event_dir in all_event_dir:
        ev_name = event_dir.split('/')[-1]
        if ev_name not in unigue_events:
            unigue_events.append(ev_name)
unigue_events.sort()
print('Number of tele-events used: ', str(len(unigue_events)))


unique_traces = []
for path in path_wavs_list_part1:
    all_event_dir = glob.glob(path + 'P*')
    for event_dir in all_event_dir:
        traces = glob.glob(event_dir + '/*Z.SAC')
        for tr in traces:
            tr_name = tr.split('/')[-1]
            if tr_name not in unique_traces :
                unique_traces.append(tr_name)

print('Number of tele-events used: ', str(len(unigue_events)))
print('Number of traces used: ', str(len(unique_traces)))


stations = []
for path in path_wavs_list_part1:
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
            # if station_name not in stations:
            stations.append(station_name)

unique_all_sta = []
for s in stations:
    if s not in unique_all_sta:
        unique_all_sta.append(s)
# number of RFs on each station
for station in unique_all_sta:
    print(station, stations.count(station))
# using this for making the gmt plot



unique_traces = []
for path in path_wavs_list_part1:
    all_event_dir = glob.glob(path + 'P*')
    for event_dir in all_event_dir:
        traces = glob.glob(event_dir + '/*.SAC')
        for tr in traces:
            tr_name = tr.split('/')[-1]
            if tr_name not in unique_traces and tr_name.split('.')[-2][-1] is not 'E' and tr_name.split('.')[-2][-1] is not 'N' and tr_name.split('.')[-2][-1] is not 'Z':
                unique_traces.append(tr_name)

stations_without_NE = []
for trace in unique_traces:
    sta = trace.split('.')[-3]
    cha = trace.split('.')[-2]
    net = trace.split('.')[-4]
    stanam = net + '.' + sta + '.' + cha
    if stanam not in stations_without_NE:
        stations_without_NE.append(stanam)

e_trace_name = 'CH.WOLEN.BH2'
for net in inv:
    for sta in net:
        for cha in sta:
            cha_name = net.code + '.' + sta.code + '.' + cha.code
            if e_trace_name == cha_name:
                print(cha_name, cha.azimuth, cha.dip)


# CREATE a quakeml catalog for the teleseismic events... from the cut waveforms...
from obspy import UTCDateTime, read
from obspy.core.event import Catalog, Event, Origin, Magnitude
import obspy
import numpy as np
cat = Catalog()
cat.description = "Teleseismic events for RF calculations with AlpArray data"
unigue_events = []
for path in path_wavs_list_part1:
    all_event_dir = glob.glob(path + 'P*')
    for event_dir in all_event_dir:
        ev_name = event_dir.split('/')[-1]
        print(ev_name)
        if ev_name not in unigue_events:
            unigue_events.append(ev_name)
            traces = glob.glob(event_dir + '/*SAC')
            trace = read(traces[0])
            tr = trace[0]
            e = Event()
            o = Origin()
            m = Magnitude()

            o.time = tr.stats.starttime + 120
            o.latitude = tr.stats.sac.evla
            o.longitude = tr.stats.sac.evlo
            o.depth = tr.stats.sac.evdp * 1000
            e.comments.append(obspy.core.event.Comment(text="%s" % ev_name))
            m.mag = tr.stats.sac.mag
            m.magnitude_type = "M"

            e.origins = [o]
            e.magnitudes = [m]
            m.origin_id = o.resource_id

            cat.append(e)

cat.events.sort(key=lambda e: e.origins[-1].time)
cat.plot('global')
unigue_events.sort()
print('Number of tele-events used: ', str(len(unigue_events)))
cat.write('Teleseismic_events_RF.xml', format='QUAKEML')






# Compare RFs

