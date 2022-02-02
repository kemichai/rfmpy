"""
Number of traces...
"""
import glob
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
                        hard_drive_dir + 'RF_data/DATA_RFAA_part_1/EASI/easi_data/',
                        hard_drive_dir + 'RF_data/DATA_RFAA_part_1/FRANCE/data_sort/',
                        hard_drive_dir + 'RF_data/DATA_RFAA_part_1/FRANCE/data/',
                        hard_drive_dir + 'RF_data/DATA_RFAA_part_1/North_Italy/events_fri_ven/',
                        hard_drive_dir + 'RF_data/DATA_RFAA_part_2/Austria/data_AAA_corrected/',
                        # hard_drive_dir + 'RF_data/DATA_RFAA_part_2/CIFAlps/data_YP2012/',
                        hard_drive_dir + 'RF_data/DATA_RFAA_part_2/data_DINAR/',
                        hard_drive_dir + 'RF_data/DATA_RFAA_part_2/HU_SK/data/',
                        hard_drive_dir + 'RF_data/DATA_RFAA_part_3/AARF/DATA_MOBST/data/',
                        hard_drive_dir + 'RF_data/DATA_RFAA_part_3/AARF/DATA_PERMST/data/',
                        hard_drive_dir + 'RF_data/DATA_RFAA_part_3/GERMANY/DE_AA_RF/DATA/data/',
                        hard_drive_dir + 'RF_data/CIFALPS/data_YP2012/',
                        hard_drive_dir + 'RF_data/INGV-Permanent-data/',
                        hard_drive_dir + 'RF_data/INGV-Temporary-data/data/',
]

unigue_events = []
for path in path_wavs_list_part1:
    all_event_dir = glob.glob(path + 'P*')
    for event_dir in all_event_dir:
        ev_name = event_dir.split('/')[-1]
        if ev_name not in unigue_events:
            unigue_events.append(ev_name)
unigue_events.sort()

unique_traces = []
for path in path_wavs_list_part1:
    all_event_dir = glob.glob(path + 'P*')
    for event_dir in all_event_dir:
        traces = glob.glob(event_dir + '/*SAC')
        for tr in traces:
            tr_name = tr.split('/')[-1]
            if tr_name not in unique_traces :
                unique_traces.append(tr_name)

print('Number of tele-events used: ', str(len(unigue_events)))
print('Number of traces used: ', str(len(unique_traces)))

unique_traces = []
for path in path_wavs_list_part1:
    all_event_dir = glob.glob(path + 'P*')
    for event_dir in all_event_dir:
        traces = glob.glob(event_dir + '/*SAC')
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
