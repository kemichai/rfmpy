"""
...
"""
import glob

sta_name = []
number_of_RF_discarded = []
with open('bad_RFs.txt', 'r') as f:
    for line in f:
        ln = line.split(' ')
        sta_name.append(ln[0])
        number_of_RF_discarded.append(float(ln[-1]))

all_RF_files = glob.glob('number_of_rf*.txt')

sta_name_ = []
sta_lon_ = []
sta_lat_ = []
number_of_RFs = []
for file in all_RF_files:
    print(file)
    with open(file, 'r') as f:
        for line in f:
            ln = line.split(' ')
            sta_name_.append(ln[0])
            number_of_RFs.append(float(ln[-1]))
            sta_lon_.append(float(ln[2]))
            sta_lat_.append(float(ln[1]))

for i, sta_ in enumerate(sta_name_):
    for j, sta in enumerate(sta_name):
        if sta_ == sta:
            perc = (number_of_RF_discarded[j] / (number_of_RFs[i] + number_of_RF_discarded[j]  )) * 100
            print(sta, sta_lat_[i], sta_lon_[i], perc)
    if sta_ not in sta_name:
        print(sta_, sta_lat_[i], sta_lon_[i], str(-1))

