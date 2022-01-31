"""
Get stationxml files with channel information from all the
AlpArray seismic sites

Location: Chavannes-pres-renens, CH
Date: Jan 2022
Author: Konstantinos Michailos
"""
from obspy.clients.fdsn import RoutingClient
from obspy.clients.fdsn import Client
from obspy import UTCDateTime, Inventory

client=RoutingClient('eida-routing')
networks = ['Z3', 'CH', 'IV', 'XT', 'BW', 'FR',
            'CR', 'CZ', 'GR', 'GU', 'HU', 'MN', 'NI',
            'OE', 'OX', 'RD', 'SI', 'SK', 'SL', 'TH']
inventory = Inventory()
for net in networks:
    starttime = UTCDateTime("2015-01-01")
    endtime = UTCDateTime("2019-01-02")
    try:
        inv = client.get_stations(network=net, station="*",
                                    level='channel',
                                    minlatitude=40.0, maxlatitude=51.0,
                                    minlongitude=0.0, maxlongitude=20.0,
                                    starttime=starttime,
                                    endtime=endtime)
    except:
        print(net)


client = Client("RESIF")
inv_fr = Inventory()
networks = ['FR', 'RD']
inv_fr += client.get_stations(network=net, station="*",
                          level='channel',
                          minlatitude=40.0, maxlatitude=51.0,
                          minlongitude=0.0, maxlongitude=20.0,
                          starttime=starttime, endtime=endtime)

# Plot
inv.plot(projection='local', continent_fill_color='0.9',
         water_fill_color='1.0', marker='v', size=20, label=True,
         color='#b15928', color_per_network=True,
         colormap='Paired', legend='upper left',
         show=True, outfile='AlpArray.png', fig=None)
# Store
inv.write(net + 'Alparray.xml',format='STATIONXML')
