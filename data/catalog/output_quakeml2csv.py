with open('Teleseismic_earthquakes.csv', 'a') as of:
    writer = csv.DictWriter(of, fieldnames=["ID",
                                            "Origin_time",
                                            "Latitude (degrees)",
                                            "Longitude (degrees)", "Depth (km)",
                                            "Magnitude"]
                            , delimiter=',')
    writer.writeheader()

for i, event in enumerate(cat):
    id = event.comments[-1].text
    orig = event.origins[-1].time
    lat = event.origins[-1].latitude
    lon = event.origins[-1].longitude
    dep = event.origins[-1].depth/1000
    if event.magnitudes == []:
        mag = '-'
    else:
        mag = event.magnitudes[-1].mag

    with open('Teleseismic_earthquakes.csv', 'a') as of:
        of.write('{}, {}, {}, {}, {}, {}\n'.
                  format(id, orig, lat, lon, dep, mag))


