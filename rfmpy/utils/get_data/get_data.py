"""
Script to get data from known networks on the IRIS DMC and write miniseed \
files in hour-long, single-channel files in year/julian day directories.

:Author: CJC
Modified by K. Michailos
"""


def get_data(network, starttime, endtime, outdir='.'):
    """
    Function to download all data for a given network.

    Will download all stations available between start and end times given. \
    It will create day-long files, and will download data as day-long segments.

    :type network: str
    :param network: Network code
    :type starttime: UTCDateTime
    :param starttime: Time to begin donloading from, will use the date
    :type endtime: UTCDateTime
    :param endtime: Time to end download, will use the date
    :type outdir: str
    :param outdir: Path to write to, will write in Y????/R???.01 directories \
        within this head directory.

    .. note:: This function is slow and doesn't cope with errors, suggest \
        using the mass-downloader.
    """
    import obspy
    import os
    if int(obspy.__version__.split('.')[0]) >= 1:
        from obspy.clients.fdsn import Client
    else:
        from obspy.fdsn import Client
    from obspy import UTCDateTime, Stream, read
    starttime = UTCDateTime(starttime)
    endtime = UTCDateTime(endtime)

    kdays = (endtime.datetime - starttime.datetime).days + 1
    for i in range(kdays):
        t1 = starttime + (86400 * i)
        t2 = t1 + 86400
        bulk_info = [(network, '*', '*', 'B*', t1, t2)]

        client = Client('IRIS', debug=True)
        st = client.get_waveforms_bulk(bulk_info)
        # st.plot()
        # Now you should split the data into traces by station and channel and
        # Save them into Y????/R???.01 directories, named something useful,
        # See the SAMBA_archive for details.
        for tr in st:
            chan_seq = (tr.stats.station, tr.stats.network,
                              tr.stats.location, tr.stats.channel,
                              tr.stats.starttime.strftime('%Y.%j'))
            f_name = '.'.join(chan_seq)
            path = outdir + tr.stats.starttime.strftime('/Y%Y/R%j.01/')
            if not os.path.isdir(path):
                os.makedirs(path)
            tr.write(path + f_name, format='MSEED', encoding='STEIM1')


def get_mass_data(network, starttime, endtime, outdir='.', lat=0, lon=0,
                  minrad=0, maxrad=180, providers=['IRIS']):
    """
    Use obspy's massdownloader to download lots of data, stores as day-long \
    miniseed files using IRIS DMC naming conventions.

    :type network: str
    :param network: Network code
    :type starttime: UTCDateTime
    :param starttime: Time to begin donloading from, will use the date
    :type endtime: UTCDateTime
    :param endtime: Time to end download, will use the date
    :type outdir: str
    :param outdir: Path to write to, will write in Y????/R???.01 directories \
        within this head directory.
    :type lat: float
    :param lat: Origin latitude
    :type lon: float
    :param lon: Origin longitude
    :type minrad: float
    :param minrad: Minumum radius in degrees for stations from lat/lon.
    :type maxrad: float
    :param maxrad: Maximum radius in degrees for stations from lat/lon
    :type providers: list
    :param providers: List of providers to query.  Default is IRIS. Can parse \
        an empty list and will query all providers, but slow.

    .. note:: Currently selects data using a circular domain, default is \
        set to entire globe so that a whole network can be downloaded. \
        Options left in function to allow repurposing.
    """
    def get_mseed_storage(network, station, location, channel,
                          starttime, endtime):
        """Function to define the naming for the stored file.

        .. note:: Can only have these arguments.  As such this needs to be a \
            function defined within the other function to have access to the \
            outdir variable.
        """
        import os
        # Returning True means that neither the data nor the StationXML file
        # will be downloaded.
        # If a string is returned the file will be saved in that location.
        path = os.path.join(outdir,
                            "%s/%s.%s.%s.%s.%s" % (starttime.
                                                   strftime('Y%Y/R%j.01'),
                                                   network, station, location,
                                                   channel, starttime.
                                                   strftime('%Y.%j')))
        if os.path.exists(path):
            return True
        return path

    import obspy
    from obspy.clients.fdsn.mass_downloader import CircularDomain, \
        Restrictions, MassDownloader

    domain = CircularDomain(latitude=lat, longitude=lon,
                            minradius=minrad, maxradius=maxrad)

    restrictions = Restrictions(
        # Get data for a whole year.
        starttime=starttime,
        endtime=endtime,
        # Chunk it to have one file per day.
        chunklength_in_sec=86400,
        # Considering the enormous amount of data associated with continuous
        # requests, you might want to limit the data based on SEED identifiers.
        # If the location code is specified, the location priority list is not
        # used; the same is true for the channel argument and priority list.
        network=network, station="NA390", location="*", channel="H*",
        # The typical use case for such a data set are noise correlations where
        # gaps are dealt with at a later stage.
        reject_channels_with_gaps=False,
        # Same is true with the minimum length. All data might be useful.
        minimum_length=0.0,
        # Guard against the same station having different names.
        minimum_interstation_distance_in_m=0.0,
        # Do not sanitize downloads, currently a bug
        sanitize=False,
        location_priorities=("", "01", "00", "EP", "S1", "S3", "02", "10",
                             "09", "08", "03", "04", "06", "07", "05", "20",
                             "T0", "2C", "40", "50"))


    # Restrict the number of providers if you know which serve the desired
    # data. If in doubt just don't specify - then all providers will be
    # queried.
    mdl = MassDownloader(providers=providers)
    mseed_storage = get_mseed_storage
    #  + "/{station}.{network}.{location}.{channel}.{starttime}")
    mdl.download(domain, restrictions, mseed_storage=mseed_storage,
                 stationxml_storage="stations", threads_per_client=5)


if __name__ == '__main__':
    import sys
    from obspy import UTCDateTime
    if not len(sys.argv) == 5:
        raise IOError('Insufficient arguments, need network code, ' +
                      'start-time and end-time, which need to be ' +
                      'convertable to UTCDateTime, and output directory.')
    network = sys.argv[1]
    starttime = UTCDateTime(sys.argv[2])
    endtime = UTCDateTime(sys.argv[3])
    outdir = sys.argv[4]
    get_mass_data(network, starttime, endtime, outdir)
