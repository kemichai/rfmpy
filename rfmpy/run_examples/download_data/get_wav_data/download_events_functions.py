#!/usr/bin/env python
import csv
import numpy as np
#from fdsnwsscripts import fdsnws_fetch
from subprocess import call
from obspy.core.utcdatetime import UTCDateTime
import os
from obspy.clients.fdsn import RoutingClient
from obspy.core import read
from obspy.geodetics import locations2degrees
from obspy.geodetics import gps2dist_azimuth
from obspy.taup import TauPyModel
from obspy import read_inventory
from obspy.io.sac import SACTrace
import time
from obspy.signal import rotate 
from obspy.clients.fdsn import Client
from obspy import Stream
import glob
import string

###Paste function for creating .csv
def pasteR(vector,sep=" "):
    s=str(vector[0])
    for i in np.arange(1,len(vector)):
        s+=sep
        s+=str(vector[i])
    return(s)

#####Read event csv
def read_eventcsv(path,minmag=5.5,cnames=True,useclient=False,cl="USGS",starttime=None,endtime=None):
    if useclient:
        client = Client(cl)
        cat = client.get_events(starttime=starttime,endtime=endtime,minmagnitude=minmag)
        evtimes=[]
        lats=[]
        lons=[]
        dps=[]
        mags=[]
        ids=[]
        for ev in cat:
            lats.append(ev.preferred_origin().latitude)
            lons.append(ev.preferred_origin().longitude)
            ttemp=ev.preferred_origin().time
            evtimes.append(ttemp)
            mags.append(ev.preferred_magnitude().mag)
            dps.append(ev.preferred_origin().depth/1000)
            if cnames:
                jday=("00"+str(ttemp.julday))[-3:]
                hr=("0"+str(ttemp.hour))[-2:]
                mn=("0"+str(ttemp.minute))[-2:]
                sec=("0"+str(ttemp.second))[-2:]
                line=pasteR([str(ttemp.year),jday,hr,mn,sec],sep=".")
                ids.append(line)
            else:
                ids.append(ev['extra']['eventid']['value'])
        ids=np.asarray(ids)
        ids=ids.astype("|S18")
        idi=np.arange(len(ids))
        for ind in idi:
            t=idi[ids == ids[ind]]
            tl=len(t)
            if tl > 1:
                for cnt,ind2 in enumerate(t):
                    ids[ind2]=str(ids[ind2]+string.ascii_uppercase[cnt])
        evmat=np.column_stack((ids,evtimes,lats,lons,dps,mags))

        if not len(ids) == len(np.unique(ids)):
            raise  ValueError("IDs are not unique")
        return(evmat,evtimes)


    with open(path, 'r') as csvfile:
        reader = csv.reader(csvfile, delimiter=',', quotechar='"')
        rs=[x for x in reader]
    nrs=np.asarray(rs)
    try:
        tind=np.arange(len(rs[0]))[nrs[0] == "time"][0]
        latind=np.arange(len(rs[0]))[nrs[0] == "latitude"][0]
        lonind=np.arange(len(rs[0]))[nrs[0] == "longitude"][0]
        dind=np.arange(len(rs[0]))[nrs[0] == "depth"][0]
        magind=np.arange(len(rs[0]))[nrs[0] == "mag"][0]
        idind=np.arange(len(rs[0]))[nrs[0] == "id"][0]
    except IndexError:
        raise ValueError("Required event column missing or mislabelled")
    if not cnames:
        try:
            idind=np.arange(len(rs[0]))[nrs[0] == "id"][0]
        except IndexError:
            raise ValueError("id event column missing or mislabelled. Set cnames=True")

    evtimes=np.asarray([UTCDateTime(x[tind]) for x in nrs[1:]])
    lats=np.asarray([float(x[latind]) for x in nrs[1:]])
    lons=np.asarray([float(x[lonind]) for x in nrs[1:]])
    dps=np.asarray([float(x[dind]) for x in nrs[1:]])
    mags=np.asarray([float(x[magind]) for x in nrs[1:]])
    if cnames:
        ids=[]
        for j in np.arange(len(evtimes)):
            jday=("00"+str(evtimes[j].julday))[-3:]
            hr=("0"+str(evtimes[j].hour))[-2:]
            mn=("0"+str(evtimes[j].minute))[-2:]
            sec=("0"+str(evtimes[j].second))[-2:]
            line=pasteR([str(evtimes[j].year),jday,hr,mn,sec],sep=".")
            ids.append(line)
            
    else:
        ids=np.asarray([x[idind] for x in nrs[1:]])

    evmat=np.column_stack((ids,evtimes,lats,lons,dps,mags))

    evmat=evmat[evmat[:,5] >= minmag]
    
    ids=np.asarray(ids)
    if not len(ids) == len(np.unique(ids)):
        raise  ValueError("IDs are not unique")
    return(evmat,evtimes)

#######Read station csv
def read_stationcsv(path,defaultnet="_ALPARRAY",usestatclient=False):
    if usestatclient:
        return([],[])
    if path == "*":
        return(["*"],[defaultnet])
    with open(path, 'r') as csvfile:
        reader = csv.reader(csvfile, delimiter=',', quotechar='"')
        ss=[x for x in reader]

    nss=np.asarray(ss)
    netalp=False
    statind=np.arange(len(nss[0]))[nss[0] == "station"][0]
    try:
        netind=np.arange(len(nss[0]))[nss[0] == "network"][0]
    except IndexError:
        netalp=True
        print("No network column found; defaulting to "+defaultnet)

    stations=[x[statind] for x in nss[1:]]
    if netalp:
        networks=np.repeat(defaultnet,len(stations))
    else:
        networks=[x[netind] for x in nss[1:]]
    return(stations,networks)


#####Populate wildcard
def populate(stations,networks,evtimes,routername="eida-routing",rclient=True,usestatclient=False,network=None,minlatitude=-90,minlongitude=-180,maxlatitude=90,maxlongitude=180,includeZS=False,c_inv=[]):
    outstations=[]
    outnetworks=[]
    if rclient:
        client = RoutingClient(routername)
    else:
        client= Client(routername)
    if usestatclient:
        if len(c_inv) == 0:
            inv=client.get_stations(network=network, station="*",starttime=min(evtimes),endtime=max(evtimes),includerestricted=True,level="station",minlatitude=minlatitude,minlongitude=minlongitude,maxlatitude=maxlatitude,maxlongitude=maxlongitude)
        else:
            inv=c_inv
        for net in inv:
            for stat in net:
                if (not net.code == "ZS") or includeZS:
                    outstations.append(stat.code)
                    outnetworks.append(net.code)
        return(outstations,outnetworks)

    for i in np.arange(len(stations)):
        if stations[i] == "*":
            if len(c_inv) == 0:
                try:
                    inv=client.get_stations(network=networks[i], station="*",starttime=min(evtimes),endtime=max(evtimes),includerestricted=True,level="station")
                except:
                    inv=[]
            else:
                inv=c_inv
            for net in inv:
                for stat in net:
                    if (not net.code == "ZS") or includeZS:
                        outstations.append(stat.code)
                        outnetworks.append(net.code)
        else:
            if (not networks[i] == "ZS") or includeZS:
                outstations.append(stations[i])
                outnetworks.append(networks[i])
    return(outstations,outnetworks)

####Read station metadata
def stat_meta(wd,stations,networks,evtimes,routername="eida-routing",rclient=True,mode="continue",write=True,c_inv=[]):
    if mode == "retry":
        return([],[],[],[])
    if mode == "all":
        pass
        #print("New download...")
    if mode == "continue":
        file = open(wd+"missing_stations","a+")
        file.close()
        skip=[]
        with open(wd+"missing_stations", 'r') as csvfile:
            reader = csv.reader(csvfile, delimiter=',', quotechar='"')
            ss=[skip.append(x) for x in reader]
        netstat=[x[1]+x[2] for x in skip]
        netstatin=[stations[x]+networks[x] for x in np.arange(len(stations))]
        stations=[stations[x] for x in np.arange(len(stations)) if netstatin[x] not in netstat]
        networks=[networks[x] for x in np.arange(len(networks)) if netstatin[x] not in netstat]
    if rclient:
        client = RoutingClient(routername)
    else:
        client=Client(routername)
    if len(c_inv) == 0:
        inv=read_inventory().select(station="BLANKINV12345")
        #inv=client.get_stations(network=networks[0], station=stations[0],starttime=min(evtimes),endtime=max(evtimes),includerestricted=True,level="channel")
        for i in np.arange(len(stations)):
            try:
                inv+=client.get_stations(network=networks[i], station=stations[i],starttime=min(evtimes),endtime=max(evtimes),includerestricted=True,level="channel")
            except:
                print(stations[i])
                print("No station information")
    else:
        inv=c_inv

    missing_stat=[]
    mstatlist=[]
    mnetlist=[]
    for st in np.arange(len(stations)):
        tinv=inv.select(station=stations[st])
        if len(tinv) == 0:
            missing_stat.append(["*",stations[st],networks[st],"missing_stat"])
            mstatlist.append(stations[st])
            mstatlist.append(networks[st])

    if mode == "continue":
        wm="a+"
    else:
        wm="w"
    if write:
        file = open(wd+"missing_stations",wm) 
        for l in missing_stat:
            file.write(pasteR(l,sep=",")+"\n")
        file.close()
    ostat=[]
    onet=[]
    for subnet in inv:
        for substat in subnet:
            ostat.append(substat.code)
            onet.append(subnet.code)
    return(inv,missing_stat,ostat,onet)





##Main download function

def dl_event(evline,wd,stations,networks,inv,component="CH",minepi=30,maxepi=95,ws=-10,we=50,sortby="event",mod="iasp91",phase="P",flo=0.03,fhi=2,fdsn=True,arclink_token=None,downsample=True,rotrt="ZNE->LQT",dcidpath=None,rotzne=False,znepath=None,routing=None,rclient=True):
    model = TauPyModel(model=mod)
    failure=[]
    failure3=[]
    run=[]
    ev=evline
    t=ev[1]
    lat=ev[2]
    lon=ev[3]
    d=ev[4]
    mag=ev[5]
    id=ev[0]
    bf_station_name=[]
    bf_station_network=[]
    bf_station_inc=[]
    bf_station_slowness=[]
    bf_station_baz=[]
    bf_station_epi=[]
    for j in np.arange(len(stations)):
        episkip=False
        subinv=inv.select(station=stations[j],network=networks[j],time=t)
        if len(subinv) > 0:
            net=subinv[0].code
            stat=subinv[0][0].code
            slat=subinv[0][0].latitude
            slon=subinv[0][0].longitude
            epi=locations2degrees(lat,lon,slat,slon)
            baz=gps2dist_azimuth(slat,slon,lat,lon)[1]
            if epi < minepi or epi > maxepi:
                failure.append([id,stations[j],networks[j],component,"epi_dist"])
                episkip=True
            else:
                ptime=model.get_travel_times(source_depth_in_km=d,distance_in_degree=epi,phase_list=[phase])
                if len(ptime) == 0:
                    failure.append([id,stations[j],networks[j],component,"epi_dist"])
                    episkip=True

        subinv=inv.select(station=stations[j],network=networks[j],time=t,channel=component+"*")
        if len(subinv) == 0 or episkip:
            if len(subinv) == 0 and not episkip:
                failure.append([id,stations[j],networks[j],component,"no_data"])
        else:
            if not (epi < minepi or epi > maxepi):
                #ptime=model.get_travel_times(source_depth_in_km=d,distance_in_degree=epi,phase_list=[phase])
                ptrav=ptime[0].time
                inc=ptime[0].incident_angle
                slw=ptime[0].ray_param_sec_degree
                ###add these to vectors
                bf_station_name.append(stations[j])
                bf_station_network.append(networks[j])
                bf_station_inc.append(inc)
                bf_station_slowness.append(slw)
                bf_station_baz.append(baz)
                bf_station_epi.append(epi)
                wstart=t+ptrav+ws
                wend=t+ptrav+we
                run.append([id,stat,net,wstart,wend,ptrav])
    reqname="data/"+id+str(len(stations))+str(float(UTCDateTime.now()))
    if len(run) > 0:
        ##Convert to np arrays
        bf_station_name=np.asarray(bf_station_name)
        bf_station_network=np.asarray(bf_station_network)
        bf_station_inc=np.asarray(bf_station_inc)
        bf_station_slowness=np.asarray(bf_station_slowness)
        bf_station_baz=np.asarray(bf_station_baz)
        bf_station_epi=np.asarray(bf_station_epi)
        bf_netstat=np.asarray([bf_station_network[i]+bf_station_name[i] for i in np.arange(len(bf_station_name))])

        file = open(wd+reqname,"w") 
        for st in run:
            sst=st[3].format_arclink()
            est=st[4].format_arclink()
            line=[sst,est,st[2],st[1],component+"*","*",".\n"]
            file.write(pasteR(line))
        file.close()
        failure2=[]
        rstats=[x[1] for x in run]
        rnets=[x[2] for x in run]
        if fdsn:
            cmd="fdsnws_fetch -f "+wd+reqname+" "+"-o"+" "+wd+reqname+".mseed"
        else:
            if dcidpath == None:
                cmd="arclink_fetch -k mseed4k -o "+wd+reqname+".mseed -u "+arclink_token+" -v "+wd+reqname
            else:
                cmd="arclink_fetch -k mseed4k -o "+wd+reqname+".mseed -u "+arclink_token+" -v "+wd+reqname+" -w "+dcidpath

        if routing == None:
            os.system(cmd)
        else:
            if rclient:
                client = RoutingClient(routing)
            else:
                client = Client(routing)
            strm=Stream()
            for st in run:
                try:
                    strm+=client.get_waveforms(network=st[2], station=st[1], location="*", channel=component+"*", starttime=st[3], endtime=st[4])
                except:
                    pass
            if len(strm) > 0:
                strm.write(wd+reqname+".mseed")  
        dne=False
        try:
            fsize=os.path.getsize(wd+reqname+".mseed")
        except:
            dne=True
            fsize=0
        try:
            ms=read(wd+reqname+".mseed",headeronly=True)
        except:
            dne=True
            fsize=0
        if fsize == 0:
            for rl in run:
                failure2.append([rl[0],rl[1],rl[2],component,"no_data"])
                failstats=[x[1] for x in failure2]
                failnets=[x[2] for x in failure2]
        else:
            ms=read(wd+reqname+".mseed")
            ms=merge_safe(ms)
            for tr in ms:
                if tr.stats.network == "":
                    tr2=[x for x in ms if ((not x.stats.network == "") and (x.stats.station==tr.stats.station))]
                    tr.stats.network=tr2[0].stats.network
            ustnet=np.unique([[x.stats.station,x.stats.network] for x in ms],axis=0)
            pustnet=np.unique([x.stats.station+x.stats.network for x in ms],axis=0)
            sts=ustnet[:,0]
            nets=ustnet[:,1]
            #SOMEHOW FDSN REPLACES SOME NETWORK NAMES WITH XX. THIS CHECKS AND CORRECTS BUT WILL FAIL IF TWO STATIONS HAVE THE SAME NAME AND ONE OR BOTH ARE CHANGED TO "XX" -__-
            nets=np.asarray(nets)
            rnets=np.asarray(rnets)
            rstats=np.asarray(rstats)
            sts=np.asarray(sts)
            for cn,el in enumerate(nets):
                if el == 'XX':
                    nets[cn]=rnets[rstats == sts[cn]][0]
            for tr in ms:
                if tr.stats.network == "XX":
                    tr.stats.network = rnets[rstats == tr.stats.station][0]
            failstats=[rstats[x] for x in np.arange(len(rstats)) if (rstats[x]+rnets[x] not in pustnet)]
            failnets=[rnets[x] for x in np.arange(len(rstats)) if (rstats[x]+rnets[x] not in pustnet)]
            for i in np.arange(len(failstats)):
                failure2.append([id,failstats[i],failnets[i],component,"no_data"])

            for l in np.arange(len(sts)):
                subst=sts[l]
                subms=ms.select(station=subst,network=nets[l])
                ##Get the proper inc and baz values
                bf_tf=bf_netstat == (nets[l]+subst)
                inc=bf_station_inc[bf_tf][0]
                epi=bf_station_epi[bf_tf][0]
                baz=bf_station_baz[bf_tf][0]
                slw=bf_station_slowness[bf_tf][0]
                #subms.merge()
                runline2=[x for x in run if x[2] == nets[l] and x[1] == subst][0]
                stt=runline2[3]
                ett=runline2[4]
                a=ws*-1
                o=a-runline2[5]
                if downsample:
                    for trz in subms:
                        srate=trz.stats['sampling_rate']
                        cdm=np.lcm(np.int(srate),20)
                        if srate > 20:
                            if cdm > srate:
                                trz.resample(cdm,no_filter=True)
                            fac=cdm/20
                            intfac=int(fac)
                            if intfac == fac:
                                trz.detrend()
                                trz.filter("lowpass",freq=10,zerophase=True)
                                trz.decimate(intfac,no_filter=True)
                #subms.trim(starttime=stt,endtime=ett)
                #subms._trim_common_channels()
                if not rotzne:
                    znes=False
                    rota=False
                if rotzne:
                    znes=False
                    try:
                        subms._trim_common_channels()
                        trmt=True
                    except:
                        trmt=False
                    with open(znepath, 'r') as csvfile:
                        reader = csv.reader(csvfile, delimiter=',', quotechar='"')
                        r=[x for x in reader]
                    rs=np.asarray(r)
                    crl=[x for x in rs if x[0] == nets[l] and x[1] == subst]
                    if len(crl) > 0:
                        rota=True
                    if len(crl) > 0 and trmt:
                        nazm=float(crl[0][2])
                        nev=float(crl[0][5])
                        if nev > 4 and len(subms) == 3:
                            try:
                                onetwo=False
                                zo=subms.select(channel="*Z")[0].data
                                no=subms.select(channel="*N")[0].data
                                eo=subms.select(channel="*E")[0].data
                            except:
                                try:
                                    onetwo=True
                                    zo=subms.select(channel="*Z")[0].data
                                    no=subms.select(channel="*1")[0].data
                                    eo=subms.select(channel="*2")[0].data
                                except:
                                    pass
                            try:
                                z,n,e=rotate.rotate2zne(zo,0,-90,no,nazm,0,eo,nazm+90,0)
                                if onetwo:
                                    subms.select(channel="*Z")[0].data=z
                                    subms.select(channel="*1")[0].data=n
                                    subms.select(channel="*1")[0].stats.channel=subms.select(channel="*1")[0].stats.channel[:2]+"N"
                                    subms.select(channel="*2")[0].data=e
                                    subms.select(channel="*2")[0].stats.channel=subms.select(channel="*2")[0].stats.channel[:2]+"E"
                                    znes=True
                                else:
                                    subms.select(channel="*Z")[0].data=z
                                    subms.select(channel="*N")[0].data=n
                                    subms.select(channel="*E")[0].data=e
                                    znes=True
                            except:
                                pass
                if not rotrt == None:
                    try:
                        subms._trim_common_channels()
                        trmt=True
                    except:
                        trmt=False
                    if trmt:
                        try:
                            subms.rotate(rotrt,back_azimuth=baz,inclination=inc)
                        except ValueError:
                            zstart=subms[0].stats.starttime
                            subms[1].stats.starttime=zstart
                            subms[2].stats.starttime=zstart
                            subms.rotate(rotrt,back_azimuth=baz,inclination=inc)
                subms.detrend()
                if not (flo == None or fhi == None):
                    subms.filter(type='bandpass',freqmin=flo,freqmax=fhi,zerophase=True,corners=2)
                subms.merge()
                subms.trim(starttime=stt,endtime=ett)
                subunstat=np.unique(np.asarray([x.stats.station for x in subms]))
                for subunst in subunstat:
                    try:
                        subms.select(station=subunst)._trim_common_channels()
                    except:
                        pass
                ts=[True for x in subms if x.stats.npts < .9*(we-ws)*x.stats.sampling_rate]
                nan=[True for x in subms if sum(np.isnan(x.data)) > 0]
                if sum(ts) > 0:
                    failure2.append([id,subst,nets[l],component,"missing_vals"])
                else:
                    if sum(nan) > 0:
                        failure2.append([id,subst,nets[l],component,"missing_vals"])
                if sortby == "event":
                    wp=wd+"data/"+phase+"_"+id+"/"
                if sortby == "station":
                    wp=wd+"data/"+phase+"_"+nets[l]+"."+subst+"/"
                if not os.path.exists(wp):
                    os.makedirs(wp)
                if len(subms) < 3:
                    failure2.append([id,subst,nets[l],component,"missing_vals"])
                for tr in subms:
                    sactr = SACTrace.from_obspy_trace(tr)
                    if not rotrt==None:
                        trinv=inv.select(station=tr.stats.station,network=tr.stats.network,time=t)
                    else:
                        trinv=inv.select(station=tr.stats.station,network=tr.stats.network,time=t,channel=tr.stats.channel)
                    if rotrt == None:
                        try:
                            sactr.cmpaz=trinv[0][0][0].azimuth
                            sactr.cmpinc=trinv[0][0][0].dip+90 ##convert to degrees from vertical
                        except:
                            pass
                    if znes:
                        sactr.user2=1
                        if tr.stats.channel[2] == "Z":
                            sactr.cmpaz=0
                            sactr.cmpinc=0
                        if tr.stats.channel[2] == "N":
                            sactr.cmpaz=0
                            sactr.cmpinc=90 ##convert to degrees from vertical
                        if tr.stats.channel[2] == "E":
                            sactr.cmpaz=90
                            sactr.cmpinc=90
                    else:
                        sactr.user2=0
                    if rota:
                        sactr.user3=1
                    else:
                        sactr.user3=0
                    sactr.stla=trinv[0][0].latitude
                    sactr.stlo=trinv[0][0].longitude
                    sactr.evdp=d
                    sactr.evlo=lon
                    sactr.evla=lat
                    sactr.mag=mag
                    sactr.a=a
                    sactr.ka=phase
                    sactr.kuser1=phase
                    sactr.o=o
                    sactr.user0=inc
                    sactr.user1=slw
                    sactr.baz=baz
                    sactr.gcarc=epi
                    sactr.stel=trinv[0][0].elevation
                    sactr.kevnm=runline2[0][2:18]
                    sactr.write(wp+id+"."+tr.stats.network+"."+tr.stats.station+"."+tr.stats.channel+".SAC")
                failure.append([id,subst,nets[l],component,"completed"])
        try:
            os.remove(wd+reqname+".mseed")
        except:
            pass
        os.remove(wd+reqname)

        for i in np.arange(len(failure2)):
            stat=failure2[i][1]
            subst=stat
            net=failure2[i][2]
            rline=[x for x in np.arange(len(rstats)) if stat == rstats[x] and net == rnets[x]]
            rline2=run[rline[0]]
            reqname="data/"+id+stat+net+str(float(UTCDateTime.now()))
            file = open(wd+reqname,"w") 
            sst=rline2[3].format_arclink()
            est=rline2[4].format_arclink()
            line=[sst,est,rline2[2],rline2[1],component+"*","*",".\n"]
            file.write(pasteR(line))
            file.close()
            if fdsn:
                cmd="fdsnws_fetch -f "+wd+reqname+" "+"-o"+" "+wd+reqname+".mseed"
            else:
                if dcidpath == None:
                    cmd="arclink_fetch -k mseed4k -o "+wd+reqname+".mseed -u "+arclink_token+" -v "+wd+reqname
                else:
                    cmd="arclink_fetch -k mseed4k -o "+wd+reqname+".mseed -u "+arclink_token+" -v "+wd+reqname+" -w "+dcidpath
            if routing == None:
                os.system(cmd)
            else:
                if rclient:
                    client = RoutingClient(routing)
                else:
                    client=Client(routing)
                try:
                    strm = client.get_waveforms(network=rline2[2], station=rline2[1],location= "*",channel= component+"*", starttime=rline2[3], endtime=rline2[4])
                    strm.write(wd+reqname+".mseed")
                except:
                    pass
            dne=False
            try:
                fsize=os.path.getsize(wd+reqname+".mseed")
            except:
                dne=True
                fsize=0
            try:
                subms=read(wd+reqname+".mseed",headeronly=True)
            except:
                dne=True
                fsize=0

            if fsize == 0:
                failure3.append([id,stat,net,component,"no_data"])
            else:
                subms=read(wd+reqname+".mseed")
                ##Get the proper inc and baz values
                bf_tf=bf_netstat == (net+subst)
                inc=bf_station_inc[bf_tf][0]
                epi=bf_station_epi[bf_tf][0]
                baz=bf_station_baz[bf_tf][0]
                slw=bf_station_slowness[bf_tf][0]
                subms=merge_safe(subms)
                for tr in subms:
                    #tr.stats.sampling_rate=np.round(tr.stats.sampling_rate,1)
                    if tr.stats.network == "":
                        tr2=[x for x in ms if ((not x.stats.network == "") and (x.stats.station==tr.stats.station))]
                        tr.stats.network=tr2[0].stats.network
                    if tr.stats.network == "XX":
                        tr.stats.network = net
                #subms.merge()
                runline2=[x for x in run if x[2] == net and x[1] == stat][0]
                stt=runline2[3]
                ett=runline2[4]
                a=ws*-1
                o=a-runline2[5]
                if downsample:
                    for trz in subms:
                        srate=trz.stats['sampling_rate']
                        cdm=np.lcm(np.int(srate),20)
                        if srate > 20:
                            if cdm > srate:
                                trz.resample(cdm,no_filter=True)
                            fac=cdm/20
                            intfac=int(fac)
                            if intfac == fac:
                                trz.detrend()
                                trz.filter("lowpass",freq=10,zerophase=True)
                                trz.decimate(intfac,no_filter=True)
                #subms.trim(starttime=stt,endtime=ett)
                #subms._trim_common_channels()
                if not rotzne:
                    znes=False
                    rota=False
                if rotzne:
                    znes=False
                    try:
                        subms._trim_common_channels()
                        trmt=True
                    except:
                        trmt=False
                    with open(znepath, 'r') as csvfile:
                        reader = csv.reader(csvfile, delimiter=',', quotechar='"')
                        r=[x for x in reader]
                    rs=np.asarray(r)
                    crl=[x for x in rs if x[0] == nets and x[1] == subst]
                    if len(crl) > 0:
                        rota=True
                    if len(crl) > 0 and trmt:
                        nazm=float(crl[0][2])
                        nev=float(crl[0][5])
                        if nev > 4 and len(subms) == 3:
                            try:
                                onetwo=False
                                zo=subms.select(channel="*Z")[0].data
                                no=subms.select(channel="*N")[0].data
                                eo=subms.select(channel="*E")[0].data
                            except:
                                try:
                                    onetwo=True
                                    zo=subms.select(channel="*Z")[0].data
                                    no=subms.select(channel="*1")[0].data
                                    eo=subms.select(channel="*2")[0].data
                                except:
                                    pass
                            try:
                                z,n,e=rotate.rotate2zne(zo,0,-90,no,nazm,0,eo,nazm+90,0)
                                if onetwo:
                                    subms.select(channel="*Z")[0].data=z
                                    subms.select(channel="*1")[0].data=n
                                    subms.select(channel="*1")[0].stats.channel=subms.select(channel="*1")[0].stats.channel[:2]+"N"
                                    subms.select(channel="*2")[0].data=e
                                    subms.select(channel="*2")[0].stats.channel=subms.select(channel="*2")[0].stats.channel[:2]+"E"
                                    znes=True
                                else:
                                    subms.select(channel="*Z")[0].data=z
                                    subms.select(channel="*N")[0].data=n
                                    subms.select(channel="*E")[0].data=e
                                    znes=True
                            except:
                                pass
                if not rotrt == None:
                    try:
                        subms._trim_common_channels()
                        trmt=True
                    except:
                        trmt=False
                    if trmt:
                        try:
                            subms.rotate(rotrt,back_azimuth=baz,inclination=inc)
                        except ValueError:
                            zstart=subms[0].stats.starttime
                            subms[1].stats.starttime=zstart
                            subms[2].stats.starttime=zstart
                            subms.rotate(rotrt,back_azimuth=baz,inclination=inc)
                subms.detrend()
                if not (flo == None or fhi == None):
                    subms.filter(type='bandpass',freqmin=flo,freqmax=fhi,zerophase=True,corners=2)
                subms.merge()
                subms.trim(starttime=stt,endtime=ett)
                subunstat=np.unique(np.asarray([x.stats.station for x in subms]))
                for subunst in subunstat:
                    try:
                        subms.select(station=subunst)._trim_common_channels()
                    except:
                        pass
                ts=[True for x in subms if x.stats.npts < 0.9*(we-ws)*x.stats.sampling_rate]
                nan=[True for x in subms if sum(np.isnan(x.data)) > 0]
                if sum(ts) > 0:
                    failure3.append([id,stat,net,component,"missing_vals"])
                else:
                    if sum(nan) > 0:
                        failure3.append([id,stat,net,component,"missing_vals"])
                if sortby == "event":
                    wp=wd+"data/"+phase+"_"+id+"/"
                if sortby == "station":
                    wp=wd+"data/"+phase+"_"+net+"."+subst+"/"
                if not os.path.exists(wp):
                    os.makedirs(wp)
                if len(subms) < 3:
                    failure3.append([id,stat,net,component,"missing_vals"])
                for tr in subms:
                    if not rotrt==None:
                        trinv=inv.select(station=tr.stats.station,network=tr.stats.network,time=t)
                    else:
                        trinv=inv.select(station=tr.stats.station,network=tr.stats.network,time=t,channel=tr.stats.channel)
                    sactr = SACTrace.from_obspy_trace(tr)
                    if rotrt == None:
                        try:
                            sactr.cmpaz=trinv[0][0][0].azimuth
                            sactr.cmpinc=trinv[0][0][0].dip+90 ##convert to degrees from vertical
                        except:
                            pass
                    if znes:
                        sactr.user3=1
                        if tr.stats.channel[2] == "Z":
                            sactr.cmpaz=0
                            sactr.cmpinc=0
                        if tr.stats.channel[2] == "N":
                            sactr.cmpaz=0
                            sactr.cmpinc=90 ##convert to degrees from vertical
                        if tr.stats.channel[2] == "E":
                            sactr.cmpaz=90
                            sactr.cmpinc=90
                    else:
                        sactr.user3=0
                    sactr.stla=trinv[0][0].latitude
                    sactr.stlo=trinv[0][0].longitude
                    sactr.evdp=d
                    sactr.evlo=lon
                    sactr.evla=lat
                    sactr.mag=mag
                    sactr.a=a
                    sactr.ka=phase
                    sactr.kuser1=phase
                    sactr.o=o
                    sactr.baz=baz
                    sactr.user0=inc
                    sactr.gcarc=epi
                    sactr.stel=trinv[0][0].elevation
                    sactr.kevnm=runline2[0][2:18]
                    sactr.user1=slw
                    sactr.write(wp+id+"."+tr.stats.network+"."+tr.stats.station+"."+tr.stats.channel+".SAC")
                failure.append([id,subst,net,component,"completed"])
            try:
                os.remove(wd+reqname+".mseed")
            except:
                pass
            os.remove(wd+reqname)
    for l in failure3:
        failure.append(l)
    return(failure)


def dl_BH_HH(evmat,wd,stations,networks,inv,component="CH",minepi=35,maxepi=95,ws=-10,we=50,sortby="event",mod="iasp91",phase="P",flo=0.03,fhi=2,mode="continue",fdsn=True,arclink_token=None,downsample=True,rotrt="ZNE->LQT",dcidpath=None,rotzne=False,znepath=None,routing=None,client_name="eida-routing",rclient=True,retry_network="*",includeZS=True):
    if mode == "retry":
        evtimes=np.asarray([x[1] for x in evmat])
        completed_list,failure_list=retry_download(wd,evmat,evtimes,minepi=minepi,maxepi=maxepi,ws=ws,we=we,sortby=sortby,flo=flo,fhi=fhi,mod=model,fdsn=fdsn,arclink_token=arclink_token,phase=phase,downsample=downsample,rotrt=rotrt,dcidpath=dcidpath,rotzne=rotzne,znepath=znepath,client_name=client_name,rclient=rclient,retry_network=retry_network,includeZS=includeZS)
        return(completed_list,failure_list)
    if mode == "all":
        print("Downloading all events...")
        file = open(wd+"completed_events","w")
        file.close()
        file = open(wd+"missing_events","w")
        file.close()
    if not os.path.exists(wd+"/data"):
        os.makedirs(wd+"/data")
    if mode == "continue":
        print("Continuing download...")
        skip=[]
        file = open(wd+"completed_events","a+")
        file.close()
        file = open(wd+"missing_events","a+")
        file.close()
        with open(wd+"completed_events", 'r') as csvfile:
            reader = csv.reader(csvfile, delimiter=',', quotechar='"')
            ss=[skip.append(x) for x in reader]
        with open(wd+"missing_events", 'r') as csvfile:
            reader = csv.reader(csvfile, delimiter=',', quotechar='"')
            ss=[skip.append(x) for x in reader]

    for evl in evmat:
        #evl=evmat[44]
        completed_list=[]
        failure_list=[]
        substations=stations
        subnet=networks
        if mode == "continue" and len(skip) > 0:
            subskip=np.asarray([x for x in skip if x[0] == evl[0] and (x[3] == "HH" or x[4] == "completed" or x[4] == "epi_dist")])
            if len(subskip) > 0:
                subnet=[networks[x] for x in np.arange(len(substations)) if substations[x] not in subskip[:,1]]
                substations=[x for x in substations if x not in subskip[:,1]]
        rb=dl_event(evl,wd=wd,stations=substations,networks=subnet,inv=inv,component="CH",minepi=minepi,maxepi=maxepi,ws=ws,we=we,sortby=sortby,flo=flo,fhi=fhi,fdsn=fdsn,arclink_token=arclink_token,phase=phase,downsample=downsample,rotrt=rotrt,dcidpath=dcidpath,rotzne=rotzne,znepath=znepath,routing=routing,rclient=rclient)
        blank=[completed_list.append(x) for x in rb if x[4] == "completed"]
        blank=[failure_list.append(x) for x in rb if not x[4] == "completed"]
        restat=[x[1] for x in rb if x[4] == "no_data" or x[4] == "missing_vals"]
        renet=[x[2] for x in rb if x[4] == "no_data" or x[4] == "missing_vals"]
        rh=dl_event(evl,wd=wd,stations=restat,networks=renet,inv=inv,component="HH",minepi=minepi,maxepi=maxepi,ws=ws,we=we,sortby=sortby,flo=flo,fhi=fhi,fdsn=fdsn,arclink_token=arclink_token,phase=phase,downsample=downsample,rotrt=rotrt,dcidpath=dcidpath,rotzne=rotzne,znepath=znepath,routing=routing,rclient=rclient)
        blank=[completed_list.append(x) for x in rh if x[4] == "completed"]
        blank=[failure_list.append(x) for x in rh if not x[4] == "completed"]
        wm="a+"
        file = open(wd+"missing_events",wm) 
        for l in failure_list:
            file.write(pasteR(l,sep=",")+"\n")
        file.close()

        file = open(wd+"completed_events",wm) 
        for l in completed_list:
            file.write(pasteR(l,sep=",")+"\n")
        file.close()
    return(completed_list,failure_list)

def retry_download(wd,evmat,evtimes,minepi=35,maxepi=95,ws=-10,we=50,sortby="event",mod="iasp91",phase="P",flo=0.03,fhi=2,fdsn=True,arclink_token=None,downsample=True,rotrt="ZNE->LQT",dcidpath=None,rotzne=False,znepath=None,routing=None,client_name="eida-routing",rclient=True,retry_network="*",includeZS=True):
    ###Attempt missing stations
    missing=[]
    with open(wd+"missing_stations", 'r') as csvfile:
        reader = csv.reader(csvfile, delimiter=',', quotechar='"')
        ss=[missing.append(x) for x in reader]
    stations_all=[x[1] for x in missing]
    networks_all=[x[2] for x in missing]
    if (retry_network == "*" or retry_network == "_ALPARRAY"):
        if includeZS:
            stations=stations_all
            networks=networks_all
        else:
            stations=[stations_all[x] for x in np.arange(len(stations_all)) if networks_all[x] != "ZS"]
            networks=[x for x in networks_all if x != "ZS"]
    else:
        stations=[stations_all[x] for x in np.arange(len(stations_all)) if networks_all[x] == retry_network]
        networks=[x for x in networks_all if x == retry_network]
##Read station metadata
    try:
        inventory,missing_stat,stations,networks=stat_meta(wd,stations,networks,evtimes=evtimes,mode="all",routername=client_name,rclient=rclient)
    except:
        stations=[]

    if len(stations) > 0:
        comp,fail=dl_BH_HH(evmat,wd=wd,stations=stations,networks=networks,inv=inventory,minepi=minepi,maxepi=maxepi,ws=ws,we=we,sortby=sortby,flo=flo,fhi=fhi,mode="continue",mod=model,phase=phase,fdsn=fdsn,arclink_token=arclink_token,downsample=downsample,rotrt=rotrt,dcidpath=dcidpath,rotzne=rotzne,znepath=znepath,routing=routing)
    
    comb=[pasteR([networks[x],stations[x]],sep=",")+"\n" for x in np.arange(len(stations))]
    comb_all=[pasteR([networks_all[x],stations_all[x]],sep=",")+"\n" for x in np.arange(len(stations_all))]
    wrt=[x for x in comb_all if not x in comb]
    file = open(wd+"missing_stations",'a+') 
    for l in wrt:
        file.write(l)
    file.close()
###Attempt missing events

    missing_events=[]
    with open(wd+"missing_events", 'r') as csvfile:
        reader = csv.reader(csvfile, delimiter=',', quotechar='"')
        if (retry_network == "*" or retry_network == "_ALPARRAY"):
            if includeZS:
                ss=[missing_events.append(x) for x in reader if x[4] != "epi_dist"]
            else:
                ss=[missing_events.append(x) for x in reader if (x[2] != "ZS" and x[4] != "epi_dist")]
        else:
            ss=[missing_events.append(x) for x in reader if (x[2] == retry_network and x[4] != "epi_dist")]

    epi_events=[]
    with open(wd+"missing_events", 'r') as csvfile:
        reader = csv.reader(csvfile, delimiter=',', quotechar='"')
        if (retry_network == "*" or retry_network == "_ALPARRAY"):
            if includeZS:
                ss=[epi_events.append(x) for x in reader if x[4] == "epi_dist"]
            else:
                ss=[epi_events.append(x) for x in reader if (x[4] == "epi_dist" or x[2] == "ZS")]
        else:
            ss=[epi_events.append(x) for x in reader if (x[4] == "epi_dist" or x[2] != retry_network)]

    mevname=np.unique(np.asarray([x[0] for x in missing_events]))

    evmat=[x for x in evmat if x[0] in mevname]


    stations=np.asarray([x[1] for x in missing_events])
    unstat=np.unique(np.asarray(stations),return_index=True)[1]
    networks=np.asarray([x[2] for x in missing_events])
    inventory,missing_stat,stations,networks=stat_meta(wd,stations[unstat],networks[unstat],evtimes=evtimes,mode="all",write=False,routername=client_name,rclient=rclient)

    missing_BH=[x for x in missing_events if x[3] == "CH"]
    missing_HH=np.asarray([x for x in missing_events if x[3] == "HH"])
    HHmerged=[x[0]+x[1]+x[2] for x in missing_HH]

    new_missing_list=[]
    for line in missing_BH:
        evl1=[x for x in evmat if x[0] == line[0]]
        if len(evl1) > 0:
            evl=evl1[0]
            restat=[line[1]]
            renet=[line[2]]
            gls=glob.glob(wd+"data/*"+evl[0]+"/*"+renet[0]+"."+restat[0]+"*")
            if not len(gls) == 3:
                rt=dl_event(evl,wd=wd,stations=restat,networks=renet,inv=inventory,component="CH",minepi=minepi,maxepi=maxepi,ws=ws,we=we,sortby=sortby,flo=flo,fhi=fhi,arclink_token=arclink_token,phase=phase,downsample=downsample,rotrt=rotrt,dcidpath=dcidpath,rotzne=rotzne,znepath=znepath,routing=routing,rclient=rclient,fdsn=fdsn)[0]
                new_missing_list.append(rt)
                if rt[4] == "completed":
                    file = open(wd+"completed_events","a+") 
                    file.write(pasteR(rt,sep=",")+"\n")
                    file.close()
                if (rt[4] == 'no_data' or rt[4] == 'missing_vals') and (rt[0]+rt[1]+rt[2]) in HHmerged:
                    rtt=dl_event(evl,wd=wd,stations=restat,networks=renet,inv=inventory,component="HH",minepi=minepi,maxepi=maxepi,ws=ws,we=we,sortby=sortby,flo=flo,fhi=fhi,arclink_token=arclink_token,phase=phase,downsample=downsample,rotrt=rotrt,dcidpath=dcidpath,rotzne=rotzne,znepath=znepath,routing=routing,rclient=rclient,fdsn=fdsn)[0]
                    new_missing_list.append(rtt)
                    if rtt[4] == "completed":
                        file = open(wd+"completed_events","a+") 
                        file.write(pasteR(rtt,sep=",")+"\n")
                        file.close()
            else:
                print("Event already downloaded!")
                file = open(wd+"completed_events","a+") 
                file.write(evl[0]+","+restat[0]+","+renet[0]+","+gls[0][-7:-5]+",completed\n")
                file.close()

    failure_list = [x for x in new_missing_list if not x[4] == 'completed']
    completed_list = [x for x in new_missing_list if x[4] == 'completed']

    file = open(wd+"missing_events","w") 
    for l in failure_list:
        file.write(pasteR(l,sep=",")+"\n")
    file.close()
    file = open(wd+"missing_events","a+") 
    for l in epi_events:
        file.write(pasteR(l,sep=",")+"\n")
    file.close()
    
    return(failure_list,completed_list)



#removes traces with bad sampling rates before merging
def merge_safe(ms):
    for tr in ms:
        tr.stats.sampling_rate=np.round(tr.stats.sampling_rate,1)
    statnet=np.asarray([x.stats.network+x.stats.station for x in ms])
    u_statnet=np.unique(statnet)
    for u_substat in u_statnet:
        u_tf=statnet == u_substat
        p_stream=Stream()
        for itf,tf in enumerate(u_tf):
            if tf:
                p_stream+=ms[itf]
        n_s=np.asarray([x.stats.npts for x in p_stream])
        s_r=np.asarray([x.stats.sampling_rate for x in p_stream])
        l_sr=[]
        for u_s_r in np.unique(s_r):
            l_sr.append(np.sum(n_s[s_r == u_s_r]))
        l_sr=np.asarray(l_sr)
        best_sr=np.unique(s_r)[l_sr == np.max(l_sr)][0]
        for s_tr in p_stream:
            if not s_tr.stats.sampling_rate == best_sr:
                s_tr.stats.delete=True
    for tr in ms:
        if "delete" in tr.stats:
            if tr.stats.delete:
                ms.remove(tr)
    ms.merge(fill_value="interpolate")
    return(ms)


def verify_missing(wd):
    missing_events=[]
    with open(wd+"missing_events", 'r') as csvfile:
        reader = csv.reader(csvfile, delimiter=',', quotechar='"')
        ss=[missing_events.append(x) for x in reader]

    completed_events=[]
    with open(wd+"completed_events", 'r') as csvfile:
        reader = csv.reader(csvfile, delimiter=',', quotechar='"')
        ss=[completed_events.append(x) for x in reader]

    comp_merge=np.asarray([completed_events[x][0]+completed_events[x][1]+completed_events[x][2]+completed_events[x][3] for x in np.arange(len(completed_events))])
    miss_merge=np.asarray([missing_events[x][0]+missing_events[x][1]+missing_events[x][2]+missing_events[x][3] for x in np.arange(len(missing_events))])

    new_miss=[]
    for i,x in enumerate(miss_merge):
        if not x in comp_merge:
            ln=missing_events[i]
            new_miss.append(ln)

    file = open(wd+"missing_events","w") 
    for l in new_miss:
        file.write(pasteR(l,sep=",")+"\n")
    file.close()




###And run the program

#No longer use arclink_fetch
fdsn=True

#Depreciated arclink parameters
#arclink_fetch token used for accessing restricted data (used if fdsn=False) or email for other networks (requires dcidpasswords.txt to be set)
arclink_token=None
##dcidpasswords location (if data requires decryption
dcidpath='/data/home/mroczek/dcidpasswords.txt'


#Check that none of the missing_events haven't already been completed
#This should only be run if you suspect that missing_events is out of date (e.g. because retry mode failed before finishing)
#if mode == "retry":
#    verify_missing(wd)

###Read events csv
##
evmat,evtimes=read_eventcsv(eventcsv,minmag=minmag,cnames=cnames,useclient=useclient,cl=cl,starttime=starttime,endtime=endtime)

###Read stations csv
stations,networks=read_stationcsv(stationcsv,usestatclient=usestatclient)
##

###Populate * wild card
stations,networks=populate(stations,networks,evtimes,usestatclient=usestatclient,network=network,minlatitude=minlatitude,minlongitude=minlongitude,maxlatitude=maxlatitude,maxlongitude=maxlongitude,includeZS=includeZS,routername=client_name,rclient=rclient,c_inv=c_inv)
##

###Read station metadata
inventory,missing_stat,stations,networks=stat_meta(wd,stations,networks,evtimes=evtimes,mode=mode,routername=client_name,rclient=rclient,c_inv=c_inv)
##


###Begin download
comp,fail=dl_BH_HH(evmat,wd=wd,stations=stations,networks=networks,inv=inventory,minepi=minepi,maxepi=maxepi,ws=ws,we=we,sortby=sortby,flo=flo,fhi=fhi,mode=mode,mod=model,fdsn=fdsn,arclink_token=arclink_token,phase=phase,downsample=downsample,rotrt=rotrt,dcidpath=dcidpath,rotzne=rotzne,znepath=znepath,client_name=client_name,rclient=rclient,retry_network=network,includeZS=includeZS)
##

