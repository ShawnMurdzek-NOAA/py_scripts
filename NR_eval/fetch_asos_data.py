"""
Fetch ASOS Data from Iowa State University Iowa Environmental Mesonet

Based on this script from GitHub: 
https://github.com/akrherz/iem/blob/main/scripts/asos/iem_scraper_example.py

Command-line arguments:
argv[1] = Start time (YYYYMMDD)
argv[2] = End time (YYYYMMDD)
argv[3] = Number of years to go back from the start and end time provided
argv[4] = Text file containing surface station IDs (1 per line)

NOTE: This script does not work on Hera or Jet. Instead, I run it locally on my work laptop and
upload the data to Hera or Jet.

"""

import datetime
import json
import os
import sys
import time

from urllib.request import urlopen

# Number of attempts to download data
MAX_ATTEMPTS = 6
# HTTPS here can be problematic for installs that don't have Lets Encrypt CA
SERVICE = "http://mesonet.agron.iastate.edu/cgi-bin/request/asos.py?"

def download_data(uri):
    """Fetch the data from the IEM
    The IEM download service has some protections in place to keep the number
    of inbound requests in check.  This function implements an exponential
    backoff to keep individual downloads from erroring.
    Args:
      uri (string): URL to fetch
    Returns:
      string data
    """
    attempt = 0
    while attempt < MAX_ATTEMPTS:
        try:
            data = urlopen(uri, timeout=300).read().decode("utf-8")
            if data is not None and not data.startswith("ERROR"):
                return data
        except Exception as exp:
            print(f"download_data({uri}) failed with {exp}")
            time.sleep(5)
        attempt += 1

    print("Exhausted attempts to download, returning empty data")
    return ""


def get_stations_from_filelist(filename):
    """Build a listing of stations from a simple file listing the stations.
    The file should simply have one station per line.
    """
    stations = []
    if not os.path.isfile(filename):
        print(f"Filename {filename} does not exist, aborting!")
        sys.exit()
    with open(filename, encoding="ascii") as fh:
        for line in fh:
            stations.append(line.strip())
    return stations


def get_stations_from_networks():
    """Build a station list by using a bunch of IEM networks."""
    stations = []
    states = """AK AL AR AZ CA CO CT DE FL GA HI IA ID IL IN KS KY LA MA MD ME
     MI MN MO MS MT NC ND NE NH NJ NM NV NY OH OK OR PA RI SC SD TN TX UT VA VT
     WA WI WV WY"""
    networks = []
    for state in states.split():
        networks.append(f"{state}_ASOS")

    for network in networks:
        # Get metadata
        uri = (
            "https://mesonet.agron.iastate.edu/"
            f"geojson/network/{network}.geojson"
        )
        data = urlopen(uri)
        jdict = json.load(data)
        for site in jdict["features"]:
            stations.append(site["properties"]["sid"])
    return stations


def download_alldata():
    """An alternative method that fetches all available data.
    Service supports up to 24 hours worth of data at a time."""
    # timestamps in UTC to request data for
    startts = datetime.datetime(2012, 8, 1)
    endts = datetime.datetime(2012, 9, 1)
    interval = datetime.timedelta(hours=24)

    service = SERVICE + "data=all&tz=Etc/UTC&format=comma&latlon=yes&"

    now = startts
    while now < endts:
        thisurl = service
        thisurl += now.strftime("year1=%Y&month1=%m&day1=%d&")
        thisurl += (now + interval).strftime("year2=%Y&month2=%m&day2=%d&")
        print(f"Downloading: {now}")
        data = download_data(thisurl)
        outfn = f"{now:%Y%m%d}.txt"
        with open(outfn, "w", encoding="ascii") as fh:
            fh.write(data)
        now += interval


def main(startts, endts, station_file_name):
    """Our main method"""
    # timestamps in UTC to request data for
    #startts = datetime.datetime(2012, 8, 1)
    #endts = datetime.datetime(2012, 9, 1)

    service = SERVICE + "data=all&tz=Etc/UTC&format=comma&latlon=yes&elev=yes&"

    service += startts.strftime("year1=%Y&month1=%m&day1=%d&")
    service += endts.strftime("year2=%Y&month2=%m&day2=%d&")

    # Two examples of how to specify a list of stations
    stations = get_stations_from_filelist(station_file_name)
    # stations = get_stations_from_filelist("mystations.txt")
    for station in stations:
        uri = f"{service}&station={station}"
        print(f"Downloading: {station}")
        data = download_data(uri)
        outfn = f"{station}_{startts:%Y%m%d%H%M}_{endts:%Y%m%d%H%M}.txt"
        with open(outfn, "w", encoding="ascii") as fh:
            fh.write(data)

start_first_yr = int(sys.argv[1])
end_first_yr = int(sys.argv[2])
for dyr in range(int(sys.argv[3])):
    startts = datetime.datetime.strptime(str(start_first_yr-(dyr*10000)), '%Y%m%d')
    endts = datetime.datetime.strptime(str(end_first_yr-(dyr*10000)), '%Y%m%d')
    print('start time = %s, end time = %s' % (startts.strftime('%Y%m%d'), endts.strftime('%Y%m%d')))
    main(startts, endts, sys.argv[4])
