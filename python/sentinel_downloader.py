#this code is used to download the sentinel 2 data and copied from
#https://github.com/jgomezdans/grabba_grabba_hey/blob/master/grabba_grabba_hey/sentinel_downloader.py
# and change a little bit for the qi  and whether data
from functools import partial
import hashlib
import os
import datetime
import sys
import xml.etree.cElementTree as ET
import re

import requests
from concurrent import futures


# hub_url = "https://scihub.copernicus.eu/dhus/search?q="
hub_url = "https://scihub.copernicus.eu/apihub/search?q="
MGRS_CONVERT = "http://legallandconverter.com/cgi-bin/shopmgrs3.cgi"
aws_url = 'http://sentinel-s2-l1c.s3.amazonaws.com/?delimiter=/&prefix=tiles/'
aws_url_dload = 'http://sentinel-s2-l1c.s3.amazonaws.com/'
requests.packages.urllib3.disable_warnings()


def get_mgrs(longitude, latitude):
    """A method that uses a website to infer the Military Grid Reference System
    tile that is used by the Amazon data buckets from the latitude/longitude

    Parameters
    -------------
    longitude: float
        The longitude in decimal degrees
    latitude: float
        The latitude in decimal degrees
    Returns
    --------
    The MGRS tile (e.g. 29TNJ)
    """
    r = requests.post(MGRS_CONVERT,
                      data=dict(latitude=longitude,
                                longitude=latitude, xcmd="Calc", cmd="gps"))
    for liner in r.text.split("\n"):
        if liner.find("<title>") >= 0:
            mgrs_tile = liner.replace("<title>", "").replace("</title>", "")
            mgrs_tile = mgrs_tile.replace(" ", "")
    try:
        return mgrs_tile[:5]  # This should be enough
    except NameError:
        return None

def parse_aws_xml(xml_text):

    tree = ET.ElementTree(ET.fromstring(xml_text))
    files_to_get = []
    for elem in tree.iter():
        for k in elem.getchildren():
            if k.tag.find ("Key") >= 0:
                if k.text.find ("tiles") >= 0:
                    files_to_get.append( k.text )
    return files_to_get

def aws_grabber(url, output_dir):
    output_fname = os.path.join(output_dir, url.split("tiles/")[-1])
   
    if not os.path.exists(os.path.dirname (output_fname)):
        # We should never get here, as the directory should always exist 
        # Note that in parallel, this can sometimes create a race condition
        # Groan
        os.makedirs (os.path.dirname(output_fname))
    with open(output_fname, 'wb') as fp:
        while True:
            try:
                #print url
                r = requests.get(url, stream=True)
                break
            except requests.execeptions.ConnectionError:
                time.sleep ( 240 )
        for block in r.iter_content(8192):
            fp.write(block)
    print "Done with %s" % output_fname
    return output_fname


def download_sentinel_amazon(longitude, latitude, start_date, output_dir,
                             end_date=None, n_threads=15):
    """A method to download data from the Amazon cloud """
    # First, we get hold of the MGRS reference...
    mgrs_reference = get_mgrs(longitude, latitude)
    utm_code = mgrs_reference[:2]
    lat_band = mgrs_reference[2]
    square = mgrs_reference[3:]

    front_url = aws_url + "%s/%s/%s" % (utm_code, lat_band, square)
    this_date = start_date
    one_day = datetime.timedelta(days=1)
    files_to_download = []
    if end_date is None:
        end_date = datetime.datetime.today()
    while this_date <= end_date:

        the_url = "{0}{1}".format(front_url, "/{0:d}/{1:d}/{2:d}/0/".format(
            this_date.year, this_date.month, this_date.day))
        # change to contain the qi file and auxilary file
        r1 = requests.get(the_url)
        r2 = requests.get(the_url+'qi/')
        r3 = requests.get(the_url+'auxiliary/')
        more_files = parse_aws_xml(r1.text) + parse_aws_xml(r2.text) + parse_aws_xml(r3.text)
        #print more_files
        if len(more_files) > 0:
            files_to_download.extend ( more_files )
        this_date += one_day
    the_urls = []
    
    for fich in files_to_download:
        the_urls.append(aws_url_dload + fich)
       
        ootput_dir = os.path.dirname ( os.path.join(output_dir, 
                                                    fich.split("tiles/")[-1]))
        if not os.path.exists ( ootput_dir ):
            os.makedirs ( ootput_dir )
    
    ok_files = []
    download_granule_patch = partial(aws_grabber, output_dir=output_dir)
    with futures.ThreadPoolExecutor(max_workers=n_threads) as executor:
        for fich in executor.map(download_granule_patch, the_urls):
            #print fich
            ok_files.append(fich)

if __name__ == "__main__":    # location = (43.3650, -8.4100)
    # input_start_date = "2015.01.01"
    # input_end_date = None

    # username = "guest"
    # password = "guest"

    # input_sensor = "S2"


    # output_dir = "/data/selene/ucfajlg/tmp/"
    # granules, retfiles = download_sentinel ( location, input_start_date,
    # input_sensor, output_dir )

    download_sentinel_amazon(43.3650, -8.4100, datetime.datetime(2016, 1, 1),
                             "/tmp/", end_date=datetime.datetime(2016, 1, 25) )
