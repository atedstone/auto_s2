#!/usr/bin/env python
"""

Requires setup_classifier.sh or equivalent to be run first.

Automatically download Sentinel-2 imagery and subset appropriately.

1. download_process_s2.py, just for download
2. gdalwarp spatial subset (hopefully, driver-dependent!)
3. gdalwarp convert to multiband compressed tif, NIR-RGB bands only
4. quick-look compressed pngs to check if worth downloading the real deal.

"""


import subprocess
import numpy as np
import os
import sys
import json
import configparser
import datetime as dt
import calendar
import matplotlib.pyplot as plt

import geoutils as gu

from sentinelsat import SentinelAPI

import sentinel2_tools

import argparse

parser = argparse.ArgumentParser(description='Automatic processing of specific Sentinel-2 imagery')
parser.add_argument('-t', '--tile', dest='tile', default=None, 
        help='Tile ID')
parser.add_argument('-n', '--today', dest='today', default=False, action='store_true', 
        help='Try download of today')

args = parser.parse_args()

tile = '22WFV'
# set start and end dates
today = dt.datetime.now()
startDate = today - dt.timedelta(days=10)
endDate = today
# xmin ymin xmax ymax in UTM22N
area = (620630, 7423106, 686230, 7436626)

print('Starting auto_sentinel2.py')

# Open API to Copernicus SciHub
cscihub_cred = configparser.ConfigParser()
cscihub_cred.read_file(open(os.environ['CSCIHUB_SECRET']))
chub_api = SentinelAPI(cscihub_cred.get('account','user'), 
        cscihub_cred.get('account','password'),
        'https://scihub.copernicus.eu/dhus')


# dates creates a tuple from the user-defined start and end dates
dates = (startDate, endDate)
print(dates)

# set path to save downloads
L1Cpath = os.environ['PROCESS_DIR']


products = sentinel2_tools.query_L1C(chub_api, L1Cpath, tile, dates, 50)
print(products)

for ix, row in products.iterrows():

    print('Filename', row.filename)

    if not os.path.exists(os.path.join(L1Cpath, row.filename)):

        prod_info = chub_api.download(row.uuid, directory_path=L1Cpath)
        print('\t prod_info:', prod_info)
        os.chdir(L1Cpath)
        unzip = ('unzip -o {file}').format(file=prod_info['path'])
        print(unzip)
        subprocess.call(unzip, shell=True)

        if prod_info['path'][:-3] != 'SAFE':
            prod_dir = prod_info['path'][:-4] #remove .zip from folder end
            prod_dir = prod_dir + '.SAFE'
        else:
            prod_dir = prod_info['path']

        # sen2cor = ('{pl}/L2A_Process {file} --resolution=10').format(pl=os.environ['SEN2COR_BIN'], 
        #     file=prod_dir) 
        # print(sen2cor)
        # subprocess.call(sen2cor, shell=True)

        gdal_sds = ('SENTINEL2_L1C:{file}/MTD_MSIL1C.xml:10m:EPSG_32622').format(file=prod_dir)
        #full_path = os.path.join(L1Cpath, prod_dir, gdal_sds)

        file_out = os.path.join(L1Cpath, 'processed', prod_info['title'] + '.tif')
        #subprocess.call('mkdir -p {f}'.format(f=file_out), shell=True)
        gdal = ('gdalwarp -of GTiff -co "COMPRESS=LZW" -tr 10 10 -te {xmin} {ymin} {xmax} {ymax} {file_in} {file_out}').format(
                xmin=area[0], ymin=area[1],
                xmax=area[2], ymax=area[3],
                file_in=gdal_sds, 
                file_out=file_out)
        print(gdal)
        subprocess.call(gdal, shell=True)

        gdaladdo = ('gdaladdo --config COMPRESS_OVERVIEW JPEG {file}').format(file=file_out)
        print(gdaladdo)
        subprocess.call(gdaladdo, shell=True)

        im = gu.Raster(file_out)
        plt.figure()
        plt.imshow(im.data[3,:,:], cmap='Greys_r')
        plt.title(file_out)
        plt.savefig(file_out + '.jpg', dpi=200)


        # Still need to make a quick look jpg ideally