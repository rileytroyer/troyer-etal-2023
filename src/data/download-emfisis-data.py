""" Script to download Van Allen Probe EMFISIS data for the times during injections.

@author Riley Troyer
science@rileytroyer.com
"""

####################### Initialize Program #######################
# Libraries
from bs4 import BeautifulSoup
from datetime import datetime
import logging
import numpy as np
import multiprocessing
import os
import pickle
import requests
import wget



# Initiate logging
logging.basicConfig(filename = f'logs/download-emfisis-data-{datetime.today().date()}.log',
                    encoding='utf-8',
                    format='%(asctime)s %(levelname)-8s %(message)s',
                    level=logging.INFO,
                    datefmt='%Y-%m-%d %H:%M:%S')

####################### End Initializing #######################


####################### Local Functions #######################
def get_url_paths(url:str, ext:str='', params:dict={}) -> list:
    """ Function to extract file names from https directory
    Gets files in url directory with ext extension
    Does this by parsing the html text from the webpage.
    INPUT
    url- url of directory to get files from 
    ext- extension of the files
    OUTPUT
    parent- list of all file pathnames within directory
    """

    response = requests.get(url, params=params)
    if response.ok:
        response_text = response.text
    else:
        return response.raise_for_status()
    soup = BeautifulSoup(response_text, 'html.parser')
    parent = [url + node.get('href') for node in soup.find_all('a')
            if node.get('href').endswith(ext)]
    return parent
####################### End of Local Functions #######################


####################### START OF PROGRAM #######################

# Define directories to save the various data
mag_save_dir = 'data/raw/mag-waveform/'
psd_save_dir = 'data/raw/l4-mag/'

# Check if these exists
if not os.path.exists(mag_save_dir):
    os.makedirs(mag_save_dir)
if not os.path.exists(psd_save_dir):
    os.makedirs(psd_save_dir)

# Read in the pickle file with times to download data for
with open('data/interim/rbsp-quiet-time-location.pickle',
          'rb') as handle:
    passby_dict = pickle.load(handle) 

for key in list(passby_dict.keys())[0:10]:

    logging.info(f'Downloading data for: {key}.')
    
    for probe in passby_dict[key].keys():
        
        # Get unique days for event, might be more than one
        dates = np.unique([d[0].date() for d in
                           passby_dict[key][probe]['Time']])
        
        for date in dates:
        
            if date > datetime(2019, 7, 16).date():
                continue
                
            # Web directory where l4 psd files are stored
            l4_data_url = ('https://emfisis.physics.uiowa.edu/Flight/RBSP-'
                            + probe[-1].upper() + '/L4/' 
                            + str(date.year) + '/' + str(date.month).zfill(2)
                            + '/' + str(date.day).zfill(2) + '/')
            # Find all files in html directory
            l4_files = get_url_paths(l4_data_url, '.cdf')
                
            # See if there is a sheath corrected file
            e_corrected_filebase = ('rbsp-' + probe[-1].lower()
                                    + '_wna-survey-sheath-corrected-e_emfisis-L4_'
                                    + str(date.year) + str(date.month).zfill(2) 
                                    + str(date.day).zfill(2))
            # l4_psd_filebase = ('rbsp-' + probe[-1].lower()
            #                     + '_wna-survey_emfisis-L4_'
            #                     + str(date.year) + str(date.month).zfill(2) 
            #                     + str(date.day).zfill(2))
            
            # If there isn't a sheath corrected file notify and create a 'no-data' file
            try:
                e_corrected_filepathname = [f for f in l4_files if 
                                            e_corrected_filebase in f][-1]
                e_corrected_filename = e_corrected_filepathname.split('/')[-1]
            except Exception as e:
                logging.warning(f'No e sheath corrected file for {date}.'
                                f' Creating file and continuing with error {e}')
                # If file doesn't exist for density create a proxy file
                proxy_file  = open(psd_save_dir + f'nodata-{date}-{probe}', 'w')
                proxy_file.close()
                continue
                
            # And download it if the file doesn't already exists
            if not os.path.exists(psd_save_dir + e_corrected_filename):
                try:
                    wget.download(e_corrected_filepathname, psd_save_dir, bar=None)
                except Exception as e:
                    logging.warning(f'Unable to download: {e_corrected_filepathname}'
                                    f' with error: {e}.')
                    continue
                
            # # Get the specific L4 PSD file we are looking for
            # l4_psd_filepathname = [f for f in l4_files if l4_psd_filebase in f][-1]
            # l4_psd_filename = l4_psd_filepathname.split('/')[-1]
            
            # # And download it if the file doesn't already exists
            # if not os.path.exists(psd_save_dir + l4_psd_filename):
            #     try:
            #         wget.download(l4_psd_filepathname, psd_save_dir, bar=None)
            #     except Exception as e:
            #         logging.warning(f'Unable to download: {l4_psd_filepathname} with'
            #                         f' error: {e}.')
            #         continue
            
            # Web directory where magnetic field (for gyrofrequency) files are stored
            # mag_file_url = ('https://emfisis.physics.uiowa.edu/Flight/RBSP-'
            #                 + probe[-1].upper() + '/L3/' 
            #                 + str(date.year) + '/' + str(date.month).zfill(2)
            #                 + '/' + str(date.day).zfill(2) + '/')
            mag_file_url = (f'https://cdaweb.gsfc.nasa.gov/pub/data/rbsp/{probe}/'
                            f'l3/emfisis/magnetometer/4sec/gei/{date.year}/')

            mag_filebase = ('rbsp-' + probe[-1].lower() +
                            '_magnetometer_4sec-gei_emfisis-l3_'
                            + str(date.year) + str(date.month).zfill(2) 
                            + str(date.day).zfill(2))
            
            # Find all files in html directory
            mag_files = get_url_paths(mag_file_url, '.cdf')

            # Get the specific file we are looking for
            mag_filepathname = [f for f in mag_files if mag_filebase in f][0]
            mag_filename = mag_filepathname.split('/')[-1]

            # And download it
            # but only if the file doesn't already exists
            if not os.path.exists(mag_save_dir + mag_filename):
                try:
                    wget.download(mag_filepathname, mag_save_dir, bar=None)
                except Exception as e:
                    logging.warning(f'Unable to download: {mag_filepathname} with'
                                    f' error: {e}.')
                    continue

logging.info('All finished.')