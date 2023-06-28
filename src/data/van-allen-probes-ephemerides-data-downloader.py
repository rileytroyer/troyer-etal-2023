"""Script to download ephemerides data from the Van Allen Probes

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
import requests
import wget



# Initiate logging
logging.basicConfig(filename = f'logs/van-allen-probes-ephemerides-download-{datetime.today().date()}.log',
                    encoding='utf-8',
                    format='%(asctime)s %(levelname)-8s %(message)s',
                    level=logging.INFO,
                    datefmt='%Y-%m-%d %H:%M:%S')

####################### End Initializing #######################


####################### Local Functions #######################
def get_url_paths(url, ext='', params={}):
    """ Function to extract file names from https directory
    Gets files in url directory with ext extension
    Does this by parsing the html text from the webpage. I did not
    write this function. Requires libraries requests,
    bs4.BeautifulSoup
    DEPENDENCIES
        bs4.BeautifulSoup, requests
    INPUT
    url
        type: string
        about: url of directory to get files from
    ext
        type: string
        about: extension of the files
    OUTPUT
    parent
        type: list
        about: list of all file pathnames within directory
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

def download_rbsp_files(year:int):
    """Function to download all the magnetic ephemeris data from 
    the Van Allen probes for a list of years and the specified
    probe.
    INPUT
    years - downloads all data for this year
    OUTPUT
    logging information
    """
    
    # Where to save the data to
    save_dir='data/raw/rbsp-magephem/'

    # Create directory if it doesn't exists
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)

    # Loop through all the specified years
    for probe in ['rbspa', 'rbspb']:

        # Look for the files here
        # url = ('https://rbsp-ect.newmexicoconsortium.org/data_pub/' +
        #                 probe + '/MagEphem/definitive/' + year + '/')
        url = (f'https://cdaweb.gsfc.nasa.gov/pub/data/rbsp/{probe}/'
               f'ephemeris/ect-mag-ephem/hdf5/def-1min-t89q/{year}/')
        ext = 'h5' #...with this file type

        logging.info(f'Downloading data from: {url}.')

        # Run function to get all file pathnames 
        files = get_url_paths(url, ext)

        logging.info(f'Response received from: {url}.')

        # Only use T89Q magnetic field model
        files_t89q = np.array([s for s in files if 'T89Q' in s])

        counter = 0

        for file_path in files_t89q: 
            
            # Download the file, sometimes this doesn't work at first
            result = None
            counter = 0
            while result is None:

                # Break out of loop after 10 iterations
                if counter > 10:
                    logging.warning(f'Unable to get data from: {file_path}')
                    break

                try:
                    # Before downloading check if file already exists
                    if not os.path.exists(save_dir + file_path[-42:]):
                        wget.download(file_path, save_dir + file_path[-42:], bar=None)
                    
                    result = True

                except:
                    pass

                counter = counter + 1

            
        logging.info(f'Finished collecting data for: {year} and {probe}.')
####################### End of Local Functions #######################


####################### START OF PROGRAM #######################

if __name__ == '__main__':

    # Download files for specified years and probe
    years = ['2012']#['2012', '2013', '2014', '2015', '2016', '2017', '2018', '2019']

    # Specify where to save files
    save_dir = 'data/raw/rbsp-magephem/'

    # Create the directory if it doesn't already exist
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)

    logging.info('Download starting.')


    try: 

        # Get number of threads for multiprocessing
        num_workers = 8#multiprocessing.cpu_count()
        pool = multiprocessing.get_context("spawn").Pool(processes=num_workers)
        logging.info(f'Multiprocessing pool generated, num_workers = {num_workers}.')

        # Start downloading in multiple processes
        pool.map(download_rbsp_files, years)

        # Make sure to close then join afterwards
        pool.close()
        pool.join()

        logging.info(f'Pool joined. All downloads finished.')

        # for year in years:
        #     download_rbsp_files(year)

    except KeyboardInterrupt():
        logging.info('Program killed.')
        try:
            pool.close()
            pool.join()
        except:
            logging.info('No pools to kill.')