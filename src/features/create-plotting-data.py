""" Script to compile Van Allen Probe EMFISIS data for the times during injections.

@author Riley Troyer
science@rileytroyer.com
"""

####################### Initialize Program #######################
# Libraries
from datetime import datetime
import h5py
import logging
import numpy as np


# Initiate logging
logging.basicConfig(filename = f'logs/create-plotting-data-{datetime.today().date()}.log',
                    encoding='utf-8',
                    format='%(asctime)s %(levelname)-8s %(message)s',
                    level=logging.INFO,
                    datefmt='%Y-%m-%d %H:%M:%S')

####################### End Initializing #######################


####################### Local Functions #######################

####################### End of Local Functions #######################


####################### START OF PROGRAM #######################
# Read in the h5 data with the processed and curated data
data_file = h5py.File('data/processed/chorus-delay-data.h5')

psd_type = 'integrated'

if psd_type == 'integrated':
    extension = ''
if psd_type == 'max':
    extension = '_max'

# Initialize lists to store data in
mlt, l, mlat = [], [], []
delay = []
ubc_b, lbc_b = [], []
ubc_e, lbc_e = [], []

# Loop through each group (key) in h5 file and combine all data
for n, group in enumerate(data_file):

    # Read in data for a single event
    event_data = data_file[group]

    # Extend lists
    mlt.extend(data_file[group]['mlt'])
    l.extend(data_file[group]['l'])
    mlat.extend(data_file[group]['mlat'])
    delay.extend(data_file[group]['delay'])
    ubc_b.extend(data_file[group]['b_ubc' + extension])
    lbc_b.extend(data_file[group]['b_lbc' + extension])
    ubc_e.extend(data_file[group]['e_ubc' + extension])
    lbc_e.extend(data_file[group]['e_lbc' + extension])

    if n%100 == 0:
        logging.info(f'Finished with {n} of {len(data_file)} events.')

# Convert to arrays and change nan to 0
ubc_b = np.nan_to_num(ubc_b)
lbc_b = np.nan_to_num(lbc_b)
ubc_e = np.nan_to_num(ubc_e)
lbc_e = np.nan_to_num(lbc_e)

# Create full chorus arrays
if psd_type == 'integrated':
    chorus_b = ubc_b + lbc_b
    chorus_e = ubc_e + lbc_e
if psd_type == 'max':
    chorus_b = np.maximum(ubc_b, lbc_b)
    chorus_e = np.maximum(ubc_e, lbc_e)

delay = np.nan_to_num(delay)
mlt = np.nan_to_num(mlt)
l = np.nan_to_num(l)
mlat = np.nan_to_num(mlat)

# Write to h5 file
with h5py.File(f'data/processed/analysis-data-{psd_type}.h5', 'w') as h5f:

    # Create datasets for each variable
    ubc_b_ds = h5f.create_dataset('ubc_b', shape=ubc_b.shape,
                                  dtype=float, data=ubc_b)
    lbc_b_ds = h5f.create_dataset('lbc_b', shape=lbc_b.shape,
                                  dtype=float, data=lbc_b)
    chorus_b_ds = h5f.create_dataset('chorus_b', shape=chorus_b.shape,
                                  dtype=float, data=chorus_b)
    ubc_e_ds = h5f.create_dataset('ubc_e', shape=ubc_e.shape,
                                  dtype=float, data=ubc_e)
    lbc_e_ds = h5f.create_dataset('lbc_e', shape=lbc_e.shape,
                                  dtype=float, data=lbc_e)
    chorus_e_ds = h5f.create_dataset('chorus_e', shape=chorus_e.shape,
                                  dtype=float, data=chorus_e)

    delay_ds = h5f.create_dataset('delay', shape=delay.shape,
                                  dtype=float, data=delay)
    mlt_ds = h5f.create_dataset('mlt', shape=mlt.shape,
                                  dtype=float, data=mlt)
    l_ds = h5f.create_dataset('l', shape=l.shape,
                                  dtype=float, data=l)
    mlat_ds = h5f.create_dataset('mlat', shape=mlat.shape,
                                  dtype=float, data=mlat)    

    # Create some global informational attributes
    h5f.attrs['PSD Selection Type'] = psd_type

    if psd_type == 'integrated':
        units = ('B: nT^2, E: mV^2/m^2, delay: s, mlt: hr, mlat: deg')
    if psd_type == 'max':
        units = ('B: nT^2/Hz, E: mV^2/m^2/Hz, delay: s, mlt: hr, mlat: deg')

    h5f.attrs['Units'] = units

logging.info(f'Finished. Data stored at: data/processed/analysis-data-{psd_type}.h5')