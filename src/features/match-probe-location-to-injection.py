""" Script to look at all Van Allen Probe locations and only find those that are during injection periods as estimated from SME.

@author Riley Troyer
science@rileytroyer.com
"""

####################### Initialize Program #######################
# Libraries
# Libraries for notebook
from datetime import datetime
from dateutil import parser
import h5py
import numpy as np
import os
import pandas as pd
import pickle
import pytz



# Initiate logging
logging.basicConfig(filename = f'logs/match-probe-location-to-injection-{datetime.today().date()}.log',
                    encoding='utf-8',
                    format='%(asctime)s %(levelname)-8s %(message)s',
                    level=logging.INFO,
                    datefmt='%Y-%m-%d %H:%M:%S')

####################### End Initializing #######################


####################### Local Functions #######################

####################### End of Local Functions #######################


####################### START OF PROGRAM #######################

# Filtering parameters
small_l = 3
large_mlat = 30 #degrees
start_mlt = 0
end_mlt = 24
bad_values = -1e-31

# Read in the SME data file
quiet_times = pd.read_csv('data/interim/sme-injections-quiet-times.txt', delimiter=',',
                          names=['Quiet Start', 'Quiet End', 'Injection Start',
                                 'Injection Length', 'Quiet Length'])

# Get a list of RBSP magnetic ephemerides
# Get a list of all files in directory
footpoint_dir = 'data/raw/rbsp-magephem/'
footpoint_files = os.listdir(footpoint_dir)

# I like to filter to only the file types I want 
#...in case there are others
footpoint_files = [f for f in footpoint_files if f.endswith('.h5')]

# Dictionary to store locations in
rbsp_matched_quiet_times = {}

start_i = 0
# Loop through each quiet and get RBSP probe positions
for df_index in quiet_times.index:
    
    # Figure out which files we will need to read in
    # May need more than one if quiet period is in 2 UTC days
    start_time = parser.isoparse(quiet_time.loc[df_index, 'Quiet Start'])
    end_time = parser.isoparse(quiet_time.loc[df_index, 'Quiet End'])

    file_dates = pd.date_range(start_time.date(), end_time.date(),
                               freq='d').strftime('%Y%m%d').tolist()

    # Select the magnetic ephemerides files
    rbsp_filenames = [f for f in footpoint_files if any(d in f for d in file_dates)]
    
    if len(rbsp_filenames) == 0:
        continue
    
    # Set dictionary
    rbsp_matched_quiet_times[start_time] = {'rbspa' : {},
                                               'rbspb' : {}}

    for probe in ['rbspa', 'rbspb']:
        
        # Lists to store location in
        mlat = []
        l = []
        mlt = []
        isotime = []
        footpoint = []
        rgeo = []
        
        for rbsp_filename in [f for f in rbsp_filenames if probe in f]:
            
            # Read in the file
            file = h5py.File(footpoint_dir + rbsp_filename, 'r')

            # Get times and add to date to produce a datetime
            isotime.extend(list(file['IsoTime']))

            # Read data into arrays
            mlat.extend(list(file['CDMAG_MLAT']))
            l.extend(list(file['L'][:, -1]))
            mlt.extend(list(file['CDMAG_MLT']))
            
            # Add footprint data too
            footpoint.extend(list(file['Pfn_geod_LatLon']))
            
            # Add geographic location vector
            rgeo.extend(np.array(file['Rgeo']))
            
        # Convert time to datetime and array
        # Need an extra function to account for some 
        # seconds being 60
        def parse_func(t):
            if t[-3:-1] == b'60':
                t = t[0:-3] + b'59Z'
            return parser.isoparse(t)
        
        isotime = np.array([parse_func(t) for t in isotime])
        mlat = np.array(mlat)
        l = np.array(l)
        mlt = np.array(mlt)
        footpoint = np.array(footpoint)
        rgeo = np.array(rgeo)
        
        # Only get times during the quiet period
        # And when probe was in desired location
        selected_i = np.argwhere((isotime > start_time.replace(tzinfo=pytz.UTC)) 
                                 & (isotime < end_time.replace(tzinfo=pytz.UTC))
                                 & ((mlt < end_mlt) | (mlt > start_mlt))
                                 & (np.abs(mlat) < large_mlat)
                                 & (l > small_l)
                                 & (l != bad_values))
        
        # Skip if no times that fit
        if len(selected_i) == 0:
            del rbsp_matched_quiet_times[start_time][probe]
            continue
        
        # Otherwise write to dictionary
        rbsp_matched_quiet_times[start_time][probe] = {'Time' : isotime[selected_i],
                                                          'MLT' : mlt[selected_i],
                                                          'L' : l[selected_i],
                                                          'MLat' : mlat[selected_i],
                                                          'Footpoint' : footpoint[selected_i],
                                                          'Rgeo' : rgeo[selected_i]}

# Loop through dictionary and remove entries with out any data
keys = list(rbsp_matched_quiet_times.keys())
for key in keys:
    
    if rbsp_matched_quiet_times[key] == {}:
        del rbsp_matched_quiet_times[key]

# Write the dictionary with conjunction times to a pickle file
#...if we need it again we don't have to calculate it all out
with open('data/interim/rbsp-quiet-time-location.pickle',
                  'wb') as handle:
    pickle.dump(rbsp_matched_quiet_times, handle, protocol=pickle.HIGHEST_PROTOCOL)

