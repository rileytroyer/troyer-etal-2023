""" Script to find injection periods and quiet time following injections
using the SME index data.

@author Riley Troyer
science@rileytroyer.com
"""

####################### Initialize Program #######################
# Libraries
from datetime import datetime
import logging
import numpy as np
import os
from pathlib import Path
from scipy.ndimage import uniform_filter1d
from scipy.signal import savgol_filter
import sys

# Add root to path
path_root = Path(__file__).parents[2]
sys.path.append(str(path_root))

# Function to read in PFISR data
from src.data.sme_functions import sme_read_process



# Initiate logging
logging.basicConfig(filename = f'logs/find-injections-with-sme-{datetime.today().date()}.log',
                    encoding='utf-8',
                    format='%(asctime)s %(levelname)-8s %(message)s',
                    level=logging.INFO,
                    datefmt='%Y-%m-%d %H:%M:%S')

####################### End Initializing #######################


####################### Local Functions #######################
def sme_find_quiet_times(sme, sme_dates, quiet_threshold=150,
                         high_threshold=250):
    """Function to find quiet times with no injection as seen in the
    SME data. This is defined as when SME is less than threshold
    INPUT
    sme
        type: array of ints
        about: sme values
    sme_dates
        type: array of datetimes
        about: datetimes associated with each sme value
    quiet_threshold=150:
        type: int
        about: below this is considered quiet time
    high_threshold=250:
        type: int
        about: above this is considered most likely injection
    OUTPUT
    quiet_times
        type: list of lists
        about: start and stop times of each quiet time period

    """
    
    # Find where sme is below threshold
    low_sme_i = np.argwhere(sme < quiet_threshold)[:, 0]
    
    # Group by connected times
    quiet_times_i = [[low_sme_i[0]]]
    
    for i in range(1, len(low_sme_i)):
        if low_sme_i[i-1]+1 == low_sme_i[i]:
            quiet_times_i[-1].append(low_sme_i[i])
            
        else:
            quiet_times_i.append([low_sme_i[i]])
            
    # Loop through each group and get start and stop times of quiet periods
    quiet_times = []
    prev_end = 0
    for group in quiet_times_i:

        # Remove if period is less than specified length
        if len(group) < 10:
            continue
            
        # Remove if no higher SME prior to last period end
        if np.max(sme[prev_end:group[-1]]) < high_threshold:
            continue
        
        # Try to find where the start of the injection happens, where the spike above high threshold starts
        # First get the index of the peak value in index
        peak_index = np.argmax(sme[prev_end:group[-1]])

        # Now get all index values that are less than threshold
        non_active_times = np.argwhere(sme[prev_end:group[-1]] < high_threshold)

        # Of these get the largest that is less than the peak index
        try:
            injection_start_index = sorted(non_active_times[non_active_times < peak_index])[-1]
        except:
            injection_start_index = prev_end
        injection_start_time = sme_dates[prev_end:group[-1]][injection_start_index]

        # Get the approximate injection length, or at least between end of last and start of current
        injection_len = (sme_dates[group[0]] - injection_start_time).total_seconds()

        # Get the length of a quiet period
        quiet_period_len = (sme_dates[group[-1]] - sme_dates[group[0]]).total_seconds()
        
        # Reset previous end point
        prev_end = group[-1]
        
        # Get first and last index and pull datetime from this
        quiet_times.append([sme_dates[group[0]], sme_dates[group[-1]],
                            injection_start_time, injection_len, quiet_period_len])
    
    return quiet_times
####################### End of Local Functions #######################


####################### START OF PROGRAM ####################### 

# Files with SME data
sme_dir = 'data/raw/sme/'
sme_files = sorted(os.listdir(sme_dir))

all_quiet_times = []

# Loop through each file
for file in sme_files:

    logging.info(f'Finding injections for {file}')
    # Read in and smooth the SME data
    sme_smooth, sme_dates = sme_read_process(sme_dir + file)

    # Find quiet times
    quiet_times = sme_find_quiet_times(sme_smooth, sme_dates,
                                       quiet_threshold=150,
                                       high_threshold=250)
    
    all_quiet_times.extend(quiet_times)
    
# Convert to array
all_quiet_times = np.array(all_quiet_times)

logging.info('Finished')

# Save times as a text file
np.savetxt('data/interim/sme-injections-quiet-times.txt',
           all_quiet_times.astype(str),
           fmt='%s', delimiter=',')

logging.info('Wrote data to: data/interim/sme-injections-quiet-times.txt')