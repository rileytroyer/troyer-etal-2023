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
from pathlib import Path
import pickle
from scipy.integrate import simpson
from scipy.interpolate import interp1d
import sys

# Add root to path
path_root = Path(__file__).parents[2]
sys.path.append(str(path_root))

# Function to read in PFISR data
from src.data.van_allen_probe_functions import read_process_rbsp_data


# Initiate logging
logging.basicConfig(filename = f'logs/van-allen-probe-injection-data-compiler-{datetime.today().date()}.log',
                    encoding='utf-8',
                    format='%(asctime)s %(levelname)-8s %(message)s',
                    level=logging.INFO,
                    datefmt='%Y-%m-%d %H:%M:%S')

####################### End Initializing #######################


####################### Local Functions #######################
def save_dict_to_h5(event:str, dic:dict, path:str, h5_file:h5py.File):
    """Convert dictionary into h5py file.
    INPUT
    event - isotime of event that dictionary has data for.
    dic - dictionary of data for event
    path - base directory within h5 file, usually us '/'.
    h5_file - the h5 file to write data to.
    OUTPUT
    Writes to h5 file. Otherwise raises error.
    """
    
    # Create new group within file for event
    group_path = path + event + '/'
    
    # Loop through each key in dictionary
    for key, item in dic.items():
        
        # Key type of dictionary item to make sure it is writable to h5
        if isinstance(item, (np.ndarray, np.int64, np.float64, str, bytes)):
            h5_file[group_path + key] = item
        
        else:
            raise ValueError('Cannot save %s type'%type(item))

def filter_write_to_dict(i, threshold = 10**-7):
    """Function to filter chorus data by psd threshold,
    LBC/UBC.
    INPUT
    i
        type: integer
        about: loop location
    threshold = 10**-7
        type: integer
        about: magnetic psd threshold in nT^2/Hz from Hartley et al. 2019
    OUTPUT
    none, it just write to dictionary that is already created.
    """
    
    # Set probe
    chorus_delay_dict['probe'].append(probe)
    
    # Set time
    time = ut_time[i]
    chorus_delay_dict['ut'].append(time)
    
    # Add location data to dictionary
    chorus_delay_dict['mlt'].append(mlt[i])
    chorus_delay_dict['l'].append(l[i])
    chorus_delay_dict['mlat'].append(mlat[i])
    
    # Determine delay in seconds
    chorus_delay_dict['delay'].append((time - event).total_seconds())

    # Convert time to float
    int_time = time.timestamp()

    # Get chorus frequencies
    # Calculate 0.1, 0.5 and 1 gyrofrequency
    fce = b_mag_func(int_time)*28
    fce_half = fce/2
    fce_tenth = fce/10

    # Slice by time and select just lower band chorus waves
    lbc_freq = freq[(freq > fce_tenth) & (freq < fce_half)]
    b_lbc_psd = b_power[(freq > fce_tenth) & (freq < fce_half), i]
    e_lbc_psd = e_power[(freq > fce_tenth) & (freq < fce_half), i]

    try:
        # Only continue if max magnetic chorus is > threshold
        if np.max(b_lbc_psd) >= threshold:

            # Integrate psd over LBC frequencies
            integrated_b_lbc_psd = simpson(y=b_lbc_psd[b_lbc_psd != -1e31],
                                           x=lbc_freq[b_lbc_psd != -1e31])
            integrated_e_lbc_psd = simpson(y=e_lbc_psd[e_lbc_psd != -1e31],
                                           x=lbc_freq[e_lbc_psd != -1e31])

            # Get max value
            chorus_delay_dict['b_lbc_max'].append(np.max(b_lbc_psd[b_lbc_psd 
                                                                   != -1e31]))
            chorus_delay_dict['e_lbc_max'].append(np.max(e_lbc_psd[e_lbc_psd 
                                                                   != -1e31]))

        else:
            integrated_b_lbc_psd = np.nan
            chorus_delay_dict['b_lbc_max'].append(np.nan)
            integrated_e_lbc_psd = np.nan
            chorus_delay_dict['e_lbc_max'].append(np.nan)

    except:
        integrated_b_lbc_psd = np.nan
        chorus_delay_dict['b_lbc_max'].append(np.nan)
        integrated_e_lbc_psd = np.nan
        chorus_delay_dict['e_lbc_max'].append(np.nan)


    # Add to list
    chorus_delay_dict['b_lbc'].append(integrated_b_lbc_psd)
    chorus_delay_dict['e_lbc'].append(integrated_e_lbc_psd)

    # Slice by time and select just upper band chorus waves
    ubc_freq = freq[(freq > fce_half) & (freq < fce)]
    b_ubc_psd = b_power[(freq > fce_half) & (freq < fce), i]
    e_ubc_psd = e_power[(freq > fce_half) & (freq < fce), i]

    try:
        # Only continue if max chorus is > threshold
        if np.max(b_ubc_psd) >= threshold:

            # Integrate psd over UBC frequencies
            integrated_b_ubc_psd = simpson(y=b_ubc_psd[b_ubc_psd != -1e31],
                                           x=ubc_freq[b_ubc_psd != -1e31])
            integrated_e_ubc_psd = simpson(y=e_ubc_psd[e_ubc_psd != -1e31],
                                           x=ubc_freq[e_ubc_psd != -1e31])

            # Get max value
            chorus_delay_dict['b_ubc_max'].append(np.max(b_ubc_psd[b_ubc_psd 
                                                                   != -1e31]))
            chorus_delay_dict['e_ubc_max'].append(np.max(e_ubc_psd[e_ubc_psd 
                                                                   != -1e31]))

        else:
            integrated_b_ubc_psd = np.nan
            chorus_delay_dict['b_ubc_max'].append(np.nan)
            integrated_e_ubc_psd = np.nan
            chorus_delay_dict['e_ubc_max'].append(np.nan)

    except:
        integrated_b_ubc_psd = np.nan
        chorus_delay_dict['b_ubc_max'].append(np.nan)
        integrated_e_ubc_psd = np.nan
        chorus_delay_dict['e_ubc_max'].append(np.nan)

    # Add to list
    chorus_delay_dict['b_ubc'].append(integrated_b_ubc_psd)
    chorus_delay_dict['e_ubc'].append(integrated_e_ubc_psd)
####################### End of Local Functions #######################


####################### START OF PROGRAM #######################

mag_save_dir = 'data/raw/mag-waveform/'
psd_save_dir = 'data/raw/l4-mag/'

# Threshold, more than this is chorus
threshold = 10**-7 # from Hartley et al. 2019

# Create dictionary to store data
chorus_delay_dict = {'delay' : [],
                     'b_ubc' : [],
                     'e_ubc' : [],
                     'b_ubc_max' : [],
                     'e_ubc_max' : [],
                     'b_lbc' : [],
                     'e_lbc' : [],
                     'b_lbc_max' : [],
                     'e_lbc_max' : [],
                     'mlt' : [],
                     'l' : [],
                     'mlat' : [],
                     'probe' : [],
                     'ut':[]}

# Create H5 file to store data in
h5_data_filename = 'data/processed/chorus-delay-data.h5'

# Read in the pickle file
with open('data/interim/rbsp-quiet-time-location.pickle',
          'rb') as handle:
    passby_dict = pickle.load(handle)

# Loop  through each quiet time event

logging.info(f'Starting program. {len(passby_dict.keys())} events to process.')

offset = 0
offset_end = 45
for i, event in enumerate(list(passby_dict.keys())[0:]):
    
    logging.info(f'Starting event {event}.')
    
    i = i + offset
    
    # Loop through each probe in event
    for probe in passby_dict[event].keys():
        
        #logging.info(f'Starting probe {probe}.')
        
        # Get the unique dates in event
        dates = np.unique([d[0].date() for d in
                           passby_dict[event][probe]['Time']])
        
        if dates[0] > datetime(2019, 7, 16).date():
            logging.warning(f'Date {date} after 2019-07-16.')
            continue
            
        for j, date in enumerate(dates):
            
            if date > datetime(2019, 7, 16).date():
                logging.warning(f'Date {date} after 2019-07-16.')
                continue
                
            if j==0:
                # Read in data files
                try:
                    (ut_time, freq,
                     b_power, e_power, density,
                     ut_time_b_mag, b_mag,
                     l, mlt, mlat) = read_process_rbsp_data(probe, date, psd_save_dir, mag_save_dir)
                except Exception as e:
                    logging.warning(f'Unable to read in rbsp data for {probe} and {date}.'
                                    f' Returned error {e}.')
                    ut_time = [np.nan]
                    continue
                
            else:
                try:
                    # Read in data files
                    (ut_time_tmp, freq_tmp,
                     b_power_tmp, e_power_tmp, density_tmp,
                     ut_time_b_mag_tmp, b_mag_tmp,
                     l_tmp, mlt_tmp, mlat_tmp) = read_process_rbsp_data(probe, date, psd_save_dir, mag_save_dir)
                    
                    # Check if density is nan, if so skip
                    if np.nan in ut_time_tmp:
                        logging.warning(f'NaN in time for {probe} and {date}.')
                        continue
                except Exception as e:
                    logging.warning(f'Unable to read in rbsp data for {probe} and {date}.'
                                    f' Returned errror {e}.')
                    ut_time_tmp = [np.nan]
                    continue
                    
                # If no density with first date don't append
                if np.nan in ut_time:
                    
                    logging.warning(f'NaN in density for date prior to {date} for {probe}.')
                    
                    # Append to first data
                    ut_time = ut_time_tmp
                    b_power = b_power_tmp
                    e_power = e_power_tmp
                    ut_time_b_mag = ut_time_b_mag_tmp
                    b_mag = b_mag_tmp
                    density = density_tmp
                    l = l_tmp
                    mlt = mlt_tmp
                    mlat = mlat_tmp
                
                else:
                    # Append to first data
                    ut_time = np.concatenate((ut_time, ut_time_tmp), axis=0)
                    b_power = np.concatenate((b_power, b_power_tmp), axis=1)
                    e_power = np.concatenate((e_power, e_power_tmp), axis=1)
                    ut_time_b_mag = np.concatenate((ut_time_b_mag, ut_time_b_mag_tmp),
                                                     axis=0)
                    b_mag = np.concatenate((b_mag, b_mag_tmp), axis=0)
                    density = np.concatenate((density, density_tmp), axis=0)
                    
                    l = np.concatenate((l, l_tmp), axis=0)
                    mlt = np.concatenate((mlt, mlt_tmp), axis=0)
                    mlat = np.concatenate((mlat, mlat_tmp), axis=0)
                    
        # Check if there is any data for event
        if np.nan in ut_time:
            logging.warning(f'No data for entire event {event} for {probe}.')
            del ut_time
            
            # Not entirely sure this is needed, but I don't think it hurts
            try:
                del b_power, e_power, b_mag, ut_time_b_mag, density, l, mlt, mlat
            except:
                logging.warning('No psd variables to delete')
            continue
                
        # Select only data within desired times
        start_time = sorted(passby_dict[event][probe]['Time'])[0][0].replace(tzinfo=None)
        end_time = sorted(passby_dict[event][probe]['Time'])[-1][0].replace(tzinfo=None)
        
        # Select data between specified times for data in psd file
        psd_selector = (ut_time >= start_time) & (ut_time <= end_time)
        b_power = b_power[:, psd_selector]
        e_power = e_power[:, psd_selector]
        density = density[psd_selector]
        l = l[psd_selector]
        mlt = mlt[psd_selector]
        mlat = mlat[psd_selector]
        ut_time = ut_time[psd_selector]
        
        # Select data between specified times for B-field magnitude file
        b_field_selector = (ut_time_b_mag >= start_time) & (ut_time_b_mag <= end_time)
        b_mag = b_mag[b_field_selector]
        ut_time_b_mag = ut_time_b_mag[b_field_selector]
        
        # If there isn't enough data, skip
        if len(ut_time_b_mag) < 1:
            
            logging.warning(f'Not enough b_mag data for {probe} and {event}.')
            
            del ut_time, b_power, e_power, b_mag, ut_time_b_mag, density, l, mlt, mlat
            continue
        
        # If there isn't enough density data, skip
        if len(density) < 1:
            
            logging.warning(f'Not enough density data for {probe} and {event}.')
            
            del ut_time, b_power, e_power, b_mag, ut_time_b_mag, density, l, mlt, mlat
            continue
        
        # Convert b_field and density times into floats
        int_time_b_mag = np.array([t.timestamp() for t in ut_time_b_mag])
        
        # Create a linearly interpolated model of b_field
        # Time needs to be in float format
        b_mag_func = interp1d(int_time_b_mag, b_mag,
                                fill_value='extrapolate')
        
        for k, time in enumerate(ut_time):
            
            # Check if density is low enough
            # Based on Li et al. 2010
            den_check = density[k]
            l_check = l[k]
            
            # Smaller of 10(6.6/L)**4 or 50 cm^-3
            small_den = 10*(6.6/l_check)**4
            if small_den > 50:
                small_den = 50
            
            # If density isn't low enough skip
            if den_check > small_den:
                #logging.warning(f'Density too large for {probe} and {event}.')
                continue
            
            # Filter data based on threshold and frequency
            filter_write_to_dict(k, threshold = 10**-7)
            
        # Clear variables
        del ut_time, b_power, e_power, b_mag, ut_time_b_mag, density, l, mlt, mlat
        
        #logging.info(f'Finished processing {event} for {probe}')
        
    # Write dictionary to h5 file
    
    # First change time to iso format string
    chorus_delay_dict['ut'] = np.array([t.isoformat() + 'Z' for t in 
                                        chorus_delay_dict['ut']]).astype('S27')
    
    # Turn dictionary items into arrays
    for key, item in chorus_delay_dict.items():
        chorus_delay_dict[key] = np.array(item)
        
    # Convert probe to bytes type
    chorus_delay_dict['probe'] = chorus_delay_dict['probe'].astype('S5')
    
    with h5py.File(h5_data_filename, 'a') as h5_file:
        try:
            save_dict_to_h5(event.isoformat() + 'Z', chorus_delay_dict, '/', h5_file)
        except Exception as e:
            logging.warning(f'Unable to write event for {date} into h5 file.'
                            f' Returned error {e}.')
        
    # Clear data dictionary
    for key in chorus_delay_dict:
        chorus_delay_dict[key] = []
    
    logging.info(f'Finished processing {event}.')

# Lastly add information to h5 file
with h5py.File(h5_data_filename, 'a') as h5_file:
    h5_file.attrs['about'] = ('Magnetic and electric chorus data from RBSP EMFISIS. '
                              'Organized by event -> individual measurements. '
                              'Times are in ut datasets and stored in an '
                              'ISO format byte string. '
                              'To convert times to datetime run: '
                              'datetime.datetime.fromisoformat(ISOTIME.decode("utf-8")')


# # Write the dictionary with conjunction times to a pickle file
# #...if we need it again we don't have to calculate it all out
# with open('data/processed/chorus-delay-data.pickle',
#                   'wb') as handle:
#     pickle.dump(chorus_delay_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)

logging.info(f'All finished. H5 file saved at: {h5_data_filename}')