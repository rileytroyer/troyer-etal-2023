""" Functions for working with the EMFISIS Van Allen Probe data.

@author Riley Troyer
science@rileytroyer.com
"""

# Libraries
import cdflib
from datetime import datetime
import numpy as np
import os
from scipy.ndimage.filters import uniform_filter1d
from scipy.signal import savgol_filter


def read_rbsp_sheath_corrected_psd(file_name:str) -> 'np.ndarray x 8, str':
    """Function to read in a Van Allen EMFISIS sheath corrected 
    datafile and convert the data to numpy arrays
    INPUT
    file_name - name of file, needs to be a .cdf file
    OUTPUT
    ut_time -  datetimes for each data point
    freq - requencies for each data point
    b_power - total magnetic field power values for each data point
    e_power - total sheath corrected electric field power values for each data point
    density - electron density measurements
    l, mlt, mlat - locations of instruments
    instrument - which instrument data is from
    """
    
    cdf_file = cdflib.CDF(file_name)
    
    # Get the time and frequency data
    freq = cdf_file.varget('WFR_frequencies')
    time = cdf_file.varget('Epoch')
    #...convert time to datetime format
    ut_time = cdflib.cdfepoch.to_datetime(time).astype(datetime)
    
    # Get the total power data
    b_power = np.transpose(cdf_file.varget('bsum'))
    e_power = np.transpose(cdf_file.varget('esum'))
    
    # Electric field has noise at 1781hz and 3555hz
    e_power[(freq==3555)|(freq==1781), :] = -1e31
    
    # Get density data
    density = cdf_file.varget('density')
    
    # Get L, MLT, and MLAT
    l = cdf_file.varget('l')
    mlt = cdf_file.varget('mlt')
    mlat = cdf_file.varget('maglat')
    
    # # Get additional information, like instrument
    # instrument = cdf_file.globalattsget()['Source_name']
    # instrument = instrument['Source_name'][0][0][0:5]
    
    return ut_time, freq, b_power, e_power, density, l, mlt, mlat#, instrument

def read_rbsp_emfisis_b_field(file_name:str) -> "np.ndarray, np.ndarray":
    """Function to read in a Van Allen L3 EMFISIS datafile and
    retrieve the magnetic b_field strength from it.
    INPUT
    file_name- name of file, needs to be a .cdf file
    OUTPUT
    ut_time - datetimes for each data point
    b_field_mag - magnetic b_field magnitude for each data point
    """
    
    cdf_file = cdflib.CDF(file_name)
    
    # Get the time and magnetic b_field data
    time = cdf_file.varget('Epoch')
    #...convert time to datetime format
    ut_time = cdflib.cdfepoch.to_datetime(time).astype(datetime)

    
    # Get the magnetic b_field magnitude
    b_field_mag = cdf_file.varget('Magnitude')
    
    return ut_time, b_field_mag

def read_process_rbsp_data(probe:str, date:datetime,
                           psd_save_dir:str,
                           mag_save_dir:str) -> 'np.ndarray x 10':
    """ Function to read and return smoothed b-field and wave power data
    from EMFISIS instruments. Also returns associated time and location of spacecraft
    INPUT
    probe - which probe to get data for rbspa or rbspb
    date - date to get data for
    psd_save_dir - where are wave power data files stored
    mag_save_dir - where are magnetic field power data files stored
    OUTPUT
    ut_time - times of measurements
    freq - frequency bins of psd
    b_power - magnetic field power measurements
    e_power - electric field power measurments
    density - electron density measurements
    ut_time_b_mag - times associated with magnetic field measurements
    b_mag - magnetic field measurements
    l, mlt, mlat - locations of instruments
    """ 
    
    # Check if a no data file for density exists    
    if os.path.exists(psd_save_dir + f'nodata-{date}-{probe}'):
        raise Exception(f'File for {date} and {probe} does not exist')
    
    else: 
        # Get sheath corrected e field psd data
        e_corrected_filebase = ('rbsp-' + probe[-1].lower()
                                + '_wna-survey-sheath-corrected-e_emfisis-L4_'
                                + str(date.year) + str(date.month).zfill(2) 
                                + str(date.day).zfill(2))

        # All of the density files
        psd_files = os.listdir(psd_save_dir)

        # Find just the density file we need
        psd_filename = [f for f in psd_files if e_corrected_filebase in f][0]

        # Read in psd data
        (ut_time, freq,
         b_power, e_power,
         density, l, mlt, mlat) = read_rbsp_sheath_corrected_psd(psd_save_dir + psd_filename)

        if len(density) < 1:
            raise Exception(f'Not enough density data for {date} and {probe}.')
    
    # Read in b-field mag (for gyrofrequency) files
    mag_files = os.listdir(mag_save_dir)

    # This is the file we want
    mag_filebase = ('rbsp-' + probe[-1].lower() 
                    + '_magnetometer_4sec-gei_emfisis-L3_'
                    + str(date.year) + str(date.month).zfill(2) 
                    + str(date.day).zfill(2))

    # Get the specific filename
    mag_filename = [f for f in mag_files if mag_filebase in f][0]

    # Smooth power over 6 min
    b_power = uniform_filter1d(b_power, size=6, axis=1)
    e_power = uniform_filter1d(e_power, size=6, axis=1)
    #power = savgol_filter(power, window_length=7, polyorder=3, axis=1)

    # Read in the B-b_field data
    (ut_time_b_mag,
     b_mag) = read_rbsp_emfisis_b_field(mag_save_dir + mag_filename)
    
    return ut_time, freq, b_power, e_power, density, ut_time_b_mag, b_mag, l, mlt, mlat