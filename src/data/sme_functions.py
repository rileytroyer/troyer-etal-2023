"""Functions to analyze SME index

@author Riley Troyer
science@rileytroyer.com
"""

from datetime import datetime
import numpy as np
from scipy.ndimage import uniform_filter1d

def sme_read_process(filepath:str, smooth_size:int=6) -> 'np.ndarray, np.ndarray':
    """Function to read in downloaded SME data and smooth it
    using a basic rolling average scipy uniform_filter1D.
    INPUT
    filepath- filepath where the sme datafile is stored.
    smooth_size - how big should the smoothing window be
    OUTPUT
    sme - smoothed sme data.
    sme_dates - timestamp of each sme value
    """
    
    # Read in a file
    sme_data = np.loadtxt(filepath, dtype='str', skiprows=105)
    
    # Parse out datetimes
    sme_dates = np.array([datetime(int(d[0]), int(d[1]), int(d[2]),
                             int(d[3]), int(d[4]), int(d[5]))
                         for d in sme_data])

    # Get the SME data
    sme = sme_data[:, 6].astype(int)

    # Smooth SME
    sme = uniform_filter1d(sme, size=6)
    #sme = savgol_filter(sme, 7, 4)
    
    return sme, sme_dates