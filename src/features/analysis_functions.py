""" Functions to help with analysis of post injection chorus with EMFISIS.

@author Riley Troyer
science@rileytroyer.com
"""

# Libraries
from datetime import datetime
from matplotlib import cm
from matplotlib import pyplot as plt
from matplotlib import colors as mcolors
from matplotlib import ticker
import numpy as np
import pickle
from scipy.stats import pearsonr
import scipy.stats as stats
import statsmodels.graphics.gofplots as sm
from sklearn.utils import resample

def log_log_transformation(distribution:np.ndarray) -> 'np.ndarray, float':

    # Log transform the distribution
    distribution = np.log10(distribution)

    # Scale to ranges to be above zero and do another log transform
    scale1 = abs(np.min(distribution)) + 1
    distribution = np.log10(distribution + scale1)

    return distribution, scale1

def mean_of_distribution_transform(distribution:np.ndarray) -> 'float, float, float':
    """Transform distribution to something that is normal, get mean, and return
    what mean would be in original.
    INPUT
    distribution - original distribution
    OUTPUT
    value - original value that mean of normal distribution corresponds to (most likely point)
    error_low, error_high - standard error of the mean in original distribution.
    """

    # Transform the distribution
    distribution, scale1 = log_log_transformation(distribution)

    # test comment
    
    # Perform fit
    params = stats.norm.fit(distribution)

    # Get detransformed value and errors
    def detransformed(value):
        return 10**(10**(value) - scale1)

    value =detransformed(params[0])
    # error_low = detransformed(params[0] - np.std(distribution)/np.sqrt(len(distribution)))
    # error_high = detransformed(params[0] + np.std(distribution)/np.sqrt(len(distribution)))
    error_low = detransformed(params[0] - params[1]/np.sqrt(len(distribution)))
    error_high = detransformed(params[0] + params[1]/np.sqrt(len(distribution)))

    return value, error_low, error_high

def average_bins(delay:np.ndarray, chorus:np.ndarray, method:str='peak',
                 bin_size:int=10, cutoff:float=10**-7):
    """Function to take large arrays of delays and chorus integration
    and get averages over delay bins.
    """
    
    # Average altitudes up to 5 hour after
    delay_bins = np.arange(0, 5*60*60, bin_size*60)
    chorus_bins = np.zeros(len(delay_bins))
    chorus_bins_q1 = np.zeros(len(delay_bins))
    chorus_bins_q3 = np.zeros(len(delay_bins))

    statistics = np.zeros(len(delay_bins))

    def avg_func(data, method='median'):

        if method=='median':
            return np.nanmedian(data)
        if method=='mean':
            return np.nanmean(data)
        if method=='gmean':
            return stats.gmean(data)
        if method=='peak':
            value, error_low, error_high = mean_of_distribution_transform(data)
            return value, error_low, error_high

    # Loop through each bin and average
    for n, delay_bin in enumerate(delay_bins):

        # For last bin average all the rest
        if n == len(delay_bins)-1:
            chorus_bin = chorus[delay >= delay_bin]

        else:
            chorus_bin = chorus[(delay >= delay_bin) 
                                & (delay < delay_bin + bin_size*60)]
            
        # Remove data that wouldn't be considered chorus waves
        #cutoff = 10**-8.5
        #cutoff = 10**-7
        #chorus_bin = chorus_bin[chorus_bin > cutoff]
        # Remove the lowest half of data
        #chorus_bin = chorus_bin[chorus_bin > np.nanmedian(chorus_bin)]
        
        # If no data for bin make undefined
        if len(chorus_bin) < 2:
            chorus_bins[n] = np.nan
            chorus_bins_q1[n] = np.nan
            chorus_bins_q3[n] = np.nan
            continue
        
        # How much data in bin
        statistics[n] = len(chorus_bin[np.isfinite(chorus_bin)])
        
        # Set values where statistics indicate less than 3 events to nan
        min_stat = ((bin_size*60)/6) * 0
        if statistics[n] < min_stat:
            chorus_bins[n] = np.nan
            chorus_bins_q1[n] = np.nan
            chorus_bins_q3[n] = np.nan
            continue

        chorus_bins[n], chorus_bins_q1[n], chorus_bins_q3[n]  = avg_func(chorus_bin, method=method)

        if method=='median':
            # Get 1st (25%) quartile
            chorus_bins_q1[n] = np.nanmedian(chorus_bin[chorus_bin 
                                                        < np.nanmedian(chorus_bin)])
            # Get 3rd (75%) quartile
            chorus_bins_q3[n] = np.nanmedian(chorus_bin[chorus_bin 
                                                        > np.nanmedian(chorus_bin)])

        if method=='mean':
            # Get standard deviation
            chorus_bins_q1[n] = chorus_bins[n] - np.nanstd(chorus_bin)
            chorus_bins_q3[n] = chorus_bins[n] + np.nanstd(chorus_bin)
   
        if method=='gmean':
            # Get geometric mean
            chorus_bins_q1[n] = chorus_bins[n]/stats.gstd(chorus_bin)
            chorus_bins_q3[n] = chorus_bins[n]*stats.gstd(chorus_bin)
            
    return delay_bins, chorus_bins, chorus_bins_q1, chorus_bins_q3, statistics

def create_chorus_selector(chorus:np.ndarray, lbc:np.ndarray, ubc:np.ndarray,
                           general_selector:np.ndarray) -> dict:
    """Function to select chorus based on general selector and 
    make sure it is larger than zero and not infinite/nan
    """

    # Lower and upper band together
    chorus_selector = (np.isfinite(chorus) & (chorus > 0)
                       & general_selector)

    # Lower band
    lbc_selector = (np.isfinite(lbc) & (lbc > 0)
                       & general_selector)

    # Upper band
    ubc_selector = (np.isfinite(ubc) & (ubc > 0)
                       & general_selector)
    
    return_dict = {'chorus_selector' : chorus_selector,
                   'lbc_selector' : lbc_selector,
                   'ubc_selector' : ubc_selector}
    
    return return_dict

def create_plotting_data(delay:np.ndarray, chorus:np.ndarray,
                         chorus_selector:np.ndarray,
                         lbc:np.ndarray, lbc_selector:np.ndarray,
                         ubc:np.ndarray, ubc_selector:np.ndarray,
                         bin_size:int, cutoff:int, method:str='peak'):
    """Get statistical plotting data for both electic and magnetic fields
    INPUT
    delay - delay from start of quiet period
    chorus - full ubc + lbc chorus measurements
    chorus_selector - indices of good chorus measurements
    lbc - lower-band chorus measurements
    lbc_selector - indices of good lbc measurements
    ubc - upper-band chorus measurements
    ubc_selector - indices of good ubc measurements
    bin_size - how many minutes to bin delay by
    cutoff - chorus measurements less than this will be ignored
    method - how to get "average" of bin_size
    OUTPUT
    return_dict - dictionary with data to plot
    """

    # Calculate average over bin for full chorus
    (delay_chorus_bins, chorus_bins,
     chorus_bins_q1, chorus_bins_q3,
     statistics) = average_bins(delay[chorus_selector], chorus[chorus_selector],
                                method=method, bin_size=bin_size,
                                cutoff=cutoff)

    # Calculate for LBC
    (delay_lbc_bins, lbc_bins,
     lbc_bins_q1, lbc_bins_q3,
     lbc_statistics) = average_bins(delay[lbc_selector], lbc[lbc_selector],
                                    method=method, bin_size=bin_size,
                                    cutoff=cutoff)

    # Calculate for UBC
    (delay_ubc_bins, ubc_bins,
     ubc_bins_q1, ubc_bins_q3,
     ubc_statistics) = average_bins(delay[ubc_selector], ubc[ubc_selector],
                                    method=method, bin_size=bin_size,
                                    cutoff=cutoff)
    
    return_dict = {'delay_chorus' : delay_chorus_bins,
                   'chorus_bins' : chorus_bins,
                   'chorus_bins_q1_q3' : [chorus_bins_q1, chorus_bins_q3],
                   'chorus_statistics' : statistics,
                   'delay_lbc' : delay_lbc_bins,
                   'lbc_bins' : lbc_bins,
                   'lbc_bins_q1_q3' : [lbc_bins_q1, lbc_bins_q3],
                   'ubc_statistics' : ubc_statistics,
                   'delay_ubc' : delay_ubc_bins,
                   'ubc_bins' : ubc_bins,
                   'ubc_statistics' : ubc_statistics}
    
    return return_dict

def bootstrap_slope_error(x:np.ndarray, y:np.ndarray, n_samples:int) -> float:
    """ Function to calculate the error in a slope
    estimated via a linear regression using the bootstrap method. 
    Resample data many times and find the slope each time.
    Standard deviations in these slopes will be the error.
    Requires sklearn.utils.resample
    INPUT
    x - independent variable
    y - dependent variable
    n_samples -  number of resamplings to perform
    OUTPUT
    std - standard deviation in the slope
    """
    
    # Array to store slopes in
    sampled_slopes = np.zeros(n_samples)

    for i in range(n_samples):

        # Resample data
        new_x, new_y = resample(x, y, #n_samples=int(len(y)/2),
                                replace=True)

        # Get slope
        slope, intercept, r_value, p_value, std_err = stats.linregress(new_x, new_y)

        # Write to array
        sampled_slopes[i] = slope
        
    # # Get the confidence interval of the distribution
    # confidence = 0.95
    # low, high = np.sort(sampled_slopes)[[int(len(sampled_slopes)*1-confidence),
    #                                      int(len(sampled_slopes)*confidence)]]
        
    return np.std(sampled_slopes)