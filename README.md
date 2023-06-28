# Read Me File for Project
see document-struture.md for a layout of the project structure and where files are located
you will need to create reports/figures/, logs/, and data/ directories as these are included in my .gitignore file.

## Setup
1. When running code make sure you are in the base directory as this will ensure that all the code runs as expected.
2. I recommend creating a new virtual environment and installing all the dependencies through pip3 with the requirements.txt file.

## 1. Data
If you are running this for the first time you will want to download the data.

### 1.1 SME Index
The SME index is the supermag implementation of the AE index. You can download this data here: https://supermag.jhuapl.edu/indices/?layers=SME.E&fidelity=low&start=2012-01-30T06%3A00%3A00.000Z&step=14400&tab=download

This data is stored under data/raw/sme

### 1.2 Van Allen Probes Magnetic Ephemerides
We need to know the location of the Van Allen Probes satellites, thus we need the ephemerides data from the mission. This can be found here: https://cdaweb.gsfc.nasa.gov/pub/data/rbsp/{probe}/ephemeris/ect-mag-ephem/hdf5/def-1min-t89q/{year}/

where {probe} is either rbspa or rbspb and {year} is 2012, 2013, etc.

To download all ephemerides for the entire mission run the script at: src/data/van-allen-probes-ephemerides-data-downloader.py

I've tried to transition my code to use the official NASA CDAWeb data repository, however this was acting quite slow at the time of testing. If you do have issues the data may also be available via: https://emfisis.physics.uiowa.edu/Flight/RBSP-A/LANL/MagEphem/ or https://rbsp-ect.newmexicoconsortium.org/data_pub/

### 1.3 Van Allen Probes EMFISIS data


## 2. Data processing

### 2.1 Injection finding
After you've downloaded the needed data, the next step is to find all injection/quiet periods using the SME index. To do this run the code in src/features/find_injections_with_sme.py

This will create a text file with quiet period start, quiet end, injection start, quiet length, and injection length. This file is located at data/interim/sme-injections-quiet-times.txt

There is also an associated notebook under notebooks/reports/plot-sme-injections.ipynb that creates some plots with this file.

### 2.2 Matching injections to probe location
After you've identified all of the injection/quiet periods we need to find the location of the Van Allen Probes during these periods. Using this information we can pick out which times we should download the EMFISIS data for. To perform this analysis run the script at: src/features/match-probe-location-to-injection.py

### 2.3 Download EMFISIS data
When you have the list of good probe times during injections you should download this data, see data section 1.3 on doing this.

### 2.4 Process and compile the data
With all of the raw data downloaded you can now process the data and store just the parts needed for the analysis. The process is done in the script located at: src/features/van-allen-probe-injection-data-compiler.py. The output of this script is a .h5 file located at data/processed/chorus-delay-data.h5.

The following is a brief overview of what the code does:

This script combines all the data that we've downloaded and processed. It does this through several functions. read_rbsp_sheath_corrected_psd reads in a local CDF file from the EMIFISIS instrument and extracts the power spectral density (PSD), density, location information, and associated time and frequency of the specified measurement. read_rbsp_emfisis_b_field reads in a magnetic field amplitude CDF file and extracts the magnetic field strength and associated time. read_process_rbsp_data gets the PSD, magnetic field strength, and density for a specified probe and date using the previously described functions. If there is no density data file or there is a file, but no data in it the function returns nan values. filter_write_to_dict filters the PSD to lower and upper band and only selects times that meet the specified threshold value.

The notebook uses these functions to loop through every quiet time and extract the PSD, magnetic field strength, density, and magnetic ephemerides for each probe. If the quiet period spans more than 1 day it will read in both days and concatenate the data.

The code then filters the data to the times during the quiet period and checks to make sure there is data after this filtering. It then creates interpolated functions for the magnetic field strength, density, and emphemerides data. 

The code then loops through each time within the a quiet period and first checks if the density is low enough. If it is, the code uses the function filter_write_to_dict to find the fce using the B-field magnitude for the time. Using fce it filters the PSD for that time to lower and upper band chorus. It then checks if the max PSD value within these frequency ranges is larger than the specified threshold. If it is, the function integrates the PSD over the frequency ranges and writes the integrated values and ephemerides information to a dictionary.

## 3. Data analysis

### 3.1 Create the data to visualize
Before creating any plots or analysis you need to compile just the data you need and write this to another h5 file. This isn't entirely neccessary, but speeds up and similifies the analysis/plotting process. To do this run the script at src/features/create-plotting-data.py

You may want to run this for a psd type of both integrated and max. This just specifies the chorus measurement is based on the max psd value for a timestep or the integrated psd. 

### 3.2




