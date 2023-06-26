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

where {probe} is either rbspa or rbspb and year is 2012, 2013, etc.

To download all ephemerides for the entire mission run the script at: src/data/van-allen-probes-ephemerides-data-downloader.py

## Data processing

### Injection finding
After you've downloaded the needed data, the next step is to find all injection/quiet periods using the SME index. To do this run the code in src/features/find_injections_with_sme.py

This will create a text file with quiet period start, quiet end, injection start, quiet length, and injection length. This file is located at data/interim/sme-injections-quiet-times.txt

There is also an associated notebook under notebooks/reports/plot-sme-injections.ipynb that creates some plots with this file.

### Matching injections to probe location



