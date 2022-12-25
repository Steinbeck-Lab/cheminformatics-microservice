## Purpose
This directory has python files to help:
- downloading NMRShiftDB database as a collection of NMReData files.
- downloading the experimental data from the links in the NMReData.
- Unzipping and structuring the files and folders in a manner complying with nmrXiv submission requirements.

## Running
After navigating to the current folder "nmrshiftdb", run: python3 nmrshiftdb-download.py

## Structure
The folder contains two python files:
- helpers.py has the functions
- nmrshiftdb-download.py is the script which we run.
- output folder generated after running the script, containing the NMReData, along with two folders:
  - without_raw: has the NMReData of the molecules that doesn't have raw data
  - with_raw: has the NMReData of the molecules that have raw data along with the experimental data.

### with_raw folder:
This is the main folder where the data of interest for submission to nmrXiv exist. it has subfolders named after the molecules' ids in NMRShiftDB, and within the subfolders, there are more subfolders named after the spectra ids in NMRShiftDB.
