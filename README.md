# **NYC-Chem**

WRF-GHG and STILT project involving HALO CH4 mission over NYC Jul 20 - Aug 12 2023.

Validation code and helper scripts (submit, writescript, etc.) written by Steve Jester.

WRF-GHG 4.5.2 (?) by NOAA/NCAR/UCAR with modifications by Xiao-Ming Hu and Steve Jester.

For more information please email Steve Jester at steve@belumenus.com

## *Major files*

### plot-atmosphere.ipynb

This Jupyter Notebook is used to provide a sanity check on the WRF simulation. Most options and procedures are described in the notbook. Works only for meteorology variables as currently written (10/19/2024). Support for GHG variables planned.

### validate.ipynb

This Jupyter Notebook is used to validate a WRF simulation against observations. Currently, surface and upper air tests are working properly, however users are advised to make sure their upper air observations are complete prior to testing. Satelite classes are in progress (see *validate.py*). Proper procedures are in progress.

## *Helper scripts/notebooks*

### writenamelist.py

This Python script is setup to write namelists for a multi-stage WRF run. All the primary options are included in the first block of statements before going into the for loop.

### writescript.py

This Python script writes out submission scripts to submit a multi-stage WRF run to a cluster HPC using SLURM. Currently is is written for use with Stampede3 hosted at the University of Texas however this can be edited. Most options are included in the first block of statements before going into the for loop. Some statements in the for loop need to be edited based on your setup and HPC cluster.

### validate.py

This Python script is a prototype for *validate.ipynb*. Some classes from this file have not been transfered to the notebook and have not been tested yet.

### domainplot.ipynb

This Jupyter Notebook is a minimal notebook to plot the domain of a WRF file and superimpose a STILT domain on top of it. Setup and procedures are described in the notebook.

### submit.sh

This is a Bash script used for systematically submitting a multipart WRF simulation to Stampede3. The '623' mentioned is the number of characters needed to be skipped in order to get the proper job id for setting up dependencies. This will need to be edited for any HPC cluster other than Stampede3.