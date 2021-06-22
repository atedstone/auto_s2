# README

At time of writing (June 2021), this pipeline is set up to run on beo05 at UniFr.


## Installation

Put all scripts on the system.

Requires a conda environment (`auto_s2`) with, amongst other packages:

* SentinelSat
* Geoutils

Set up a cron job on the master login node:

    crontab -e
    0 10 * * * /home/tedstona/scripts/auto_s2/setup_qsub_s2.sh


## Overview of files

* setup_qsub_s2.sh : this sets up the paths needed for the cron scheduler to find qsub.
* run_auto_s2.sh : this is the queued job script. It takes the form of a regular job submission script.
* setup_sentinel2.sh : env. variables for working with Sentinel-2 data.
* auto_sentinel2.py : the Python script which processes the data.
* sentinel2_tools.py : originally developed with/for Joe Cook, this has tools relating to downloading.

