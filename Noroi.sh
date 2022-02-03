#!/bin/bash

# filelist contains a list of the full path to each of the weekly
# Fermi datafiles to process
# 
# gteselect outputs a file named gtselectE1_E2GeVp202.fits
# gtmktime takes this file and produces gtmktimeE1_E2GeVp202.fits
# both files are written in the directory where the script is run


#########Activating conda environment############################
source ~/anaconda3/etc/profile.d/conda.sh
#########Activating fermitools###################################
conda activate fermi


gtselect infile=@filelist.txt outfile=gtseclect10_20GeVp202.fits \
	ra=INDEF dec=INDEF rad=INDEF evclass=3 tmin=INDEF tmax=INDEF \
	emin=10e3 emax=20e3 zmax=100
gtselect infile=@filelist.txt outfile=gtseclect20_30GeVp202.fits \
	ra=INDEF dec=INDEF rad=INDEF evclass=3 tmin=INDEF tmax=INDEF \
	emin=20e3 emax=30e3 zmax=100
gtselect infile=@filelist.txt outfile=gtseclect30_40GeVp202.fits \
	ra=INDEF dec=INDEF rad=INDEF evclass=3 tmin=INDEF tmax=INDEF \
	emin=30e3 emax=40e3 zmax=100
gtselect infile=@filelist.txt outfile=gtseclect40_50GeVp202.fits \
	ra=INDEF dec=INDEF rad=INDEF evclass=3 tmin=INDEF tmax=INDEF \
	emin=40e3 emax=50e3 zmax=100
gtselect infile=@filelist.txt outfile=gtseclect50_60GeVp202.fits \
	ra=INDEF dec=INDEF rad=INDEF evclass=3 tmin=INDEF tmax=INDEF \
	emin=50e3 emax=60e3 zmax=100

SPACECRAFT=/srv/data/fermi/spacecraft/p202/lat_spacecraft_merged.fits

gtmktime scfile=$SPACECRAFT \
	filter="DATA_QUAL==1 && LAT_CONFIG==1 && ABS(ROCK_ANGLE)<52" \
	roicut=no evfile=gtseclect10_20GeVp202.fits \
	outfile=gtmktime10_20GeVp202.fits

gtmktime scfile=$SPACECRAFT \
	filter="DATA_QUAL==1 && LAT_CONFIG==1 && ABS(ROCK_ANGLE)<52" \
	roicut=no evfile=gtseclect20_30GeVp202.fits \
	outfile=gtmktime20_30GeVp202.fits
gtmktime scfile=$SPACECRAFT \
	filter="DATA_QUAL==1 && LAT_CONFIG==1 && ABS(ROCK_ANGLE)<52" \
	roicut=no evfile=gtseclect30_40GeVp202.fits \
	outfile=gtmktime30_40GeVp202.fits
gtmktime scfile=$SPACECRAFT \
	filter="DATA_QUAL==1 && LAT_CONFIG==1 && ABS(ROCK_ANGLE)<52" \
	roicut=no evfile=gtseclect40_50GeVp202.fits \
	outfile=gtmktime40_50GeVp202.fits
gtmktime scfile=$SPACECRAFT \
	filter="DATA_QUAL==1 && LAT_CONFIG==1 && ABS(ROCK_ANGLE)<52" \
	roicut=no evfile=gtseclect50_60GeVp202.fits \
	outfile=gtmktime50_60GeVp202.fits


