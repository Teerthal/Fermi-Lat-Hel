#!/bin/bash

SPACECRAFT=/srv/data/fermi/spacecraft/p202/lat_spacecraft_merged.fits

gtltcube evfile=gtmktime10_20GeVp202.fits scfile=$SPACECRAFT \
	outfile=cube10_20GeVp202.fits dcostheta=0.025 binsz=1

gtltcube evfile=gtmktime20_30GeVp202.fits scfile=$SPACECRAFT \
	outfile=cube20_30GeVp202.fits dcostheta=0.025 binsz=1

gtltcube evfile=gtmktime30_40GeVp202.fits scfile=$SPACECRAFT \
	outfile=cube30_40GeVp202.fits dcostheta=0.025 binsz=1

gtltcube evfile=gtmktime40_50GeVp202.fits scfile=$SPACECRAFT \
	outfile=cube40_50GeVp202.fits dcostheta=0.025 binsz=1

gtltcube evfile=gtmktime50_60GeVp202.fits scfile=$SPACECRAFT \
	outfile=cube50_60GeVp202.fits dcostheta=0.025 binsz=1


gtexpmap evfile=gtmktime10_20GeVp202.fits scfile=$SPACECRAFT \
	expCube=cube10_20GeVp202.fits outfile=exp10_20GeVp202.fits \
	irfs=P7CLEAN_V6 srcrad=180 nlong=720 nlat=360 nenergies=5


gtexpmap evfile=gtmktime20_30GeVp202.fits scfile=$SPACECRAFT \
	expCube=cube20_30GeVp202.fits outfile=exp20_30GeVp202.fits \
	irfs=P7CLEAN_V6 srcrad=180 nlong=720 nlat=360 nenergies=5


gtexpmap evfile=gtmktime30_40GeVp202.fits scfile=$SPACECRAFT \
	expCube=cube30_40GeVp202.fits outfile=exp30_40GeVp202.fits \
	irfs=P7CLEAN_V6 srcrad=180 nlong=720 nlat=360 nenergies=5


gtexpmap evfile=gtmktime40_50GeVp202.fits scfile=$SPACECRAFT \
	expCube=cube40_50GeVp202.fits outfile=exp40_50GeVp202.fits \
	irfs=P7CLEAN_V6 srcrad=180 nlong=720 nlat=360 nenergies=5


gtexpmap evfile=gtmktime50_60GeVp202.fits scfile=$SPACECRAFT \
	expCube=cube50_60GeVp202.fits outfile=exp50_60GeVp202.fits \
	irfs=P7CLEAN_V6 srcrad=180 nlong=720 nlat=360 nenergies=5







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

python ./cutblat.py

