[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

# juno-waves_flux_density_calculation
Routines for calculating the Juno/Waves flux density

# Determination of the calibration gain

To obtain the calibration gain:

Please run the OBTAIN_CALIBRATION_GAIN.pro routine which will do the following steps: 
	
	OBTAIN_CALIBRATION_GAIN,/creation_spdyn_sav,/all_calibration, version=v
	
(with v=1 to reproduce the results from the [Louis et al., 2021](https://doi.org/10.1029/2021JA029435) paper)

1) Creation of spdyn_sav data using:
	PLOT_SPDYN_SURVEY,yyyydddb
		-> contains linear and db data: 	rdb or rlin: raw data
								zdb or zlin: FFT filtered data
								zdb2 or zlin2: FFT filtered daily background subtract data
	it will save the daily data under */data_n1>$#years/$#yyyyddd_spdyn_v$#v.sav*
2) Make time series for each frequency. It will used :
		ALL_TIMESERIES,’all’,delta_t=60,/linear,/noback_subtract
		it will call MAKE_TIMESERIES and used zlin data (FFT filtered)
		and save the time series under:
			juno>calibration>make_timeseries>ALL_linear_noback_subtract>ALL_timeseries_d60_channels_0-109_zlin.sav
3) Create a background using:
		MAKE_BACKGROUND_PJ. It will call MAKE_BACKGROUND & BACKGROUND
4) Subtract background, using:
		SUBTRACT_BACKGROUND
5) Create the time series for calibration, thus containing the ephemeris, using:
		CALIBRATION_READ,0,109,60,/all,/linear
			restore,>make_timeseries>ALL_linear_noback_subtract>ALL_timeseries_d60_channels_0-109_zlin.sav
6) Select only interval where Dist_Juno > 30 R_Jupiter and abs(Magnetic_latitude)<15°
		CALIBRATION_WRITE2,60,/all,/linear
7) Determination of the calibration gain, using:
		INTERCAL,freq,gain

# Calibrating the data

Please use the CALIBRATION_PROCESS.PRO routine
	
	calibration_process,YYYYDDDb, YYYYDDDe,/reso1sec,/create_savefile
 
keyword `/create_savefile` needs to be added if you want to calibrate the data of a new day whose saveset has not yet been calculated\\

keyword `version=v`(with `v` the number of the version) is 1 by default ([Louis et al., 2021](https://doi.org/10.1029/2021JA029435)). Note that currently, the "[JUNO E/J/S/SS WAVES CALIBRATED SURVEY FULL RESOLUTION V1.0, JNO-E/J/SS-WAV-3-CDR-SRVFULL-V1.0](https://doi.org/10.17189/1519710)", which are the PDS file used for `version=01`. Only the "[JUNO E/J/S/SS WAVES CALIBRATED SURVEY FULL RESOLUTION V2.0, JNO-E/J/SS-WAV-3-CDR-SRVFULL-V2.0](https://doi.org/10.17189/1520498)" are currently available through PDS (i.e., `version=02` in our pipeline).


# Creation of cdf calibrated file

To create cdf file of the processed data, please used the `make_cdf_juno_waves_calibrated.py` python routine

# Requirements

Required folder tree structure:
Data need to be stored in directories following this template: '*../../data_PDS/WAVES_SURVEY/YYYYDOY_ORBIT_##*' with '*##*' the orbit number

idl 8.5

Python 3.7.4

Python packages:

os

datetime

scipy==1.3.1

numpy==1.17.2

argparse == 1.1

astropy == 3.2.2   	

json == 2.0.9

requests == 2.22.0

spacepy == 0.2.2

pycdf
