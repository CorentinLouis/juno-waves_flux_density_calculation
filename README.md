[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

# juno-waves_flux_density_calculation
Routines for calculating the Juno/Waves flux density


To obtain the calibration gain:
1) You first need to have daily dynamic spectrum:
	Creation of spdyn_sav data using:
	PLOT_SPDYN_SURVEY,yyyydddb
		-> contains linear and db data: 	rdb or rlin: raw data
								zdb or zlin: FFT filtered data
								zdb2 or zlin2: FFT filtered daily background subtract data
	it will save the daily data under juno>stage>spdyn_sav_data>yyyyddd.sav


2) Then please use the OBTAIN_CALIBRATION_GAIN.pro routine which will do the following steps: 
	1) Make time series for each frequency. It will used :
		ALL_TIMESERIES,’all’,delta_t=60,/linear,/noback_subtract
		it will call MAKE_TIMESERIES and used zlin data (FFT filtered)
		and save the time series under:
			juno>calibration>make_timeseries>ALL_linear_noback_subtract>ALL_timeseries_d60_channels_0-109_zlin.sav
	2) Create a background using:
		MAKE_BACKGROUND_PJ. It will call MAKE_BACKGROUND & BACKGROUND
	3) Subtract background, using:
		SUBTRACT_BACKGROUND
	4) Create the time series for calibration, thus containing the ephemeris, using:
		CALIBRATION_READ,0,109,60,/all,/linear
			restore,>make_timeseries>ALL_linear_noback_subtract>ALL_timeseries_d60_channels_0-109_zlin.sav
	5) Select only interval where Dist_Juno > 30 R_Jupiter and abs(Magnetic_latitude)<15°
		CALIBRATION_WRITE2,60,/all,/linear
	6) Determination of the calibration gain, using:
		INTERCAL,freq,gain
