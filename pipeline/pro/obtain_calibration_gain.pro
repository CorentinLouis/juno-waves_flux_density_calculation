PRO OBTAIN_CALIBRATION_GAIN,all_calibration=all_calibration,$
							creation_spdyn_save=creation_spdyn_save,$
							make_timeseries=make_timeseries,$
							make_background=make_background,subtract_background=subtract_background,$
							calibration_timeseries=calibration_timeseries,calibration_data=calibration_data, $
							calibration_gain=calibration_gain,$
							yyyydddb=yyyydddb,yyyyddde=yyyyddde,$
							raw_data=raw_data,obtain_daily_min_max_values=obtain_daily_min_max_values,$
							linear=linear,dB=dB,version=version

if ~keyword_set(yyyydddb) then yyyydddb=2016100
if ~keyword_set(yyyyddde) then yyyyddde=2019174

if ~keyword_set(linear) and ~keyword_set(dB) then linear=1

if ~keyword_set(version) then version=1
version_name='_v'+string(format='(I02)',version)
if version eq 1 then version_dailysaveset = '' else version_dailysaveset = version_name


if keyword_Set(obtain_daily_min_max_values) then begin
	rdbmin=[]
	zdbmin=[]
	dbzmin=[]
	rdbmin=[]
	rlinmin=[]
	zlinmin=[]
	zdbmin=[]

	rdbmax=[]
	zdbmax=[]
	dbzmax=[]
	rdbmax=[]
	rlinmax=[]
	zlinmax=[]
	zdbmax=[]  
	tt_evo_signal=[]
	for i=100L,366L do begin
 		file_day = file_search('spdyn_sav_data/','2016'+strtrim(i,2)+'_spdyn'+version_dailysaveset+'.sav')
 		if file_day[0] ne '' then begin
 			print,file_day[0]
 			restore,file_day[0]
 			wn0=where(rdb gt 0)
 			rdbmin=[rdbmin,min(rdb[wn0])]
			zdbmin=[zdbmin,min(zdb[wn0]+200)]
			rlinmin=[rlinmin,min(rlin[wn0])]
			zlinmin=[zlinmin,min(zlin[wn0])]

			rdbmax=[rdbmax,max(rdb[wn0])]
			zdbmax=[zdbmax,max(zdb[wn0])]
			rlinmax=[rlinmax,max(rlin[wn0])]
			zlinmax=[zlinmax,max(zlin[wn0])]
			tt_evo_signal=[tt_evo_signal,i]
 		endif
 	endfor
 	stop
	return
endif

if keyword_set(all_calibration) then begin
	make_timeseries=1
	make_background=1
	subtract_background=1
	calibration_timeseries=1
	calibration_data=1
	calibration_gain=1
endif

!path='/Users/serpe/Volumes/kronos/juno/stage/IDL_pro_util:'+!path





if keyword_set(raw_data) then begin

	if keyword_Set(linear) then begin
		data_type_treat='lin'
		data_type_folder='linear'
	endif else if keyword_Set(dB) then begin
		data_type_treat='db'
		data_type_folder='db'
	endif
	if keyword_set(make_timeseries) then begin
		ALL_TIMESERIES,'',yyyydddb=yyyydddb,yyyyddde=yyyyddde,delta_t=60,LINEAR=LINEAR,dB=dB,/raw,version=version
	endif

	if keyword_set(subtract_background) then begin
		print,' ### Subtraction (in linear scale) of the uniform background ###'
		restore,'./make_timeseries/ALL_'+data_type_folder+'_raw_noback_subtract/ALL_timeseries_d60_channels_16-125_r'+data_type_treat+'.sav'
		
		wfreq=where(frequencies eq 6500.00)
		plotps,'pdf/timeseries_'+strtrim(long(frequencies[wfreq]),2)+'kHz_r'+data_type_treat
			plot_io,time,timeseries[wfreq,*],/xsty,/ysty,xtit=TeXtoIDL('Time (Day of Year 2016)'),ytit=TexToIDL('Intensity (V^2/m^2/Hz)'),tit=TexToIDL('Timeseries @ '+strtrim(frequencies[wfreq],2)+' kHz');,yrange=[1e-24,1e-10]
		endps,filename='pdf/timeseries_'+strtrim(long(frequencies[wfreq]),2)+'kHz_r'+data_type_treat
	



		restore,'./background/background_ALLPJ_z'+data_type_treat+'.sav',/verb
	
		timeseries2=transpose(timeseries)
		SUBTRACT_BACKGROUND, background,sigma,0, timeseries2
		timeseries=transpose(timeseries2)
	
		save,timeseries,time,frequencies,filename='./make_timeseries/ALL_'+data_type_folder+'_raw_back_subtract/ALL_timeseries_d60_channels_16-125_r'+data_type_treat+'nb.sav'
	endif

	if keyword_set(calibration_timeseries) then $
		calibration_read,0+16,109+16,60,LINEAR=LINEAR,dB=dB,/raw,/all,/NOBACK_SUBTRACT

	if keyword_set(calibration_data) then begin
		print,' ### Select only the interval where D_Juno>Rmin=30 R_Jupiter and |Mlat|<MLatMax15° ###'
		CALIBRATION_WRITE2,60,/all,LINEAR=LINEAR,dB=dB,/raw,Rmin=30.,MLatMax=15.,/norm1AU,/NOBACK_SUBTRACT

		restore,'./calibration_data/ALL_calibration_data_freq_'+data_type_folder+'_raw_noback_subtract_d60_r'+data_type_treat+'_Rmin30_MLatmax15.sav',/verb
		plotps,'pdf/Juno-Waves_50percent_level_rawdata_r'+data_type_treat
			loadct,0
			plot_oo,fx,xmed[0,*],/xsty,yrange=[1e-23,1e-15],/ysty,/nodata,ytit=TexToIDL('Flux density @ 1 AU (W/m^2/Hz)'),xtit='Frequency (kHz)',tit='Juno/Waves 50% level - raw data'
			loadct,53
			oplot,[fx[27],fx[27]],[1e-23,1e-15],col=120,linestyle=2
			oplot,[fx[45],fx[45]],[1e-23,1e-15],col=120,linestyle=2
			oplot,[fx[71],fx[71]],[1e-23,1e-15],col=120,linestyle=2
			oplot,[fx[72],fx[72]],[1e-23,1e-15],col=120,linestyle=2
			al_legend,'LFR-Low',position=[2.5,1e-15*0.8],box=0,textcolors=120
			al_legend,'LFR-High',position=[25,1e-15*0.8],box=0,textcolors=120
			al_legend,'HFR-Low',position=[350,1e-15*0.8],box=0,textcolors=120
			al_legend,'HFR-High',position=[6000,1e-15*0.8],box=0,textcolors=120
			loadct,1
			for i=0,n_elements(xmed[*,0])-1 do oplot,fx,(xmed[i,where(xmed[i,*] gt 0.)]),col=220
			loadct,5
			xmedmed=fltarr(n_elements(xmed[0,*]))
			for i=0,n_elements(xmedmed)-1 do xmedmed[i]=median(xmed[where(xmed[*,i] gt 0.),i])
			oplot,fx,xmedmed,col=100
			loadct,0
		endps,filename='pdf/Juno-Waves_50percent_level_rawdata_r'+data_type_treat
	
	
		plotps,'pdf/Juno-Waves_1percent_level_rawdata_r'+data_type_treat
			loadct,0
			plot_oo,fx,x1[0,*],/xsty,yrange=[1e-23,1e-15],/ysty,/nodata,ytit=TexToIDL('Flux density @ 1 AU (W/m^2/Hz)'),xtit='Frequency (kHz)',tit='Juno/Waves 1% level - raw data'
			loadct,53
			oplot,[fx[27],fx[27]],[1e-23,1e-15],col=120,linestyle=2
			oplot,[fx[45],fx[45]],[1e-23,1e-15],col=120,linestyle=2
			oplot,[fx[71],fx[71]],[1e-23,1e-15],col=120,linestyle=2
			oplot,[fx[72],fx[72]],[1e-23,1e-15],col=120,linestyle=2
			al_legend,'LFR-Low',position=[2.5,1e-15*0.8],box=0,textcolors=120
			al_legend,'LFR-High',position=[25,1e-15*0.8],box=0,textcolors=120
			al_legend,'HFR-Low',position=[350,1e-15*0.8],box=0,textcolors=120
			al_legend,'HFR-High',position=[6000,1e-15*0.8],box=0,textcolors=120
			loadct,1
			for i=0,n_elements(x1[*,0])-1 do oplot,fx,(x1[i,where(x1[i,*] gt 0.)]),col=220
			loadct,5
			x1med=fltarr(n_elements(x1[0,*]))
			for i=0,n_elements(x1med)-1 do x1med[i]=median(x1[where(x1[*,i] gt 0.),i])
			oplot,fx,x1med,col=100
			loadct,0
		endps,filename='pdf/Juno-Waves_1percent_level_rawdata_r'+data_type_treat
	endif

	return
endif

;# 1) Creation of spdyn_sav data using:
;# 	PLOT_SPDYN_SURVEY,yyyydddb,/rebin15sec
;# 		-> contains linear and db data: 	rdb or rlin: raw data
;# 								zdb or zlin: FFT filtered data
;# 								zdb2 or zlin2: FFT filtered daily background subtract data
;# 	it will save the daily data under juno>stage>spdyn_sav_data>yyyyddd.sav

if keyword_Set(linear) then begin
	data_type_treat='lin'
	data_type_folder='linear'
endif else if keyword_Set(dB) then begin
	data_type_treat='db'
	data_type_folder='db'
endif

cmd = 'mkdir pdf'
spawn,cmd,resu
cmd = 'mkdir ../background'
spawn,cmd,resu

if keyword_set(creation_spdyn_save) then begin
	doy=[]
	doy=indgen(366)+2016001
	doy=doy[where(doy ge 2016100)]
	doy=[doy,indgen(365)+2017001,indgen(365)+2018001,indgen(365)+2019001]
	for idoy=0,n_elements(doy)-1 do CREATE_SAVEFILE_SPDYN_SURVEY,doy[idoy],/rebin15sec,version=version
endif

;# 2) creation of timeseries from linear FFT filtered data
if keyword_set(make_timeseries) then begin
	print,' ### Creation of timeseries from linear FFT filtered data ###'
	ALL_TIMESERIES,'',yyyydddb=yyyydddb,yyyyddde=yyyyddde,delta_t=60,LINEAR=LINEAR,dB=dB,/noback_subtract,version=version
endif
;# 3) creation of a background 
if keyword_set(make_background) then begin
	print,' ### Creation of a background ###'

	make_background_PJ,/allPJ,LINEAR=LINEAR,dB=dB,version=version
endif


;# 4) subtraction (in linear scale) of the uniform background
if keyword_set(subtract_background) then begin
	print,' ### Subtraction (in linear scale) of the uniform background ###'

	restore,'../make_timeseries/ALL_'+data_type_folder+'_noback_subtract/ALL_timeseries_d60_channels_16-125_z'+data_type_treat+version_name+'.sav'
	
 	wfreq=where(frequencies eq 6500.00)
	plotps,'pdf/timeseries_'+strtrim(long(frequencies[wfreq]),2)+'kHz'+version_name
		plot_io,time,timeseries[wfreq,*],/xsty,/ysty,xtit=TeXtoIDL('Time (Day of Year 2016)'),ytit=TexToIDL('Intensity (V^2/m^2/Hz)'),tit=TexToIDL('Timeseries @ '+strtrim(frequencies[wfreq],2)+' kHz');,yrange=[1e-24,1e-10]
	endps,filename='pdf/timeseries_'+strtrim(long(frequencies[wfreq]),2)+'kHz'+version_name
	
	restore,'../background/background_ALLPJ_z'+data_type_treat+version_name+'.sav',/verb


	plotps,'pdf/background_allmission_'+data_type_treat+version_name
		plot_oo,frequencies,background,/xsty,/ysty,xtit=TeXtoIDL('Frequency (kHz)'),ytit=TexToIDL('Intensity (V^2/m^2/Hz)'),tit=TexToIDL('Background ['+strtrim(yyyydddb,2)+'-'+strtrim(yyyyddde,2)+']'),yrange=[1e-18,1e-12]
	endps,filename='pdf/background_allmission_'+data_type_treat+version_name
	openw,u,'pdf/background_'+data_type_treat+version_name+'.txt',/get_lun
		printf,u,'Channels	Frequency	Background'
		printf,u,TexToIDL('		  (kHz)		(V^2/m^2/Hz)')
		for i=0,n_elements(frequencies)-1 do printf,u,strtrim(i+16,2)+'		'+strtrim(frequencies[i],2)+'		'+string(format='(e9.3)',background[i])
	close,u
	free_lun,u


	timeseries2=transpose(timeseries)
	SUBTRACT_BACKGROUND, background,sigma,0, timeseries2
	timeseries=transpose(timeseries2)

	plotps,'pdf/timeseries_'+strtrim(long(frequencies[wfreq]),2)+'kHz_back-subtracted_'+data_type_treat+version_name
		wn0=where(timeseries[wfreq,*] gt 0.)  
		plot_io,time,timeseries[wfreq,*],/xsty,/ysty,xtit=TeXtoIDL('Time (Day of Year 2016)'),ytit=TexToIDL('Intensity (V^2/m^2/Hz)'),tit=TexToIDL('Timeseries @ '+strtrim(frequencies[wfreq],2)+' kHz - background subtracted'),yrange=[min(timeseries[wfreq[0],wn0]),max(timeseries[wfreq[0],wn0])]
	endps,filename='pdf/timeseries_'+strtrim(long(frequencies[wfreq]),2)+'kHz_back-subtracted_'+data_type_treat+version_name
	
	cmd = 'mkdir ../make_timeseries/ALL_'+data_type_folder+'_back_subtract'
	spawn,cmd,resu
	save,timeseries,time,frequencies,filename='../make_timeseries/ALL_'+data_type_folder+'_back_subtract/ALL_timeseries_d60_channels_16-125_z'+data_type_treat+'nb'+version_name+'.sav'
endif



;# 5) Create the time series for calibration, thus containing the ephemeris, using:
if keyword_set(calibration_timeseries) then begin
	print,' ### Creation of calibration timeseries (adding ephemeris to timeseries) ###'
	CALIBRATION_READ,0+16,109+16,60,/all,LINEAR=LINEAR,dB=dB,version=version
endif


;# 6) Select only the interval where D_Juno>Rmin=30 R_Jupiter and |Mlat|<MLatMax15°
if keyword_set(calibration_data) then begin
	print,' ### Select only the interval where D_Juno>Rmin=30 R_Jupiter and |Mlat|<MLatMax15° ###'
	cmd = 'mkdir ../calibration_data'
	spawn,cmd,resu
	CALIBRATION_WRITE2,60,/all,LINEAR=LINEAR,dB=dB,Rmin=30.,MLatMax=15.,/norm1AU,version=version
	restore,'../calibration_data/ALL_calibration_data_freq_'+data_type_folder+'_back_subtract_d60_z'+data_type_treat+'nb_Rmin30_MLatmax15.sav',/verb
	plotps,'pdf/Juno-Waves_50percent_level_'+data_type_treat+version_name
		loadct,0
		plot_oo,fx,xmed[0,*],/xsty,/ysty,/nodata,ytit=TexToIDL('Flux density @ 1 AU (W/m^2/Hz)'),xtit='Frequency (kHz)',tit='Juno/Waves 50% level',yrange=[1e-24,1e-16],charsize=1.75
		loadct,53
		oplot,[fx[27],fx[27]],[1e-24,1e-16],col=120,linestyle=2
		oplot,[fx[45],fx[45]],[1e-24,1e-16],col=120,linestyle=2
		oplot,[3000,3000],[1e-24,1e-16],col=120,linestyle=2
		;oplot,[fx[71],fx[71]],[1e-24,1e-16],col=120,linestyle=2
		;oplot,[fx[72],fx[72]],[1e-24,1e-16],col=120,linestyle=2
		al_legend,'LFR-Low',position=[2.5,1e-16*0.8],box=0,textcolors=120
		al_legend,'LFR-High',position=[25,1e-16*0.8],box=0,textcolors=120
		al_legend,'HFR-Low',position=[350,1e-16*0.8],box=0,textcolors=120
		al_legend,'HFR-High',position=[6000,1e-16*0.8],box=0,textcolors=120
		loadct,1
		for i=0,n_elements(xmed[*,0])-1 do begin
			oplot,fx[0:44],xmed[i,0:44],col=220
			oplot,fx[45:-1],xmed[i,45:-1],col=220
			;wn0=where(xmed[i,0:44] gt 0.)
			;oplot,(fx[0:44])[wn0],(xmed[i,0:44])[wn0],col=220
			;wn0=where(xmed[i,72:-1] gt 0.)
			;oplot,(fx[72:-1])[wn0],(xmed[i,72:-1])[wn0],col=220
		endfor
		loadct,5
		xmedmed=fltarr(n_elements(xmed[0,*]))
		for i=0,n_elements(xmedmed)-1 do begin
			wn0=where(xmed[*,i] gt 0.)
			if wn0[0] ne -1 then xmedmed[i]=median(xmed[wn0,i])
		endfor
		oplot,fx[0:44],xmedmed[0:44],col=100
		oplot,fx[45:-1],xmedmed[45:-1],col=100
		;oplot,(fx[0:44]),(xmedmed[0:44]),col=100
		;oplot,(fx[72:-1]),(xmedmed[72:-1]),col=100
		loadct,0
	endps,filename='pdf/Juno-Waves_50percent_level_'+data_type_treat+version_name


	plotps,'pdf/Juno-Waves_1percent_level_'+data_type_treat+version_name
		loadct,0
		plot_oo,fx,x1[0,*],/xsty,/ysty,/nodata,ytit=TexToIDL('Flux density @ 1 AU (W/m^2/Hz)'),xtit='Frequency (kHz)',tit='Juno/Waves 1% level',yrange=[1e-24,1e-16],charsize=1.75
		loadct,53
		oplot,[fx[27],fx[27]],[1e-24,1e-16],col=120,linestyle=2
		oplot,[fx[45],fx[45]],[1e-24,1e-16],col=120,linestyle=2
		oplot,[3000,3000],[1e-24,1e-16],col=120,linestyle=2
		;oplot,[fx[71],fx[71]],[1e-24,1e-16],col=120,linestyle=2
		;oplot,[fx[72],fx[72]],[1e-24,1e-16],col=120,linestyle=2
		al_legend,'LFR-Low',position=[2.5,1e-16*0.8],box=0,textcolors=120
		al_legend,'LFR-High',position=[25,1e-16*0.8],box=0,textcolors=120
		al_legend,'HFR-Low',position=[350,1e-16*0.8],box=0,textcolors=120
		al_legend,'HFR-High',position=[6000,1e-16*0.8],box=0,textcolors=120
		loadct,1
		for i=0,n_elements(x1[*,0])-1 do begin
			oplot,fx[0:44],x1[i,0:44],col=220
			oplot,fx[45:-1],x1[i,45:-1],col=220
			;oplot,fx,x1[i,*],col=220
			;oplot,fx[0:44],(x1[i,0:44]),col=220
			;oplot,fx[72:-1],(x1[i,72:-1]),col=220
		endfor
		loadct,5
		x1med=fltarr(n_elements(x1[0,*]))
		for i=0,n_elements(x1med)-1 do begin
			wn0=where(x1[*,i] gt 0.)
			if wn0[0] ne -1 then x1med[i]=median(x1[wn0,i])
		endfor
		oplot,fx[0:44],x1med[0:44],col=100
		oplot,fx[45:-1],x1med[45:-1],col=100
		;oplot,fx,x1med,col=100
		;oplot,fx[0:44],x1med[0:44],col=100
		;oplot,fx[72:-1],x1med[72:-1],col=100
		loadct,0
	endps,filename='pdf/Juno-Waves_1percent_level_'+data_type_treat+version_name
endif

;# 7) Determination of the calibration gain
if keyword_set(calibration_gain) then begin
	print,'### Determination of the calibration gain ###'
	; # calculation of calibration gain for the subreceivers LFR-Low/High and HFR-High
	cmd = 'mkdir ../gain'
	spawn,cmd,resu
	INTERCAL, freq,gain,/PS,LINEAR=LINEAR,dB=dB,version=version
	
	; # calculation of calibration gain for the subreceiver HFR-Low
	make_timeseries,yyyydddb,yyyyddde,'PHC',61,87,60,'PHC',timeseries,time,LINEAR=LINEAR,dB=dB,/NOBACK_SUBTRACT,version=version
	calibration_write_pj,gain_new,temporal_window=2.,/ps,yyyydddb=yyyydddb,yyyyddde=yyyyddde,LINEAR=LINEAR,dB=dB,version=version
	
	restore,'gain_final_'+data_type_treat+version_name+'.sav'
	openw,u,'pdf/gain_final_'+data_type_treat+'.txt',/get_lun
		printf,u,'Channels	Frequency	Gain'
		printf,u,'		  (kHz)		'
		for i=0,n_elements(freq)-1 do printf,u,strtrim(i+16,2)+'		'+strtrim(freq[i],2)+'		'+string(format='(F7.3)',gain_final[i])
	close,u
	free_lun,u


	restore,'../calibration_data/ALL_calibration_data_freq_'+data_type_folder+'_back_subtract_d60_z'+data_type_treat+'nb_Rmin30_MLatmax15'+version_name+'.sav',/verb
	for i=0,n_elements(xmed[*,0])-1 do begin
		xmed[i,*]=xmed[i,*]*gain_final
		x1[i,*]=x1[i,*]*gain_final
	endfor

	plotps,'pdf/Juno-Waves_50percent_level_calibrated_'+data_type_treat+version_name
		loadct,0
		plot_oo,fx,xmed[0,*],/xsty,/ysty,/nodata,ytit=TexToIDL('Flux density @ 1 AU (W/m^2/Hz)'),xtit='Frequency (kHz)',tit='Juno/Waves 50% level after calibration',yrange=[1e-24,1e-15],charsize=1.75
		loadct,53
		oplot,[fx[27],fx[27]],[1e-24,1e-15],col=120,linestyle=2
		oplot,[fx[45],fx[45]],[1e-24,1e-15],col=120,linestyle=2
		oplot,[fx[71],fx[71]],[1e-24,1e-15],col=120,linestyle=2
		oplot,[fx[72],fx[72]],[1e-24,1e-15],col=120,linestyle=2
		al_legend,'LFR-Low',position=[2.5,1e-15*0.8],box=0,textcolors=120
		al_legend,'LFR-High',position=[25,1e-15*0.8],box=0,textcolors=120
		al_legend,'HFR-Low',position=[350,1e-15*0.8],box=0,textcolors=120
		al_legend,'HFR-High',position=[6000,1e-15*0.8],box=0,textcolors=120
		loadct,1
		for i=0,n_elements(xmed[*,0])-1 do begin
			oplot,fx[0:44],xmed[i,0:44],col=220
			oplot,fx[72:-1],xmed[i,72:-1],col=220
		endfor
		loadct,5
		xmedmed=fltarr(n_elements(xmed[0,*]))
		for i=0,n_elements(xmedmed)-1 do xmedmed[i]=median(xmed[*,i])
		oplot,fx[0:44],xmedmed[0:44],col=100
		oplot,fx[72:-1],xmedmed[72:-1],col=100
		loadct,0
	endps,filename='pdf/Juno-Waves_50percent_level_calibrated_'+data_type_treat+version_name


	plotps,'pdf/Juno-Waves_1percent_level_calibrated_'+data_type_treat+version_name
		loadct,0
		plot_oo,fx,x1[0,*],/xsty,/ysty,/nodata,ytit=TexToIDL('Flux density @ 1 AU (W/m^2/Hz)'),xtit='Frequency (kHz)',tit='Juno/Waves 1% level after calibration',yrange=[1e-24,1e-15],charsize=1.75
		loadct,53
		oplot,[fx[27],fx[27]],[1e-24,1e-15],col=120,linestyle=2
		oplot,[fx[45],fx[45]],[1e-24,1e-15],col=120,linestyle=2
		oplot,[fx[71],fx[71]],[1e-24,1e-15],col=120,linestyle=2
		oplot,[fx[72],fx[72]],[1e-24,1e-15],col=120,linestyle=2
		al_legend,'LFR-Low',position=[2.5,1e-15*0.8],box=0,textcolors=120
		al_legend,'LFR-High',position=[25,1e-15*0.8],box=0,textcolors=120
		al_legend,'HFR-Low',position=[350,1e-15*0.8],box=0,textcolors=120
		al_legend,'HFR-High',position=[6000,1e-15*0.8],box=0,textcolors=120
		loadct,1
		for i=0,n_elements(x1[*,0])-1 do begin
			oplot,fx[0:44],x1[i,0:44],col=220
			oplot,fx[72:-1],x1[i,72:-1],col=220
		endfor
		loadct,5
		x1med=fltarr(n_elements(x1[0,*]))
		for i=0,n_elements(x1med)-1 do x1med[i]=median(x1[*,i])
		oplot,fx[0:44],x1med[0:44],col=100
		oplot,fx[72:-1],x1med[72:-1],col=100
		loadct,0
	endps,filename='pdf/Juno-Waves_1percent_level_calibrated_'+data_type_treat+version_name





	plotps,'pdf/estimated_error_calibration_gain_'+data_type_treat+version_name
		restore,'gain_juno_casvg_'+data_type_treat+version_name+'.sav',/verb
		fx=freq
		g=gain

		errorPLUS=g/g50
		errorMINUS=g/g01
		errorPLUS[45:71]=0
		errormean=mean(errorplus)
		errormax=max(errorplus)
		restore,'gain_final_'+data_type_treat+version_name+'.sav',/verb
		errorPLUS[45:71]=3.
		errorMINUS[45:71]=1./3.
		
		
		loadct,13
		plot_oo,fx,g01/g50,line=2.,yrange=[min([errorplus,errorminus]),max([errorminus,errorplus])],/xstyle,ytit='Estimated error on conversion factor',xtit='Frequency (kHz)',/nodata,charsize=1.75
		oplot,[1.,1.e5],[1,1],line=1
		oplot,fx,errorminus,line=1,thick=5		
		oplot,fx,errorplus,line=1,thick=5,color=250
		loadct,53
		    oplot,[fx[27],fx[27]],[0.1,10],col=120,linestyle=2
		    oplot,[fx[45],fx[45]],[0.1,10],col=120,linestyle=2
		    oplot,[fx[71],fx[71]],[0.1,10],col=120,linestyle=2
		    oplot,[fx[72],fx[72]],[0.1,10],col=120,linestyle=2
		    al_legend,'LFR-Low',position=[2.5,10*0.9],box=0,textcolors=120
		    al_legend,'LFR-High',position=[25,10*0.9],box=0,textcolors=120
		    al_legend,'HFR-Low',position=[350,10*0.9],box=0,textcolors=120
		    al_legend,'HFR-High',position=[6000,10*0.9],box=0,textcolors=120
		loadct,0
	endps,filename='pdf/estimated_error_calibration_gain_'+data_type_treat+version_name


;print,strtrim(mean(g01[indmin:indmax]/gain[indmin:indmax]),2)+' pm' +strtrim(stddev(g01[indmin:indmax]/gain[indmin:indmax]),2)+'  ['+strtrim(min(g01[indmin:indmax]/gain[indmin:indmax]),2)+','+strtrim(max(g01[indmin:indmax]/gain[indmin:indmax]),2)+']  ;  '+strtrim(mean(g50[indmin:indmax]/gain[indmin:indmax]),2)+' pm' +strtrim(stddev(g50[indmin:indmax]/gain[indmin:indmax]),2)+'  ['+strtrim(min(g50[indmin:indmax]/gain[indmin:indmax]),2)+','+strtrim(max(g50[indmin:indmax]/gain[indmin:indmax]),2)+']'


endif



END
