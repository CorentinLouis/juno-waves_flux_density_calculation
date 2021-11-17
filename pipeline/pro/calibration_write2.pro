pro CALIBRATION_WRITE2,delta_t,period=period, Rmin=Rmin, Rmax=Rmax,MLatMin=MLatMin,MLatMax=MLatMax, $
                        norm1AU=norm1AU,$
                        QP=QP,NKOM=NKOM,BKOM=BKOM,DAM=DAM,ALL=ALL,$
                        LINEAR=LINEAR,dB=dB,RAW=RAW,NOBACK_SUBTRACT=NOBACK_SUBTRACT,CALIBRATED=CALIBRATED,version=version

    ; INPUT:
	; delta_t, in seconds (multiple of 15)

	; INPUT keywords:
 	; if /NKOM, all the 'nKOM' timeseries saves will be loaded
	; if /BKOM, all the 'bKOM' timeseries saves will be load
	; if /DAM, all the 'DAM' timeseries saves will be load

    if keyword_set(version) then version_name='_v'+string(format='(I02)',version) else version_name='_v01'

    if ~keyword_set(MLatMax) and ~keyword_set(MLatMin) then begin
        MLatMax=90.
        print,'No MLat limit has been set. Thus MLatMax=90 deg.'
    endif
    if ~keyword_set(Rmax) and ~keyword_set(Rmin) then stop,'You have to set a Distance Limit (Rmax=XX or Rmin=XX)'


    if keyword_set(Rmax) and keyword_set(MLatMax) then tit_limit='_Rmax'+strtrim(long(Rmax),2)+'_MLatmax'+strtrim(long(MLatMax),2) $
        else if keyword_set(Rmax) and keyword_set(MLatMin) then tit_limit='_Rmax'+strtrim(long(Rmax),2)+'_MLatmin'+strtrim(long(MLatMin),2) $
        else if keyword_set(Rmin) and keyword_set(MLatMin) then tit_limit='_Rmin'+strtrim(long(Rmin),2)+'_MLatmin'+strtrim(long(MLatMin),2) $
        else if keyword_set(Rmin) and keyword_set(MLatMax) then tit_limit='_Rmin'+strtrim(long(Rmin),2)+'_MLatmax'+strtrim(long(MLatMax),2)

    if keyword_set(CALIBRATED) then begin
        type='linear_calibrated_back_subtract'
        dataype='zlincal'
        print,'zlincal will be used'
    endif else begin
        if keyword_set(LINEAR) then  begin
            type='linear'
            if keyword_set(RAW) then type=type+'_raw'
            if keyword_set(NOBACK_SUBTRACT) then begin
                type=type+'_noback_subtract'
                dataype='zlin'
                if keyword_set(RAW) then dataype='rlin'
            endif else begin
                type=type+'_back_subtract'
                dataype='zlinnb'
                if keyword_set(RAW) then dataype='rlinnb'
            endelse
        endif else if keyword_set(dB) then begin
            type='db'
            if keyword_set(RAW) then type=type+'_raw'
            if keyword_set(NOBACK_SUBTRACT) then begin
                type=type+'_noback_subtract'
                dataype='zdb'
                if keyword_set(RAW) then dataype='rdb'
            endif else begin
                type=type+'_back_subtract'
                dataype='zdbnb'
                if keyword_set(RAW) then dataype='rdbnb'
            endelse
        endif
        ;print,'Calibration is only done with LINEAR data. Please set the /LINEAR keyword'
    endelse

    if keyword_set(bKOM) then emtype='bKOM' $
      else if  keyword_set(nKOM) then emtype='nKOM' $
      else if  keyword_set(DAM) then emtype='DAM' $
      else if  keyword_set(QP) then emtype='QP' $
      else if  keyword_set(ALL) then emtype='ALL'

 	;restore, '../../stage/frequences.sav'
    ;print, 'frequency file restored'
	;restore, 'Ephemerides_Juno.sav'
    ;print, "Juno's ephemeris file restored"
	;e=ephemerides
	;delvar,ephemerides
	;t=JULDAY(1,e.day,e.yy,e.hr,e.min,e.sec)-JULDAY(1,0,2016,0,0,0)
    ;period=1440l  
    if ~keyword_set(period) then period = 4*595.5 ;# pjup=595.5 #Jupiter's rotation period


	print,emtype+' selected'
	
    print,"restoring "+'../calibration_timeseries/'+ emtype+'_calibration_timeseries_'+type+'_d'+strtrim(long(delta_t),2)+'_'+dataype+'.sav'
    restore,'../calibration_timeseries/'+ emtype+'_calibration_timeseries_'+type+'_d'+strtrim(long(delta_t),2)+'_'+dataype+version_name+'.sav',/verb
    
    print,"Timeseries file restored"
    nchannels=n_elements(s[*,0])

    



    if keyword_set(Rmin) and keyword_set(MLatMax) then wint=where(R_time ge Rmin and abs(mlat_time) le MLatMax) $
            else if keyword_set(Rmin) and keyword_set(MLatMin) then wint=where(R_time ge Rmin and abs(mlat_time) ge MLatMin) $
            else if keyword_set(Rmax) and keyword_set(MLatMin) then wint=where(R_time le Rmax and abs(mlat_time) ge MLatMin) $
            else if keyword_set(Rmax) and keyword_set(MLatMax) then wlatdist=where(R_time le Rmax and abs(mlat_time) le MLatMax)
        
    ind=0d
    time_to_keep=[]
    while ind ne n_elements(wint)-period do begin
        if wint[ind+period]- wint[ind] eq period then begin
            time_to_keep=[time_to_keep,wint[ind:ind+period-1]]
            ind=ind+period
        endif else ind=ind+1d
    endwhile
    
    nt=long(n_elements(time))
    ninterval=n_elements(time_to_keep)/period
    s2=dblarr(nchannels,ninterval,period)
    time2=dblarr(ninterval,period)
    R_time2=dblarr(ninterval,period)
    mlat_Time2=dblarr(ninterval,period)
    s_tmp=dblarr(nchannels,nt)

    AU=149597870.700d

    if keyword_set(norm1AU) then $
        for i=0,nchannels-1 do begin
          s_tmp[i,where(s[i,*] gt 0)]=(s[i,where(s[i,*] gt 0)]/377.)*(r_time[where(s[i,*] gt 0)]*71492./AU)^2
        endfor
    wfreqplot=where(frequencies eq 6500.00)
    if ~keyword_set(raw) then begin
        plotps,'pdf/timeseries_'+strtrim(long(frequencies[wfreqplot]),2)+'kHz_back-subtracted_normalized'

            wn0=where(s_tmp[wfreqplot,*] gt 0.)  
            plot_io,time,s_tmp[wfreqplot,*],/xsty,/ysty,xtit=TeXtoIDL('Time (Day of Year 2016)'),ytit=TexToIDL('Normalized intensity (W/m^2/Hz)'),tit=TexToIDL('Timeseries @ '+strtrim(frequencies[wfreqplot],2)+' kHz - background subtracted - normalized @ 1 AU'),yrange=[min(s_tmp[wfreqplot[0],wn0]),max(s_tmp[wfreqplot[0],wn0])]
        endps,filename='pdf/timeseries_'+strtrim(long(frequencies[wfreqplot]),2)+'kHz_back-subtracted_normalized'
    endif

    for i=0,ninterval-1 do begin
        s2[*,i,*]=s_tmp[*,time_to_keep[period*i:period*(i+1)-1]]
        time2[i,*]=time[time_to_keep[period*i:period*(i+1)-1]]
        R_time2[i,*]=time[time_to_keep[period*i:period*(i+1)-1]]
        mlat_time2[i,*]=time[time_to_keep[period*i:period*(i+1)-1]]
    endfor

    xmoy=fltarr(ninterval,nchannels)
    xmed=fltarr(ninterval,nchannels)
    x1=fltarr(ninterval,nchannels)

    for iint=0,ninterval-1 do begin
        for ichannel=0,nchannels-1 do begin
            wn0=where(s2[ichannel,iint,*] gt 0.)
            if wn0[0] ne -1 then begin
                xmoy[iint,ichannel]=mean(s2[ichannel,iint,wn0])
                xmed[iint,ichannel]=median(s2[ichannel,iint,wn0])
                x1[iint,ichannel]=dyn_n(s2[ichannel,iint,wn0],0.99)
            endif
        endfor
    endfor

    fx=frequencies
    time=time2
    ;# saving files
    print,'Saving '+'../calibration_data/'+emtype+'_calibration_data_freq_'+type+'_d'+strtrim(long(delta_t),2)+'_'+dataype+tit_limit+'.sav'
    save,time,period,fx,xmoy,xmed,x1,filename='../calibration_data/'+emtype+'_calibration_data_freq_'+type+'_d'+strtrim(long(delta_t),2)+'_'+dataype+tit_limit+version_name+'.sav'

    print,'Save complete'


return


;# grouped data by interval
    nt=long(n_elements(time))
    ;nn=long(nt/(24.*60.))
    nn=long(nt/period)
    s2=dblarr(nchannels,nn,period)
    time2=dblarr(nn,period)
    R_time2=dblarr(nn,period)
    mlat_Time2=dblarr(nn,period)


    
    ;# intensity normalized at 1 AU
    ;# and /377 in order to retrieve z0
    if keyword_set(norm1AU) then $
        for i=0,nchannels-1 do begin
          s_tmp[i,where(s[i,*] gt 0)]=(s[i,where(s[i,*] gt 0)]/377.)*(r_time[where(s[i,*] gt 0)]*71400/1.49e8)^2
        endfor
    wfreqplot=where(frequencies eq 6500.00)
    if ~keyword_set(raw) then begin
        plotps,'timeseries_'+strtrim(long(frequencies[wfreqplot]),2)+'kHz_back-subtracted_normalized'
            plot_io,time,s_tmp[wfreqplot,*],yrange=[1e-30,1e-16],/xsty,/ysty,xtit=TeXtoIDL('Time (Day of Year 2016)'),ytit=TexToIDL('Normalized intensity (W/m^2/Hz)'),tit=TexToIDL('Timeseries @ '+strtrim(frequencies[wfreqplot],2)+' kHz - background subtracted - normalized @ 1 AU')
        endps,filename='timeseries_'+strtrim(long(frequencies[wfreqplot]),2)+'kHz_back-subtracted_normalized'
    endif

    for j=0,nn-1 do begin
        time2(j,*)=time(period*j:period*(j+1)-1)
        r_time2(j,*)=r_time(period*j:period*(j+1)-1)
        mlat_time2(j,*)=mlat_time(period*j:period*(j+1)-1)
        for i=0,nchannels-1 do begin
            s2(i,j,*)=s_tmp(i,period*j:period*(j+1)-1)
        endfor
    endfor

    xmoy=fltarr(nn,nchannels)
    xmed=fltarr(nn,nchannels)
    x1=fltarr(nn,nchannels)
    interval0=intarr(nn)
   

    for iint=0,nn-1 do begin
        
         ;# filtering the timeseries
            ;# e.g. for the calibration of LFR-L/H and HFR-H, R > Rmin=30Rj and abs(mlat) < MLatMax=15°

        if keyword_set(Rmin) and keyword_set(MLatMax) then wlatdist=where(R_time2[iint,*] ge Rmin and abs(mlat_time2[iint,*]) le MLatMax) $
            else if keyword_set(Rmin) and keyword_set(MLatMin) then wlatdist=where(R_time2[iint,*] ge Rmin and abs(mlat_time2[iint,*]) ge MLatMin) $
            else if keyword_set(Rmax) and keyword_set(MLatMin) then wlatdist=where(R_time2[iint,*] le Rmax and abs(mlat_time2[iint,*]) ge MLatMin) $
            else if keyword_set(Rmax) and keyword_set(MLatMax) then wlatdist=where(R_time2[iint,*] le Rmax and abs(mlat_time2[iint,*]) le MLatMax)
        if n_elements(wlatdist) eq period then begin
            interval0[iint]=+1
            for ichannel=0,nchannels-1 do begin
                wn0=where((s2[ichannel,iint,*])[wlatdist] gt 0.)
            ;# computing spectra of mean, 50% and 1% occurrence levels 
            ;# (i.e., at each frequency the flux density exceeded 50% and 1% of the time) 

                if wn0[0] ne -1 then begin
                  ;  if ichannel eq 6 then print,    ((s2[ichannel,iint,*])[wlatdist])[wn0]
                    xmoy[iint,ichannel]=mean(((s2[ichannel,iint,*])[wlatdist])[wn0])
                    xmed[iint,ichannel]=median(((s2[ichannel,iint,*])[wlatdist])[wn0])
                    x1[iint,ichannel]=dyn_n(((s2[ichannel,iint,*])[wlatdist])[wn0],0.99)
                endif
            endfor
        endif
    endfor

;# only save the data where d>30R_j & abs(Mlat)<15 degres
    w=where(interval0 gt 0)
    xmoy=xmoy[w,*]
    xmed=xmed[w,*]
    x1=x1[w,*]

    fx=frequencies
    time=time2[w,*]
 	;# saving files
 	print,'Saving...'
    cmd='mkdir ./calibration_data'
    spawn,cmd,resu
    save,time,period,fx,xmoy,xmed,x1,filename='./calibration_data/'+emtype+'_calibration_data_freq_'+type+'_d'+strtrim(long(delta_t),2)+'_'+dataype+tit_limit+version_name+'.sav'

   	print,'Save complete'
end








