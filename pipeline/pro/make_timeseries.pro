pro MAKE_TIMESERIES,yyyydddb,yyyyddde,names,cmin,cmax,delta_t,titre, timeseries,time, $
                    LINEAR=LINEAR,dB=dB, RAW=RAW, NOBACK_SUBTRACT=NOBACK_SUBTRACT, CALIBRATED=CALIBRATED,$
                    NOMASK=NOMASK, INVERSEMASK=INVERSEMASK, $
                    OUTPUT_PATH=OUTPUT_PATH,version=version

  ; creates save sets of intensity time series for selected parameters

  ; INPUT parameters: 
  ; yyyydddb,yyyyddde dates begin,end included
  ; names like ['NK','NK ','NK?','NB ','QP ','QP?', ...] to include or exclude (cf. below) in the time series
  ; cmin,cmax: frequency channel numbers
  ; delta_t in seconds (multiple of 15)
  ; titre, included in the name of the resulting save set

  ; INPUT keywords:
  ; if /LINEAR, integration /t,f is done with linear flux values
  ; if /RAW, integration /t,f is done with raw flux values (rlin, rdb)
  ; if /NOBACK_SUBTRACT, data without background subtraction is used (zlin, zdb)
  ; if neither /RAW nor /NOBACK_SUBTRACT, then data with background subtraction is used (zlin2, zdb2)
  ; if /NOMASK, mask is not taken into account
  ; if /INVERSEMASK, mask is replaced by 1-mask (names are excluded)
  if ~keyword_set(OUTPUT_PATH) then output_path='make_timeseries/'

  ; OUTPUT quantities: 
  ; timeseries, of intensities corresponding to input parameters and keywords
  ; time, time ramp corresponding to timeseries (doy 2016)


if keyword_set(version) then version_name='_v'+string(format='(I02)',version) else version_name='_v01'
if version eq 1 then version_dailysaveset = '' else version_dailysaveset = version_name

print,"from: ",strtrim(yyyydddb), ' to: ',strtrim(yyyyddde)

nd=aj_t16(yyyyddde)-aj_t16(yyyydddb)+1.
ntime=5760L*nd
nchannel=cmax-cmin+1
timeseries=fltarr(nchannel,ntime)
time=dindgen(ntime)/5760.d0+aj_t16(yyyydddb)	; t16
nn=n_elements(names)
if names ne '' then begin
  restore,'catalog_emissions_for_calibration.sav'
  nem=n_elements(em)
  b=bytarr(nem)
  em0=em
endif else begin
  nem = 0
  b = []
  em0 = []
endelse


if keyword_set(CALIBRATED) then begin
    print, 'ZLINCAL will be used'
endif else begin
  if keyword_set(linear) then begin
    if keyword_set(raw) then print, 'RLIN will be used'
    if keyword_set(NOBACK_SUBTRACT) then print, 'ZLIN will be used'
    if ~keyword_set(raw) and ~keyword_set(NOBACK_SUBTRACT) then print, 'ZLIN2 will be used'
  endif else if keyword_set(dB) then begin
    if keyword_set(raw) then print, 'RDB will be used'
    if keyword_set(NOBACK_SUBTRACT) then print, 'ZDB will be used'
    if ~keyword_set(raw) and ~keyword_set(NOBACK_SUBTRACT) then print, 'ZDB2 will be used'
  endif
endelse

for k=0L,nd-1L do begin
  yyyyddd=t16_aj(aj_t16(yyyydddb)+k)
  w=where(em0.yyyyddd eq yyyyddd,count)
  if count eq 0 then goto,suite
  if keyword_set(CALIBRATED) then input_path='../../data_n2/spdyn_sav_data/'+strmid(strtrim(long(yyyyddd),2),0,4)+'/' $
    else input_path='../../data_n1/spdyn_sav_data/'+strmid(strtrim(long(yyyyddd),2),0,4)+'/'
  restore,input_path+strtrim(long(yyyyddd),2)+'_spdyn'+version_dailysaveset+'.sav'
  print,'### '+input_path+strtrim(long(yyyyddd),2)+'_spdyn.sav file restored ###'
  nt=n_elements(t)
  nf=n_elements(f) 
  print,'day = '+strtrim(yyyyddd,2)+' ; start = ',strtrim(k*nt,2)+' ; end = '+strtrim((k+1)*nt-1,2)+' ; limit = '+strtrim(ntime-1,2)

  if keyword_set(CALIBRATED) then data=zlincal else begin
    if keyword_set(RAW) then if keyword_set(LINEAR) then data=rlin else data=10^(rdb/10.)
   if keyword_set(NOBACK_SUBTRACT) then if keyword_set(LINEAR) then data=zlin else data=10^((zdb*10.+min(rdb))/10.)
   if not(keyword_set(RAW)) and not(keyword_set(NOBACK_SUBTRACT)) then if keyword_set(LINEAR) then data=10^(zlin2/10.) else data=zdb2
  endelse


  mask=fltarr(nt,nf) ;# creating the empty mask
  b[*]=0    ;# emptying the table of emissions corresponding to the criteria (names and day)
  for i=0,nn-1 do $
    for j=0,nem-1 do $
    if strcmp(em(j).name, strupcase(names(i)), strlen(names(i))) and em(j).yyyyddd eq yyyyddd then b(j)=1
  
  wb=where(b eq 1)
  if wb(0) ne -1 then begin  	; fill the mask
    em2=em(wb)
    nem2=n_elements(em2)
    for i=0,nem2-1 do begin
      w=where(em2(i).x ne -1)
      npts=n_elements(w)
      x=[em2(i).x(0:npts-1),em2(i).x(0)]*5760./24.
      if keyword_Set(CALIBRATED) then y=[em2(i).y(0:npts-1),em2(i).y(0)]-16 $
        else y=[em2(i).y(0:npts-1),em2(i).y(0)]
      ind=polyfillv(x,y,nt,nf)
      mask(ind)=1.
    endfor
  endif
  if keyword_set(INVERSEMASK) then mask=1.-mask
  if keyword_set(NOMASK) then mask=fltarr(nt,nf)+1.

  timeseries(*,k*nt:(k+1)*nt-1)=transpose(data(*,cmin:cmax)*mask(*,cmin:cmax))

suite:
endfor

;w0=where(timeseries eq 0.)
;timeseries[w0]=-1

ndt=long(ntime*15/delta_t)
timeseries=rebin(timeseries(*,0:ndt*delta_t/15-1),nchannel,ndt)
time=rebin(time(0:ndt*delta_t/15-1),ndt)
; #To erase 0 points which are days without data
;if keyword_set (nozero) then $
;  for ichannel=0,nchannel do begin
;    w=where(timeseries(ichannel,*) gt 0)
;    timeseries(ichannel,*)=timeseries(w)
;    time=time(w)
;nombre d'éléments sur un delta_t d'intégration
;ntdt=nt*delta_t/86400.
;s=n_elements(timeseries)-ntdt+1
;meants=fltarr(s)
;endfor

;for i=0,s-1 do meants(i)=mean(timeseries(i:i+ntdt-1))


frequencies=f[cmin:cmax]
;#saving results
;# saving all channels
print,'Saving results for all channels'
if keyword_set(CALIBRATED) then begin
    save,timeseries,time,frequencies,filename=output_path+strtrim(titre,2) +'_timeseries_d'+strtrim(long(delta_t),2)+ '_channels_'+strtrim(long(cmin),2)+'-'+strtrim(long(cmax),2)+'_zlincal'+version_name+'.sav'
endif else begin
  if keyword_set(RAW) then begin
    if keyword_set (LINEAR) then begin
      save,timeseries,time,frequencies,filename=output_path+strtrim(titre,2)+'_timeseries_d'+strtrim(long(delta_t),2)+ '_channels_'+strtrim(long(cmin),2)+'-'+strtrim(long(cmax),2)+'_rlin'+version_name+'.sav'
      ;p=plot(time,timeseries,/xsty,/ysty,title=titre,ytitle='Intensity(linear)',xtitle='DOY 2016',title='Raw')
      ;p.Save,'./make_timeseries/'+strtrim(strmid(titre,0,4))+'_timeseries_d'+strtrim(long(delta_t),2)+ '_rlin.png'
    endif else begin
      save,timeseries,time,frequencies,filename=output_path+strtrim(titre,2) +'_timeseries_d'+strtrim(long(delta_t),2)+ '_channels_'+strtrim(long(cmin),2)+'-'+strtrim(long(cmax),2)+'_rdb'+version_name+'.sav'
      ;p=plot(time,timeseries,/xsty,/ysty,title=titre,ytitle='Intensity(dB)',xtitle='DOY 2016',title='Raw')
      ;p.Save,'./make_timeseries/'+strtrim(strmid(titre,0,4)) +'_timeseries_d'+strtrim(long(delta_t),2)+ '_rdb.png'
    endelse
  endif
  if not(keyword_set(RAW)) then begin
    if keyword_set(NOBACK_SUBTRACT) then begin
      if keyword_set (LINEAR) then begin
       save,timeseries,time,frequencies,filename=output_path+strtrim(titre,2) +'_timeseries_d'+strtrim(long(delta_t),2)+ '_channels_'+strtrim(long(cmin),2)+'-'+strtrim(long(cmax),2)+'_zlin'+version_name+'.sav'
       ;p=plot(time,timeseries,/xsty,/ysty,title=titre,ytitle='Intensity(linear)',xtitle='DOY 2016',title='FFT-filtered')
       ;p.Save,'./make_timeseries/'+strtrim(strmid(titre,0,4)) +'_timeseries_d'+strtrim(long(delta_t),2)+ '_zlin.png'
      endif else begin
        save,timeseries,time,frequencies,filename=output_path+ strtrim(titre,2) +'_timeseries_d'+strtrim(long(delta_t),2)+ '_channels_'+strtrim(long(cmin),2)+'-'+strtrim(long(cmax),2)+'_zdb'+version_name+'.sav'
        ;p=plot(time,timeseries,/xsty,/ysty,title=titre,ytitle='Intensity(dB)',xtitle='DOY 2016',title='FFT-filtered')
       ;p.Save,'./make_timeseries/'+strtrim(strmid(titre,0,4)) +'_timeseries_d'+strtrim(long(delta_t),2)+ '_zdb.png'
      endelse
    endif else begin
      if keyword_set (LINEAR) then begin
        save,timeseries,time,frequencies,filename=output_path+strtrim(strmid(titre,0,4)) +'_timeseries_d'+strtrim(long(delta_t),2)+ '_channels_'+strtrim(long(cmin),2)+'-'+strtrim(long(cmax),2)+'_zlin2'+version_name+'.sav'
      endif else begin
        save,timeseries,time,frequencies,filename=output_path+ strtrim(strmid(titre,0,4)) +'_timeseries_d'+strtrim(long(delta_t),2)+ '_channels_'+strtrim(long(cmin),2)+'-'+strtrim(long(cmax),2)+'_zdb2'+version_name+'.sav'
      endelse
    endelse
  endif
endelse
print,'Results saved'
;# saving each channel
;timeseries_tmp=timeseries
;delvar,timeseries
;for ichannel=cmin,cmax do begin
;  print,'Saving results for channel #'+strtrim(ichannel,2)
;  timeseries=reform(timeseries_tmp[ichannel-cmin,*])
;  frequency=f(ichannel)
;  ;saving results
;  if keyword_set(RAW) then begin
;    if keyword_set (LINEAR) then begin
;      save,timeseries,time,frequency,filename=OUTPUT_PATH+strtrim(strmid(titre,0,4)) +'_timeseries_d'+strtrim(long(delta_t),2)+ '_channel_'+strtrim(long(ichannel),2)+'_rlin.sav'
;      ;p=plot(time,timeseries,/xsty,/ysty,title=titre,ytitle='Intensity(linear)',xtitle='DOY 2016',title='Raw')
;      ;p.Save,'./make_timeseries/'+strtrim(strmid(titre,0,4))+'_timeseries_d'+strtrim(long(delta_t),2)+ '_rlin.png'
;    endif else begin
;      save,timeseries,time,frequency,filename=OUTPUT_PATH+strtrim(strmid(titre,0,4)) +'_timeseries_d'+strtrim(long(delta_t),2)+ '_channel_'+strtrim(long(ichannel),2)+'_rdb.sav'
;      ;p=plot(time,timeseries,/xsty,/ysty,title=titre,ytitle='Intensity(dB)',xtitle='DOY 2016',title='Raw')
;      ;p.Save,'./make_timeseries/'+strtrim(strmid(titre,0,4)) +'_timeseries_d'+strtrim(long(delta_t),2)+ '_rdb.png'
;    endelse
;  endif
;  if not(keyword_set(RAW)) then begin
;    if keyword_set(NOBACK_SUBTRACT) then begin
;      if keyword_set (LINEAR) then begin
;        save,timeseries,time,frequency,filename=OUTPUT_PATH+strtrim(strmid(titre,0,4)) +'_timeseries_d'+strtrim(long(delta_t),2)+ '_channel_'+strtrim(long(ichannel),2)+'_zlin.sav'
;        ;p=plot(time,timeseries,/xsty,/ysty,title=titre,ytitle='Intensity(linear)',xtitle='DOY 2016',title='FFT-filtered')
;        ;p.Save,'./make_timeseries/'+strtrim(strmid(titre,0,4)) +'_timeseries_d'+strtrim(long(delta_t),2)+ '_zlin.png'
;      endif else begin
;        save,timeseries,time,frequency,filename=OUTPUT_PATH+strtrim(strmid(titre,0,4)) +'_timeseries_d'+strtrim(long(delta_t),2)+ '_channel_'+strtrim(long(ichannel),2)+'_zdb.sav'
;        ;p=plot(time,timeseries,/xsty,/ysty,title=titre,ytitle='Intensity(dB)',xtitle='DOY 2016',title='FFT-filtered')
;        ;p.Save,'./make_timeseries/'+strtrim(strmid(titre,0,4)) +'_timeseries_d'+strtrim(long(delta_t),2)+ '_zdb.png'
;      endelse
;    endif else begin
;      if keyword_set (LINEAR) then begin
;        save,timeseries,time,frequency,filename=OUTPUT_PATH+strtrim(strmid(titre,0,4)) +'_timeseries_d'+strtrim(long(delta_t),2)+ '_channel_'+strtrim(long(ichannel),2)+'_zlin2.sav'
;      endif else begin
;        save,timeseries,time,frequency,filename=OUTPUT_PATH+strtrim(strmid(titre,0,4)) +'_timeseries_d'+strtrim(long(delta_t),2)+ '_channel_'+strtrim(long(ichannel),2)+'_zdb2.sav'
;      endelse
;    endelse
;  endif
;print,'Channel #'+strtrim(ichannel,2)+' saved'
;endfor
return
end