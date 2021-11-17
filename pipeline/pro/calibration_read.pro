pro CALIBRATION_READ,cmin,cmax,delta_t,$
                    QP=QP,NKOM=NKOM,BKOM=BKOM,DAM=DAM,ALL=ALL,$
                    LINEAR=LINEAR,dB=dB,RAW=RAW,NOBACK_SUBTRACT=NOBACK_SUBTRACT,CALIBRATED=CALIBRATED,version=version

 ; INPUT
 ; delta_t, in seconds (multiple of 15)

 ; INPUT keywords:
 ; if /NKOM, all the 'nKOM' timeseries saves will be loaded
 ; if /BKOM, all the 'bKOM' timeseries saves will be load
 ; if /BKOM, all the 'HOM' timeseries saves will be load

if keyword_set(version) then version_name='_v'+string(format='(I02)',version) else version_name='_v01'

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
; print,'Calibration is only done with LINEAR data. Please set the /LINEAR keyword'
endelse

if keyword_set(bKOM) then emtype='bKOM' $
  else if  keyword_set(nKOM) then emtype='nKOM' $
  else if  keyword_set(DAM) then emtype='DAM' $
  else if  keyword_set(QP) then emtype='QP' $
  else if keyword_set(ALL) then emtype='ALL'


restore, 'Ephemerides_Juno.sav' 
print,'Ephem restored'
e=ephemerides
delvar,ephemerides
t=JULDAY(1,e.day,e.yy,e.hr,e.min,e.sec) - JULDAY(1,0,2016,0,0,0)

print,emtype+' selected'
print,emtype+' writing...'
print,'restoring timeseries '+emtype+'_'+type+'/'+emtype+'_timeseries_d'+strtrim(long(delta_t),2)+'_channels_'+strtrim(cmin,2)+'-'+strtrim(cmax,2)+'_'+dataype+'.sav'
restore,'../make_timeseries/'+emtype+'_'+type+'/'+emtype+'_timeseries_d'+strtrim(long(delta_t),2)+'_channels_'+strtrim(cmin,2)+'-'+strtrim(cmax,2)+'_'+dataype+version_name+'.sav',/verb
print,'Timeseries restored'
R_time=interpol(e.dist_rj,t,time)
lat_time=interpol(e.oblat,t,time)
mlat_time=interpol(e.mlat,t,time)
s=timeseries
print,'Writing complete'

;# saving results
print,'Saving...'
cmd='mkdir ../calibration_timeseries'
spawn,cmd,resu
print,"saving "+'../calibration_timeseries/'+ emtype+'_calibration_timeseries_'+type+'_d'+strtrim(long(delta_t),2)+'_'+dataype+version_name+'.sav'
save,s,time,R_time,lat_time,mlat_time,frequencies,filename='../calibration_timeseries/'+ emtype+'_calibration_timeseries_'+type+'_d'+strtrim(long(delta_t),2)+'_'+dataype+version_name+'.sav'
print,'Save complete'
end