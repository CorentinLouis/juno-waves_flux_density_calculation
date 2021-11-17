pro ALL_TIMESERIES,names,yyyydddb=yyyydddb,yyyyddde=yyyyddde,cmin=cmin,cmax=cmax,delta_t=delta_t,titre=titre,$
                  INVERSEMASK=INVERSEMASK,LINEAR=LINEAR,dB=dB,NOBACK_SUBTRACT=NOBACK_SUBTRACT,RAW=RAW,CALIBRATED=CALIBRATED,$
                  version=version
  ; creates saves sets of intensity time series for all the dates for 2016100 to 2018222 without background subtraction with the selected parameters

  ; # INPUT parameters:
  ; # names like ['NK','NK ','NK?','NB ','QP ','QP?', ...] to include or exclude (cf. below) in the time series
  ; # cmin,cmax the canals (in #) (-16 compared to the official Waves channel number)
  ; # here channel#0=1.00100 kHz, corresponding to the official channel#16
  ; # delta_t, in seconds (multiple of 15)
  ; # titre, included in the name of the resulting save set
  
  if ~keyword_set(yyyydddb) then yyyydddb=2016100
  if ~keyword_set(yyyyddde) then yyyyddde=2018064

print,'yyyydddb=',strtrim(yyyydddb,2),' - yyyyddde=',strtrim(yyyyddde,2)
  OUTPUT_PATH='./make_timeseries/'
  cmd = 'mkdir make_timeseries'
  spawn,cmd,resu


  if strlowcase(names) eq 'nkom' then begin
    names = 'NK '
    INVERSEMASK=0
    OUTPUT_PATH=OUTPUT_PATH+'nKOM'
    if ~keyword_set(cmin) then cmin = 26 
    if ~keyword_set(cmax) then cmax = 51
    if ~keyword_Set(titre) then titre='nKOM'
  endif else $
  if strlowcase(names) eq 'bkom' then begin
    names = ['A','E','G','HO','I','N','QP','S'] 
    INVERSEMASK=1
    OUTPUT_PATH=OUTPUT_PATH+'bKOM'
    if ~keyword_set(cmin) then cmin = 14
    if ~keyword_set(cmax) then cmax=44
    if ~keyword_Set(titre) then titre='bKOM'
  endif else $
  if strlowcase(names) eq 'dam' then begin
    names = ['A','E','G','HO','I','N','QP','S'] 
    INVERSEMASK=1
    OUTPUT_PATH=OUTPUT_PATH+'DAM'
    if ~keyword_set(cmin) then cmin = 72
    if ~keyword_set(cmax) then cmax=109
    if ~keyword_Set(titre) then titre='DAM'
  endif else $
  if strlowcase(names) eq 'qp' then begin
    names = 'QP '
    INVERSEMASK=0
    OUTPUT_PATH=OUTPUT_PATH+'QP'
    if ~keyword_set(cmin) then cmin = 0
    if ~keyword_set(cmax) then cmax = 44
    if ~keyword_Set(titre) then titre='QP'
  endif
  if strlowcase(names) eq '' or strlowcase(names) eq 'all' then begin
    names=''
    INVERSEMASK=0
    NOMASK=1
    OUTPUT_PATH=OUTPUT_PATH+'ALL'
    if ~keyword_set(cmin) then cmin = 0
    if ~keyword_set(cmax) then cmax = 109
    if ~keyword_Set(titre) then titre='ALL'
  endif
  ; INPUT keywords:
  ; if /INVERSEMASK, mask is replaced by 1-mask (names are excluded)



  if keyword_set(CALIBRATED) then OUTPUT_PATH=OUTPUT_PATH+'_linear_calibrated_noback_subtract' $
    else begin
      if keyword_Set(linear) then OUTPUT_PATH=OUTPUT_PATH+'_linear' else if keyword_set(dB) then OUTPUT_PATH=OUTPUT_PATH+'_db'
      if keyword_Set(raw) then OUTPUT_PATH=OUTPUT_PATH+'_raw_noback_subtract' $
         else if keyword_set(NOBACK_SUBTRACT) then OUTPUT_PATH=OUTPUT_PATH+'_noback_subtract' else OUTPUT_PATH=OUTPUT_PATH+'_back_subtract'
    endelse



  OUTPUT_PATH=OUTPUT_PATH+'/'

  ; OUTPUT quantities:
  ; timeseries, of intensities corresponding to input parameters and keywords
  ; time, time ramp corresponding to timeseries (doy 2016)

if names ne '' then MAKE_TIMESERIES, yyyydddb,yyyyddde,names,cmin+16,cmax+16,delta_t,titre,timeseries,time,LINEAR=LINEAR,dB=dB,NOBACK_SUBTRACT=NOBACK_SUBTRACT,RAW=RAW,OUTPUT_PATH=OUTPUT_PATH,INVERSEMASK=INVERSEMASK,CALIBRATED=CALIBRATED,version=version
if names eq '' then MAKE_TIMESERIES, yyyydddb,yyyyddde,names,cmin+16,cmax+16,delta_t,titre,timeseries,time,LINEAR=LINEAR,dB=dB,NOBACK_SUBTRACT=NOBACK_SUBTRACT,RAW=RAW,OUTPUT_PATH=OUTPUT_PATH,NOMASK=NOMASK,CALIBRATED=CALIBRATED,version=version
end
;2016100
;2018064