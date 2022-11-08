PRO calibration_process,yyyydddb=yyyydddb,yyyyddde=yyyyddde,norm1AU=norm1AU,reso1sec=reso1sec,linear=linear,dB=dB,$
						create_savefile=create_savefile,version=version

;# Steps to apply the calibration to the data :
;# Use FFT-filtered times series from raw data in linear scale
;# Subtract a background (calculated from dB values) in linear scale
;# Apply calibration (gain) table -> one obtains local, calibrated flux densities
;# Optionally (norm1AU), correct for 1/R2 dependence to normalize the flux density to a constant distance (this should not be done when studying local wave E-field, for example)
;# Optionnaly, ask for 1 seconde resolution (for the sources)

if ~keyword_set(yyyydddb) then yyyydddb=2016100
if ~keyword_set(yyyyddde) then yyyyddde=2019174


if keyword_set(version) then version_name='_v'+string(format='(I02)',version) else version_name=''
if ~keyword_set(version) then $
    channelbeg = 16 $
    else channelbeg = 0
channelend = 125


if ~keyword_set(linear) and ~keyword_set(db) then begin
	linear=1
	print,"Calibration will be made using 'linear' process"
endif
if keyword_Set(linear) then begin
	data_type_treat='lin'
	data_type_folder='linear'
endif else if keyword_Set(dB) then begin
	data_type_treat='db'
	data_type_folder='db'
endif

if keyword_Set(norm1AU) then begin
	restore, 'Ephemerides_Juno.sav'
	e=ephemerides
	tephem=JULDAY(1,e.day,e.yy,e.hr,e.min,e.sec) - JULDAY(1,0,2016,0,0,0)
	print,'# Juno ephemeris restored #'
	tit1AU='_norm1AU'
endif else begin
	tit1AU=''
endelse

if keyword_Set(reso1sec) then tit1sec='_1sec' else tit1sec=''

print,"### creating directory if necessary ###"
cmd = 'mkdir ../../data_n2'
spawn,cmd,resu

year_dir = indgen(long(yyyyddde/1000)-long(yyyydddb/1000)+1,start=long(yyyydddb/1000))
for iyear=0,n_elements(year_dir)-1 do begin
        cmd = 'mkdir ../../data_n2/spdyn_sav_data'+tit1AU+version_name+'/'+strmid(strtrim(year_dir[iyear],2),0,4)
	spawn,cmd,resu
endfor

print,"### restoring background ###"
restore,'../background/background_ALLPJ_z'+data_type_treat+version_name+'.sav'

print,"### restoring calibration gain ###"
restore,'../gain/gain_final_'+data_type_treat+version_name+'.sav'
print,'# Calibration gain restored #'

gain=gain_final
if keyword_set(version) then $
	gain=[fltarr(16-channelbeg)+gain[0],gain]

nd=aj_t16(yyyyddde)-aj_t16(yyyydddb)+1.
for k=0L,nd-1L do begin
	yyyyddd=t16_aj(aj_t16(yyyydddb)+k)

	if keyword_set(create_savefile) then $
		CREATE_SAVEFILE_SPDYN_SURVEY, yyyyddd,version=version

	;# day without data to exclude of the calibration process
	skip = 0
	if ~keyword_set(reso1sec) then begin
		if ~file_test('../../data_n1/spdyn_sav_data/'+strmid(strtrim(long(yyyyddd),2),0,4)+'/'+strtrim(long(yyyyddd),2)+'_spdyn'+version_name+'.sav') then begin
			skip = 1
			goto,suite
		endif else skip = 0
	endif else begin
		if ~file_test('../../data_n1/spdyn_sav_data_1sec/'+strmid(strtrim(long(yyyyddd),2),0,4)+'/'+strtrim(long(yyyyddd),2)+'_spdyn_1sec'+version_name+'.sav') then begin
			skip = 1
			goto,suite
		endif else skip = 0
	endelse


	if ~keyword_set(reso1sec) then begin
		;# restore data FFT Filterd 
		restore,'../../data_n1/spdyn_sav_data/'+strmid(strtrim(long(yyyyddd),2),0,4)+'/'+strtrim(long(yyyyddd),2)+'_spdyn'+version_name+'.sav'
		print, '# Data Day '+strtrim(long(yyyyddd),2)+' restored #'
		t2=aj_t16(yyyyddd+t/24.)
	endif else begin
		restore,'../../data_n1/spdyn_sav_data_1sec/'+strmid(strtrim(long(yyyyddd),2),0,4)+'/'+strtrim(long(yyyyddd),2)+'_spdyn_1sec'+version_name+'.sav'
		print, '# Data Day '+strtrim(long(yyyyddd),2)+' restored #'
		t2=aj_t16(yyyyddd+t/24/60./60.)
	endelse
	

	;# background substraction
	if keyword_set(linear) then $
		zlin_new = zlin[*,channelbeg:channelend] $
	else if keyword_set(db) then $
		zlin_new = 10^((zdb[*,channelbeg:channelend]*10.+min(rdb))/10.)

	SUBTRACT_BACKGROUND,background,sigma,0,zlin_new

	print,'# Background subtracted #'
	f=f[channelbeg:channelend]
	nchannels=n_elements(f)

	;# applying calibration gain
	for i=0,nchannels-1 do begin
		zlin_new[*,i]=zlin_new[*,i]*gain[i]
	endfor
	;for i=0,44 do begin
	;	zlin_new[*,i]=zlin_new[*,i]*gain[i]
	;endfor
	;for i=72,nchannels-1 do begin
	;	zlin_new[*,i]=zlin_new[*,i]*gain[i]
	;endfor
	print,'# Calibration gain applied #'

	;# normalization @ 1 AU and retrieve physical Units (W/m^2/Hz)
	if keyword_Set(norm1AU) then begin
		time_new=interpol(tephem,tephem,t2)
		R_time=interpol(e.dist_rj,tephem,t2) 
		lat_time=interpol(e.oblat,tephem,t2) 
		mlat_time=interpol(e.mlat,tephem,t2)
	
		for i=0,nchannels-1 do $
		    zlin_new[where(zlin_new[*,i] gt 0),i]=(zlin_new[where(zlin_new[*,i] gt 0),i]/377.)*(r_time[where(zlin_new[*,i] gt 0)]*71400/1.49e8)^2
		print,'# Data normalized to 1 AU #'
	endif

	zlincal=zlin_new



	if keyword_Set(norm1AU) then tit1AU='_norm1AU' else tit1AU=''
	if keyword_Set(reso1sec) then tit1sec='_1sec' else tit1sec=''

	filename = '../../data_n2/spdyn_sav_data'+tit1AU+version_name+'/'+strmid(strtrim(long(yyyyddd),2),0,4)+'/'+strtrim(long(yyyyddd),2)+'_spdyn'+tit1AU+tit1sec+'_'+data_type_treat+version_name+'.sav'

	
	save,filename=filename,yyyyddd,t,f,zlincal
	print,'# Day '+strtrim(long(yyyyddd),2)+' calibrated and saved #'


	cmd = "python make_cdf_juno_waves_calibrated.py -i "+filename+" -data_version "+strtrim(version,2);#+" -v"
	spawn,cmd,resu
	print,resu

	

suite:	if skip eq 1 then print,'# No data for Day '+strtrim(long(yyyyddd),2) +' in version '+version_name+' #' 
endfor
END
