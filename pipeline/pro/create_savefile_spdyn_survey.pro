PRO read_jwaves_cal_data,yyyydddb,yyyyddde,data,head,mode,version=version,root_path=root_path,file=file

if ~keyword_set(version) then version='01' $
  else version=string(format='(I02)',version)
print,"Working with the PDS version "+version+" files"
year = strmid(strtrim(yyyydddb,1),0,4)
doy  = strmid(strtrim(yyyydddb,1),4,3)

if ~keyword_set(file) then begin
  if keyword_set(root_path) then begin
    file = file_search(root_path,'WAV_'+year+doy+'*E_V'+version+'.CSV')
  endif else begin
    root_path='../../data_PDS/WAVES_SURVEY/'
    file = file_search(root_path+year+'*/','WAV_'+year+doy+'*E_V'+version+'.CSV')
    if file[0] eq '' then file = file_search(root_path+strtrim(long(year)-1,2)+'*/','WAV_'+year+doy+'*E_V'+version+'.CSV')
  endelse
endif

if file[0] eq '' then print,"No file for day "+strtrim(yyyydddb,2)
if file[0] eq '' then return

for i=0L,n_elements(file)-1L do begin
  filei=file(i)
  read_pds_jnowav_survey_e,filei,datai,headi,modei
  if i eq 0 then begin
    data=datai & mode=modei & head=headi
  endif else begin
    data = [data,datai]
    mode = [mode,modei]
    w = where((headi.frequency ne head.frequency) or (headi.bandwidth ne head.bandwidth) or $
	(headi.channelID ne head.channelID) or headi.title ne head.title)
    if w[0] ne -1 then stop,'Header change between consecutive files'
  endelse
endfor

n = aj_t97(yyyyddde)-aj_t97(yyyydddb)
if n ge 1 then begin
  for j=0L,n-1L do begin
    yyyyddd = long(t97_aj(aj_t97(yyyydddb)+j))
    year = strmid(strtrim(yyyyddd,1),0,4)
    doy  = strmid(strtrim(yyyyddd,1),4,3)
    file = file_search(root_path,'WAV_'+year+doy+'*E_V'+version+'.CSV') 
    if file[0] ne '' then begin
      for i=0L,n_elements(file)-1L do begin
        filei=file(i)
        read_pds_jnowav_survey_e,filei,datai,headi,modei
        data = [data,datai]
        mode = [mode,modei]
        w = where((headi.frequency ne head.frequency) or (headi.bandwidth ne head.bandwidth) or $
	  (headi.channelID ne head.channelID) or headi.title ne head.title)
        if w[0] ne -1 then stop,'Header change between consecutive files'
      endfor
    endif
  endfor
endif

end

;--------------------------------------------------------------------------------------------------------------------------------------------------------------------
  pro CREATE_SAVEFILE_SPDYN_SURVEY, yyyydddb, save=save, restore=restore, col=col, plot_spdyn=plot_spdyn,rebin15sec=rebin15sec,plot_filtering=plot_filtering,$
                                    filetest=filetest,NHARMLIN=NHARMLIN,NHARMDB=NHARMDB,dharm=dharm,version=version
;--------------------------------------------------------------------------------------------------------------------------------------------------------------------
; /SAVE or /RESTORE => save/restore data for the day
; col=0 => N/B	col=13 or 39 => Rainbow
; /plot => summary plots

; # Initialisation #
!path=!path+':./IDL_pro_util'
!path='/Users/serpe/Volumes/kronos/juno/stage/IDL_pro_util:'+!path

if ~keyword_set(rebin15sec) then deltat='_1sec' else deltat=''
if ~keyword_set(nharm) then nharm=8. else nharm=nharm

if ~keyword_set(dharm) then dharm=0.15 else dharm=dharm
if ~keyword_set(NHARMLIN) then NHARMLIN=8.
if ~keyword_set(NHARMDB) then NHARMDB=8.
han=1.0
if not(keyword_set(col)) then col=0

; # Reading raw data #
print,"Day "+strtrim(yyyydddb,2)+" is beeing processed"
if keyword_set(restore) then restore,'./data_sav/'+strtrim(long(yyyydddb),2)+'.sav' else read_jwaves_cal_data,yyyydddb,yyyydddb,data,head,mode,version=version
if keyword_set(save) then save,data,head,mode,filename='./data_sav/'+strtrim(long(yyyydddb),2)+'.sav'
if n_elements(data) eq 0 then return

;# Reading data of day before if first element of data.flux equal to 0 #
;# solve interpolation issue when first element of flux is equal to 0 (i.e. avoid extrapolation issue) #
if min(data[0].flux) eq 0 then begin
  yyyyddd_tmp=yyyydddb
  if strmid((strtrim(yyyydddb-1,2)),4,3) ne '000' then $
    file_day_before_name = '../../data_n1/spdyn_sav_data'+deltat+'/'+strmid(strtrim(yyyydddb,2),0,4)+'/'+strtrim(long(yyyydddb-1),2)+'_spdyn'+deltat+'.sav' $
    else file_day_before_name = '../../data_n1/spdyn_sav_data'+deltat+'/'+strtrim(long(strmid(strtrim(yyyydddb,2),0,4))-1,2)+'/'+strtrim(amj_aj(strtrim(fix(yyyydddb/1000.)-1)*10000l+12*100l+31),2)+'_spdyn'+deltat+'.sav'
  file_day_before = file_search(file_day_before_name)
  if file_day_before[0] eq '' then begin
    print,"No data on previous day: no interpolation.  First element goes to 0"
    x=transpose(data.flux) > 0.
  endif else begin
    print,file_day_before_name+" is beeing restored for interpolation of the first element"
    restore,file_day_before_name,/verb
    tdaybefore=t-24.
    x=transpose(data.flux) > 0.
    x=[rlin[-1,*],x]
    delvar,yyyyddd,rdb,rlin,zlin,zdb,t,f
    yyyydddb=yyyyddd_tmp
  endelse
endif else $
  x=transpose(data.flux) > 0.

; # Constructing frequency and time tables
f=head.frequency/1.e3 ; kHz
nf=n_elements(f)
t0_second=float(strmid(data[0].time_iso,15,6))
t0_hour=float(strmid(data[0].time_iso,9,2))   
t0_minute=float(strmid(data[0].time_iso,12,2))
t0=t0_hour*60.*60.+t0_minute*60.+t0_second
t=data.TIME_SCLK-data(0).TIME_SCLK+t0
nt=n_elements(t)

if keyword_set(tdaybefore) then t=[tdaybefore[-1],t]

; # removing first elements, if they still are equal to 0
if min(x[0,*]) le 0 then begin
  i=0
  w0=[]
  while min(x[i,*]) le 0 do begin
    w0=[w0,i]
    i=i+1
  endwhile
  x=x[w0[-1]+1:-1,*]
  t=t[w0[-1]+1:-1]
endif

; #Interpolating data to 1 second (also interpolating zeroes elements)
tmin=long(min(t)) & tmax=long(max(t)) & if tmax ge 86400 then tmax=86399

nt=tmax-tmin+1		; 86400L if full day
if nt gt 86400 then nt=86400			
tt=dindgen(nt)+tmin
ydb=fltarr(nt,nf)
;zdb_plus=200.
zdb_plus=0.
for i=0,nf-1 do begin
  xx=reform(x(*,i))
  mx=DYN_N(xx,0.9995)
  wn0=where(xx gt 0. and xx le mx)

  yydb=interpol(10.*alog10(xx(wn0))+zdb_plus,t(wn0),tt)	; #+zdb_plus eventually, for having only positive dB values
  wnan=where(finite(yydb) eq 0)
  ydb(*,i)=yydb
endfor

;prime_decomp,nt,d,e,nprimes=10000
nprime=[85999,86011,86017,86027,86029,86069,86077,86083,86111,86113,86117,86131,86137,86143,86161,86171,86179,86183,86197,86201,86209,86239,86243,86249,86257,86263,86269,86287,86291,86293,86297,86311,86323,86341,86351,86353,86357,86369,86371,86381,86389,86399]
wnprime=where(nprime eq nt,nwnprime)
;while (n_elements(d) eq 1) do begin
while (nwnprime eq 1) do begin
  ydb=ydb[1:*,*]
  tt=tt[1:*]
  nt=nt-1
  wnprime=where(nprime eq nt,nwnprime)
  ;prime_decomp,nt,d,e,nprimes=10000
endwhile


; # constructing raw tables (dB and linear) #
rdb=ydb
rlin=10^((rdb-zdb_plus)/10.)

zdb=rdb
zlin=rlin
  
; # constructing frequency table for FFT
; # if yyyydddb is a periojve, then cut data in 2 hours interval (to avoid fft issues)
peri=[2016240,2016346,2016347,2017033,2017086,2017139,2017191,2017192,2017244,2017297,2017350,2018038,2018091,2018144,2018196,2018197,2018249,2018250,2018302,2018355,2019043,2019044,2019095,2019096,2019148,2019149,2019150]
if (where(peri eq yyyydddb))[0] ne -1 then begin
  
  ;iintmax=60
  iintmax=48
  print,"Perijove - FFT done on "+strtrim(long(nt/iintmax/60.),2)+" minutes window"
  for iint=0,iintmax-1 do begin
    if iint ne iintmax-1 then begin
      rdb2=rdb[nt/iintmax*iint:nt/iintmax*(iint+1)-1,*]
      rlin2=rlin[nt/iintmax*iint:nt/iintmax*(iint+1)-1,*]
    endif else begin
      rdb2=rdb[nt/iintmax*iint:-1,*]
      rlin2=rlin[nt/iintmax*iint:-1,*]
    endelse
  
    nt2=n_elements(rdb2[*,0])
    h=hanning(nt2,alpha=han,/double)
    freq=findgen(nt2)/nt2
    if nt2 mod 2 eq 0 then freq=[freq(0:nt2/2),-reverse(freq(1:nt2/2-1))] $
      else freq=[freq(0:nt2/2),-reverse(freq(1:nt2/2))]
  
    ; # Selection of harmonics
    yy=reform(rebin(rdb2,nt2,1))
    fy=abs(fft(yy*h,-1))
    w=where(freq ge 0.06 and freq le 0.07)
    ww=where(fy(w) eq max(fy(w)))
    w=w(ww(0))
    test=abs(freq/(freq(w)/2.))
    df=dharm*(freq(w)/2.)
    wwfonda=where(test le df*2./freq(w))
    wwlin=where(abs(freq) le df or ((test mod 1) gt df*2./freq(w) and (test mod 1) lt (1.-df*2./freq(w))) or test gt nharmlin+0.5)
    wwdb=where(abs(freq) le df or ((test mod 1) gt df*2./freq(w) and (test mod 1) lt (1.-df*2./freq(w))) or test gt nharmdb+0.5)
 
    ; # Filtering 
    
    for i=0,nf-1 do begin
      yydb=reform(rdb2(*,i))
      fydb=fft(yydb*h,-1)
      yylin=reform(rlin2(*,i))
      fylin=fft(yylin*h,-1)
      fydb(wwdb)=complex(0.,0.)
      fylin(wwlin)=complex(0.,0.)
      

      if iint ne iintmax-1 then begin
        zdb[nt/iintmax*iint:nt/iintmax*(iint+1)-1,i]=abs(fft(fydb*h,1));+min(rdb)
        zlin[nt/iintmax*iint:nt/iintmax*(iint+1)-1,i]=abs(fft(fylin*h,1))
      endif else begin
        zdb[nt/iintmax*iint:-1,i]=abs(fft(fydb*h,1));+min(rdb)
        zlin[nt/iintmax*iint:-1,i]=abs(fft(fylin*h,1))
      endelse
    endfor
  endfor


; # if yyyydddb is not a perijove, then FFT on the whole day
endif else begin
  h=hanning(nt,alpha=han,/double)
  freq=findgen(nt)/nt
  if nt mod 2 eq 0 then freq=[freq(0:nt/2),-reverse(freq(1:nt/2-1))] $
    else freq=[freq(0:nt/2),-reverse(freq(1:nt/2))]
  
  ; # Selection of harmonics
  yy=reform(rebin(rdb,nt,1))
  fy=abs(fft(yy*h,-1))
  w=where(freq ge 0.06 and freq le 0.07)
  ww=where(fy(w) eq max(fy(w)))
  w=w(ww(0))
  test=abs(freq/(freq(w)/2.))
  df=dharm*(freq(w)/2.)
  wwfonda=where(test le df*2./freq(w))
  wwlin=where(abs(freq) le df or ((test mod 1) gt df*2./freq(w) and (test mod 1) lt (1.-df*2./freq(w))) or test gt nharmlin+0.5)
  wwdb=where(abs(freq) le df or ((test mod 1) gt df*2./freq(w) and (test mod 1) lt (1.-df*2./freq(w))) or test gt nharmdb+0.5)

  if keyword_set(plot_filtering) then begin
    plotps,'select_harmonic_all_'+strtrim(long(yyyydddb),2)
      fyfilt=fy
      fyfilt[wwdb]=complex(0.,0.)
      plot_io,freq,fy,psym=3,/nodata,xtit='Frequency (Hz)',ytit=TexToIDL('Intensity (dB)'),tit='Fourrier Space - all channels',xrange=[min(freq),max(freq)],/xsty
      loadct,1
      oplot,freq,fy,psym=3,col=230                                                                                               
      loadct,0
      oplot,freq,fyfilt,psym=3
      
      yy=reform(rebin(rlin,nt,1))
      fy=abs(fft(yy*h,-1))
      fyfilt=fy
      fyfilt[wwlin]=complex(0.,0.)
      plot_io,freq,fy,psym=3,/nodata,xtit='Frequency (Hz)',ytit=TexToIDL('Intensity (V^2/m^2/Hz)'),tit='Fourrier Space - all channels',xrange=[min(freq),max(freq)],/xsty
      loadct,1
      oplot,freq,fy,psym=3,col=230                                                                                               
      loadct,0
      oplot,freq,fyfilt,psym=3
    endps,filename='select_harmonic_all_'+strtrim(long(yyyydddb),2)  
  endif
  ;stop
  ; # Filtering 

  if keyword_set(plot_filtering) then plotps,'select_harmonic_'+strtrim(long(yyyydddb),2)
  for i=0,nf-1 do begin
  
    yydb=reform(rdb(*,i))
    fydb=fft(yydb*h,-1)
    yylin=reform(rlin(*,i))
    fylin=fft(yylin*h,-1)
    if keyword_set(plot_filtering) then begin
      fydb_old=fydb
      fylin_old=fylin
    endif
    fydb(wwdb)=complex(0.,0.)
    fylin(wwlin)=complex(0.,0.)
    if keyword_set(plot_filtering) then begin
      if i eq 30 or i eq 50 or i eq 70 or i eq 91 then begin
      ;stop
        plot_io,freq,abs(fydb_old),psym=3,/nodata,xtit='Frequency (Hz)',ytit=TexToIDL('Intensity (dB)'),tit='Fourrier Space - Channel #'+strtrim(i+1,2)+' - '+strtrim(f[i],2)+' kHz',xrange=[min(freq),max(freq)],/xsty
        loadct,1
        oplot,freq,abs(fydb_old),psym=3,col=230
        loadct,0
        oplot,freq,abs(fydb),psym=3
        ;oplot,[1.*freq(w)/2.,1.*freq(w)/2.],[1e-12,1e2]
        ;oplot,[2.*freq(w)/2.,2.*freq(w)/2.],[1e-12,1e2]
        ;oplot,[3.*freq(w)/2.,3.*freq(w)/2.],[1e-12,1e2]
        ;oplot,[4.*freq(w)/2.,4.*freq(w)/2.],[1e-12,1e2]
        ;oplot,[5.*freq(w)/2.,5.*freq(w)/2.],[1e-12,1e2]
        ;oplot,[6.*freq(w)/2.,6.*freq(w)/2.],[1e-12,1e2]
        ;oplot,[7.*freq(w)/2.,7.*freq(w)/2.],[1e-12,1e2]
        ;oplot,[8.*freq(w)/2.,8.*freq(w)/2.],[1e-12,1e2]
        ;oplot,[-1.*freq(w)/2.,-1.*freq(w)/2.],[1e-12,1e2]
        ;oplot,[-2.*freq(w)/2.,-2.*freq(w)/2.],[1e-12,1e2]
        ;oplot,[-3.*freq(w)/2.,-3.*freq(w)/2.],[1e-12,1e2]
        ;oplot,[-4.*freq(w)/2.,-4.*freq(w)/2.],[1e-12,1e2]
        ;oplot,[-5.*freq(w)/2.,-5.*freq(w)/2.],[1e-12,1e2]
        ;oplot,[-6.*freq(w)/2.,-6.*freq(w)/2.],[1e-12,1e2]
        ;oplot,[-7.*freq(w)/2.,-7.*freq(w)/2.],[1e-12,1e2]
        ;oplot,[-8.*freq(w)/2.,-8.*freq(w)/2.],[1e-12,1e2]
  
  
        plot_io,freq,abs(fylin_old),psym=3,/nodata,xtit='Frequency (Hz)',ytit=TexToIDL('Intensity (V^2/m^2/Hz)'),tit='Fourrier Space - Channel #'+strtrim(i+1,2),xrange=[min(freq),max(freq)],/xsty
        loadct,1
        oplot,freq,abs(fylin_old),psym=3,col=230
        loadct,0
        oplot,freq,abs(fylin),psym=3
      ;  oplot,[1.*freq(w)/2.,1.*freq(w)/2.],[1e-12,1e2]
      ;  oplot,[2.*freq(w)/2.,2.*freq(w)/2.],[1e-12,1e2]
      ;  oplot,[3.*freq(w)/2.,3.*freq(w)/2.],[1e-12,1e2]
      ;  oplot,[4.*freq(w)/2.,4.*freq(w)/2.],[1e-12,1e2]
      ;  oplot,[5.*freq(w)/2.,5.*freq(w)/2.],[1e-12,1e2]
      ;  oplot,[6.*freq(w)/2.,6.*freq(w)/2.],[1e-12,1e2]
      ;  oplot,[7.*freq(w)/2.,7.*freq(w)/2.],[1e-12,1e2]
      ;  oplot,[8.*freq(w)/2.,8.*freq(w)/2.],[1e-12,1e2]
      ;  oplot,[-1.*freq(w)/2.,-1.*freq(w)/2.],[1e-12,1e2]
      ;  oplot,[-2.*freq(w)/2.,-2.*freq(w)/2.],[1e-12,1e2]
      ;  oplot,[-3.*freq(w)/2.,-3.*freq(w)/2.],[1e-12,1e2]
      ;  oplot,[-4.*freq(w)/2.,-4.*freq(w)/2.],[1e-12,1e2]
      ;  oplot,[-5.*freq(w)/2.,-5.*freq(w)/2.],[1e-12,1e2]
      ;  oplot,[-6.*freq(w)/2.,-6.*freq(w)/2.],[1e-12,1e2]
      ;  oplot,[-7.*freq(w)/2.,-7.*freq(w)/2.],[1e-12,1e2]
      ;  oplot,[-8.*freq(w)/2.,-8.*freq(w)/2.],[1e-12,1e2]
      endif
    endif
    zdb[*,i]=abs(fft(fydb*h,1));+min(rdb)
    zlin[*,i]=abs(fft(fylin*h,1))
  endfor
endelse


if keyword_set(plot_filtering) then endps,filename='select_harmonic_'+strtrim(long(yyyydddb),2)

; Background subtraction per frequency
; LIN
;MAKE_BACKGROUND, 10.*alog10(zlin),'',b,s
;s=10.^((b+s)/10.)-10.^(b/10.)
;b=10.^(b/10.)
; SUBTRACT_BACKGROUND, b,s, -1., zlin
;zlin2=zlin/rebin(reform(b,1,nf),nt,nf)
; dB
;MAKE_BACKGROUND, zdb,'',b,s
;zdb2=zdb
;SUBTRACT_BACKGROUND, b,s, -1.5, zdb2

; # filling intensity and time tables with fill values if nt  ne  86400
fillval=-1e-31
fillvaldB=-32767.0
if nt lt 86400L then begin
  if tmin gt 0 then begin
   zlin=[fltarr(tmin+1,nf)+fillval,zlin]
   rlin=[fltarr(tmin+1,nf)+fillval,rlin]
    zdb =[fltarr(tmin+1,nf)+fillvaldB,zdb]
    rdb =[fltarr(tmin+1,nf)+fillvaldB,rdb]
  endif
  if tmax lt 86399L then begin
    zlin=[zlin,fltarr(86400L-tmax,nf)+fillval]
    rlin=[rlin,fltarr(86400L-tmax,nf)+fillval]
    zdb =[zdb,fltarr(86400L-tmax,nf)+fillvaldB]
    rdb =[rdb,fltarr(86400L-tmax,nf)+fillvaldB]
  endif
endif

nt=86400L & tt=dindgen(nt)

; # rebinning to 15 seconds if asked
if keyword_Set(rebin15sec) then begin
  nt=5760
  zdb_rebin=fltarr(5760,nf) & zdb_rebin[*,*]=fillvaldB
  rdb_rebin=fltarr(5760,nf) & rdb_rebin[*,*]=fillvaldB
  zlin_rebin=fltarr(5760,nf) & zlin_rebin[*,*]=fillval
  rlin_rebin=fltarr(5760,nf) & rlin_rebin[*,*]=fillval
  for i=0l,5759l do begin
    for j=0,125 do begin
      wn0=where(zdb[15*i:15*(i+1)-1,j] gt fillvaldB)
      if wn0[0] ne -1 then begin
        zdb_rebin[i,j]=rebin(zdb[wn0+15*i,j],1,1)
        rdb_rebin[i,j]=rebin(rdb[wn0+15*i,j],1,1)
        zlin_rebin[i,j]=rebin(zlin[wn0+15*i,j],1,1)
        rlin_rebin[i,j]=rebin(rlin[wn0+15*i,j],1,1)
      endif
    endfor
  endfor
  zdb=zdb_rebin
  rdb=rdb_rebin
  zlin=zlin_rebin
  rlin=rlin_rebin

  
  t=rebin(tt/3600.,5760)
;#  zdb=rebin(zdb(*,*),5760,nf)
;#  zlin=rebin(zlin(*,*),5760,nf)
;#  rdb=rebin(rdb(*,*),5760,nf)
;#  rlin=rebin(rlin(*,*),5760,nf)
 ; #zdb2=rebin(zdb2(*,*),5760,nf)
 ; #zlin2=rebin(zlin2[*,*],5760,nf)
endif else begin
  t=tt
;#  f=f[wf]
;#  zdb=zdb[*,wf]
;#  zlin=zlin[*,wf]
;#  rdb=rdb[*,wf]
;#  rlin=rlin[*,wf]
;# zdb2=zdb2[*,wf]
;# zlin2=zlin2[*,wf]
endelse

;#to reuse if no FFT filtering on linear values (only dB)#
;stop,"decomment"
;rlin=fltarr(nt,nf)
;zlin=fltarr(nt,nf)
;rlin[*,*]=fillval
;zlin[*,*]=fillval
;wn0=where(rdb ne fillvaldB)
;rlin[wn0]=10^((rdb[wn0]-zdb_plus)/10.)
;zlin[wn0]=10^((10*zdb[wn0]-zdb_plus)/10.)
;#to reuse#

;# saving daily spdyn file
if keyword_set(filetest) then test_filename='test/' else test_filename=''
if keyword_set(version) then version_name='_v'+string(format='(I02)',version) else version_name=''
cmd = '../../data_n1/spdyn_sav_data'+deltat+'/'+strmid(strtrim(long(yyyydddb),2),0,4)
spawn,cmd,resu
filename='../../data_n1/spdyn_sav_data'+deltat+'/'+strmid(strtrim(long(yyyydddb),2),0,4)+'/'+test_filename+strtrim(long(yyyydddb),2)+'_spdyn'+deltat+version_name+'.sav'
print,"# Saving file "+filename+" #"
yyyyddd=yyyydddb
save,yyyyddd,rdb,rlin,zlin,zdb,t,f,filename=filename
print,"# File saved #"

; # plotting daily data if asked
if keyword_set(plot_spdyn) then begin
  cmd = 'mkdir pdf_spdyn_survey'
  spawn,cmd,resu

  set_ps,'./pdf_spdyn_survey/'+strtrim(long(yyyydddb),2)+'.ps',/landscape
  device,/col
  loadct,col
  !p.multi=[0,1,2]
  c=0 & if col ne 0 then c=1
  SPDYNPS,rebin(zdb2 ,5760,1100),0,24,16,126,'Time (hour)','Frequency (WAVES channel #)',strtrim(long(yyyydddb),2)+'   Intensity Log-filtered',0,0,0,0.,0.97,c,'dB'
  SPDYNPS,rebin(zlin2,5760,1100),0,24,16,126,'Time (hour)','Frequency (WAVES channel #)',strtrim(long(yyyydddb),2)+'   Intensity Lin-filtered',0,0,0,0.01,0.95,c,' '
  device,/close
  set_plot,'x'
  !p.font=-1
  ;ps_pdf,'./pdf_spdyn_survey/'+strtrim(long(yyyydddb),2)+'.ps'
  cmd='ps2pdf ./pdf_spdyn_survey/'+strtrim(long(yyyydddb),2)+'.ps ./pdf_spdyn_survey/'+strtrim(long(yyyydddb),2)+'.pdf'
  spawn,cmd
  ;cmd='rm ./pdf_spdyn_survey/'+strtrim(long(yyyydddb),2)+'.ps'
  ;spawn,cmd
endif
end
