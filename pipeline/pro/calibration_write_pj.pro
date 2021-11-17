PRO calibration_write_PJ,gain_new,temporal_window=temporal_window,ps=ps,yyyydddb=yyyydddb,yyyyddde=yyyyddde,linear=linear,db=db,version=version

version_name='_v'+string(format='(I02)',version)
if ~keyword_set(yyyydddb) then yyyydddb=2016100
if ~keyword_set(yyyyddde) then yyyyddde=2018064

if keyword_Set(linear) then begin
  data_type_treat='lin'
  data_type_folder='linear'
endif else if keyword_Set(dB) then begin
  data_type_treat='db'
  data_type_folder='db'
endif


print,"### Restoring the FFT-filtered data of HFR-Low during Perijove (using the catalog selection 'PHC') ###"
restore,'../make_timeseries/PHC_linear_noback_subtract/PHC_timeseries_d60_channels_61-87_zlin'+version_name+'.sav'
print,'### HFR-Low FFT-filtered data timeseries restored ###'

s_HFR=timeseries
time_HFR=time
freq_HFR=frequencies

print,"### Restoring background file ###"
restore,'../background/background_ALLPJ_zlin'+version_name+'.sav',/verb
print,'### background restored ###'
timeseries2=transpose(s_HFR)

print,'### Subtracting background for HFR-Low sub-receiver ###'
SUBTRACT_BACKGROUND, background[45:71],sigma,0, timeseries2
s_HFR=transpose(timeseries2)
print,'### Background subtracted ###'

print,'### Restoring the calibration timeseries ('+data_type_folder+', back subtracted, no normalized) ###'
restore,'../calibration_timeseries/ALL_calibration_timeseries_'+data_type_folder+'_back_subtract_d60_z'+data_type_treat+'nb'+version_name+'.sav',/verb
s_tmp=s
nchannels=n_elements(s_tmp[*,0])
nt=long(n_elements(time))
print,'### Calibration timeseries restored ###'

print,'### Restoring calibration gain from Juno/Cassini+Voyager comparison ###'
restore,'../gain/gain_juno_casvg_'+data_type_treat+version_name+'.sav',/verb 
print,'### Calibration gain restored ###'


s_tmp[45:71,*]=s_HFR[*,*]

temporal_window=float(temporal_window)
peri=[2016240,2016346,2017033,2017086,2017139,2017192,2017244,2017297,2017350,2018038,2018091,2018144,2018197,2018250,2018302,2018355,2019043,2019096,2019149]
peri=peri[where(peri ge yyyydddb and peri le yyyyddde)]
hhperi=[double(12.+50./60.),double(17.+04./60.),double(12.+57./60.),double(08.+52./60.),double(06.+01./60.),double(01.+54./60.),double(21.+49./60.),$
		double(17.+42./60.),double(17.+57./60.),double(13.+52./60.),double(09.+46./60.),double(05.+40./60.),double(05.+17./60.),double(01.+12./60.),$
		double(21.+06./60.),double(17.+01./60.),double(17.+34./60.),double(12.+14./60.),double(08.+08./60.)]
peri=aj_t16(peri)+hhperi/24.


nn=n_elements(peri)
period=2.*temporal_window*60.
s2=dblarr(nchannels,nn,period)
time2=dblarr(nn,period)

xmed=fltarr(nn,nchannels)
x1=fltarr(nn,nchannels)

;# apply calbratin gain to LFR-Low&High & HFR-High subreceivers
s_tmp[0:44,*]=s[0:44,*]*rebin(gain[0:44],n_elements(gain[0:44]),nt)    
s_tmp[72:nchannels-1,*]=s[72:nchannels-1,*]*rebin(gain[72:nchannels-1],n_elements(gain[72:nchannels-1]),nt)



; # sort data by time interval equal to the temporal window
for i=0,n_elements(peri)-1 do begin
	wtime=where(time ge peri[i]-temporal_window/24. and time le peri[i]+temporal_window/24.)
	s2[*,i,*]=s_tmp[*,wtime]
	time2[i,*]=time[wtime]
endfor

; # calculing the median of the 50% and 1% levels occurrence
for iint=0,nn-1 do begin
    for ichannel=0,nchannels-1 do begin
    ;# computing spectra of 50% and 1% occurrence levels 
    ;# (i.e., at each frequency the flux density exceeded 50% and 1% of the time) 
        wn0=where(s2[ichannel,iint,*] gt 0.)
        if wn0[0] ne -1 then begin
          xmed[iint,ichannel]=median(s2[ichannel,iint,wn0])
          x1[iint,ichannel]=dyn_n(s2[ichannel,iint,wn0],0.99)
        endif
    endfor
endfor

xmedmed=fltarr(110) 
x1med=fltarr(110) 
for i=0,n_elements(xmedmed)-1 do begin
  wn0=where(xmed[*,i] gt 0.)
  if wn0[0] ne -1 then xmedmed[i]=median(xmed[wn0,i])
  wn0=where(x1[*,i] gt 0.)
  if wn0[0] ne -1 then x1med[i]=median(x1[wn0,i])
endfor


xmed_tmp=xmed
xmedmed_tmp=xmedmed
x1_tmp=x1
x1med_tmp=x1med
fx=frequencies


; # step to calibrate HFR-Low:  signal must be continue between LFR-High and HFR-Low and comparison to Voyager-Cassini 50% and 1% occurrence level

 restore,'m01.sav'
  f01=xx & m01=alog10(yy)
  c1m=interpol(m01,f01,fx)
  restore,'m50.sav'
  f50=xx & m50=alog10(yy)
  cmm=interpol(m50,f50,fx)
  w5=where(abs(fx-5) eq min(abs(fx-5))) & w5=w5(0)
  cmm(where(fx lt 5))=c1m(where(fx lt 5))-c1m(w5)+cmm(w5)	; Cassini<5 kHz: 50% ajusté comme 1%
  restore,'z92.sav'
  f92=xx & m92=alog10(yy)
  c92=interpol(m92,f92,fx)
  w=where(abs(c92-c1m) eq min(abs(c92-c1m)))	; correction 1% > 10 MHz, étendue à 41 MHz, idem 50%  >30 MHz 
  w=w(0)
  cmm=[cmm(0:w-1),c92(w:*)-c1m(w)+cmm(w)]
  c1m=[c1m(0:w-1),c92(w:*)]

cmm_UP=cmm[45:71]-cmm[72]+alog10(xmedmed[72])
cmm_DOWN=cmm[45:71]-cmm[44]+alog10(xmedmed[44])
cmm_geom=cmm(45:71)-(cmm(44)+cmm(72))/2+(alog10(xmedmed[44])+alog10(xmedmed[72]))/2

c1m_UP=c1m[45:71]-c1m[72]+alog10(x1med[72])
c1m_DOWN=c1m[45:71]-c1m[44]+alog10(x1med[44])
c1m_geom=c1m(45:71)-(c1m(44)+c1m(72))/2+(alog10(x1med[44])+alog10(x1med[72]))/2

gain_hfr_low_cmm=10^cmm_geom/xmedmed[45:71]
gain_hfr_low_c1m=10^c1m_geom/x1med[45:71]


;# Calibrated using an interpolation between the last frequency channel of LFR-High and the first one of HFR-High
;gain_hfr_low=(findgen(n_elements(freq[44:72]),start=gain[44],inc=(-(gain[44]-gain[72])/(n_elements(freq[44:72])-1))))[1:-2]
gain_hfr_low=(10^(findgen(n_elements(freq[45:72]),start=alog10(gain[44]),inc=-(alog10(gain[44])-alog10(gain[72]))/27)))[0:-2]
;# Only calibrated on the 1% occurrence level
;gain_hfr_low=gain_hfr_low_c1m



;# calibration gain only based on cmm
;gain_hfr_low=10^cmm_geom/xmedmed[45:71]


;# Calibrated on the 1% &nd 50% occurrence levels
;gain_hfr_low=sqrt(gain_hfr_low_cmm*gain_hfr_low_c1m)

;# Only calibrated on the 50% occurrence level
;gain_hfr_low=gain_hfr_low_cmm
;# Calibrated by overlapped the last frequency channel of HFR-Low to the first one of HFR-High, using the 1% level occurrence
;gain_hfr_low=fltarr(n_elements(freq[45:71]))
;gain_hfr_low[*]=10^c1m_geom[-1]/x1med[71]
;gain_hfr_low[*]=sqrt(x1med[44]/x1med[45]*x1med[72]/x1med[71])
;# Calibrated by using the 1% occurrence level signal observed between 700 kHz to 2800 KHz
;gain_hfr_low=fltarr(n_elements(freq[45:71]))
;gain_hfr_low[*]=(total(gain_hfr_low_c1m[15:-1]))/n_elements(gain_hfr_low_c1m[15:-1])
;# Calibrated by using the 1% and 50% occurrence levels signal observed between 700 kHz to 2800 KHz
;gain_hfr_low=fltarr(n_elements(freq[45:71]))
;gain_hfr_low[*]=sqrt(total(gain_hfr_low_c1m[14:-1])/n_elements(gain_hfr_low_c1m[14:-1])*total(gain_hfr_low_cmm[14:-1])/n_elements(gain_hfr_low_cmm[14:-1]))

xmed_tmp[*,45:71]=xmed_tmp[*,45:71]*transpose(rebin(gain_hfr_low,n_elements(gain_hfr_low),n_elements(xmed_tmp[*,0])))
x1_tmp[*,45:71]=x1_tmp[*,45:71]*transpose(rebin(gain_hfr_low,n_elements(gain_hfr_low),n_elements(x1_tmp[*,0])))

for i=0,n_elements(xmedmed_tmp)-1 do begin
  wn0=where(xmed_tmp[*,i] gt 0.)
  if wn0[0] ne -1 then xmedmed_tmp[i]=median(xmed_tmp[wn0,i])
  wn0=where(x1_tmp[*,i] gt 0.)
  if wn0[0] ne -1 then x1med_tmp[i]=median(x1_tmp[wn0,i])
endfor
;#for i=0,109 do xmedmed_tmp[i]=median(xmed_tmp[where(xmed_tmp[*,i] gt 0.),i])


;x1_tmp=x1
;x1_tmp[*,45:71]=x1[*,45:71]*gain_hfr_low
;x1med_tmp=fltarr(110)
;for i=0,109 do x1med_tmp[i]=median(x1_tmp[*,i])

time=time2

gain_final=gain
gain_final[45:71]=gain_hfr_low
gain_first_calibration=gain
print,'Saving...'
save,time,fx,xmed,x1,filename='../calibration_data/ALL_calibration_HFR-Low_data_freq_'+data_type_folder+'_back_subtract_d60_z'+data_type_treat+'nb.sav'
freq=frequencies
save,filename='../gain/gain_final_'+data_type_treat+'.sav',freq,gain_final,gain_first_calibration,gain_hfr_low
print,'Save complete'

if keyword_set(ps) then begin
	print,'Plotting'
	plotps,'pdf/intercal_pj_'+data_type_folder
	!p.multi=[0,2,2]
		loadct,0
		plot_oo,frequencies,(xmed[0,*]),/xsty,yrange=[min((xmed)),(max(xmed))],/ysty,/nodata,xtit='Frequency(kHz)',ytit=TexToIDL('Intensity (V^2/m^2/Hz)'),tit='Juno/Waves 50% level during perijoves !c HFR-Low&High & HFR-High calibrated on Cassini & Voyager',charsize=1.75
		oplot,[frequencies[45],frequencies[45]],[min((xmed)),max((xmed))]
		oplot,[frequencies[71],frequencies[71]],[min((xmed)),max((xmed))]
		loadct,1
		for i=0,n_elements(xmed[*,0])-1 do $
			oplot,frequencies,(xmed[i,*]),col=220
		
		loadct,53
    	oplot,[fx[27],fx[27]],[1e-24,1e-6],col=120,linestyle=2
    	oplot,[fx[45],fx[45]],[1e-24,1e-6],col=120,linestyle=2
    	oplot,[fx[71],fx[71]],[1e-24,1e-6],col=120,linestyle=2
    	oplot,[fx[72],fx[72]],[1e-24,1e-6],col=120,linestyle=2
    	al_legend,'LFR-Low',position=[2.5,1e-16*0.8],box=0,textcolors=120
    	al_legend,'LFR-High',position=[25,1e-16*0.8],box=0,textcolors=120
    	al_legend,'HFR-Low',position=[350,1e-16*0.8],box=0,textcolors=120
    	al_legend,'HFR-High',position=[6000,1e-16*0.8],box=0,textcolors=120
		loadct,0
		al_legend,'HFR-Low !cnot calibrated',position=[250,1e-8*0.8],box=0
		oplot,frequencies,(xmedmed)





		plot_oo,frequencies,xmedmed,title='Step to calibrate HFR-Low: !cuse the overlap between LFR-High & LFR-Low !cand the Cassini-Voyager 50% level spectra',ytit=TexToIDL('Intensity (V^2/m^2/Hz)'), xtit='Frequency (kHz)',/xstyle,yrange=[min((xmed)),(max(xmed))],/ysty,/nodata,charsize=1.75
		loadct,1
		for i=0,n_elements(xmed[*,0])-1 do oplot,frequencies,xmed[i,*],col=220
		loadct,0
		oplot,frequencies,xmedmed
		loadct,62
		oplot,frequencies[45:71],10^(cmm[45:71]-cmm[44]+alog10(xmedmed[44])),linestyle=2,color=225
		oplot,frequencies[45:71],10^(cmm[45:71]-cmm[72]+alog10(xmedmed[72])),linestyle=2,color=225
		oplot,frequencies[45:71],10^cmm_geom,color=225
		loadct,53
    	oplot,[fx[27],fx[27]],[1e-24,1e-6],col=120,linestyle=2
    	oplot,[fx[45],fx[45]],[1e-24,1e-6],col=120,linestyle=2
    	oplot,[fx[71],fx[71]],[1e-24,1e-6],col=120,linestyle=2
    	oplot,[fx[72],fx[72]],[1e-24,1e-6],col=120,linestyle=2
    	al_legend,'LFR-Low',position=[2.5,1e-16*0.8],box=0,textcolors=120
    	al_legend,'LFR-High',position=[25,1e-16*0.8],box=0,textcolors=120
    	al_legend,'HFR-Low',position=[350,1e-16*0.8],box=0,textcolors=120
    	al_legend,'HFR-High',position=[6000,1e-16*0.8],box=0,textcolors=120
		loadct,0
		al_legend,'HFR-Low !cnot calibrated',position=[250,1e-8*0.8],box=0
		al_legend,'Reference spectrum',position=[4e2,3e-12],box=0,textcolors=120



		plot_oo,frequencies,xmedmed_tmp,yrange=[min(xmed_tmp),max(xmed_tmp)],/xsty,/ysty,/nodata,xtit='Frequency(kHz)',ytit=TexToIDL('Intensity (V^2/m^2/Hz)'),tit='Juno/Waves 50% level during PJ!cData after gain application',charsize=1.75
		loadct,1
		for i=0,n_elements(xmed_tmp[*,0])-1 do oplot,frequencies,(xmed_tmp[i,*]),col=220
		loadct,53
    	oplot,[fx[27],fx[27]],[1e-24,1e-6],col=120,linestyle=2
    	oplot,[fx[45],fx[45]],[1e-24,1e-6],col=120,linestyle=2
    	oplot,[fx[71],fx[71]],[1e-24,1e-6],col=120,linestyle=2
    	oplot,[fx[72],fx[72]],[1e-24,1e-6],col=120,linestyle=2
    	al_legend,'LFR-Low',position=[2.5,1e-16*0.8],box=0,textcolors=120
    	al_legend,'LFR-High',position=[25,1e-16*0.8],box=0,textcolors=120
    	al_legend,'HFR-Low',position=[350,1e-16*0.8],box=0,textcolors=120
    	al_legend,'HFR-High',position=[6000,1e-16*0.8],box=0,textcolors=120
		loadct,0
		
		oplot,frequencies,xmedmed_tmp



    ;# 1% occurrence level
    !p.multi=[0,2,2]
    loadct,0
    plot_oo,frequencies,(x1[0,*]),/xsty,yrange=[min((x1)),(max(x1))],/ysty,/nodata,xtit='Frequency(kHz)',ytit=TexToIDL('Intensity (V^2/m^2/Hz)'),tit='Juno/Waves 1% level during perijoves !c HFR-Low&High & HFR-High calibrated on Cassini & Voyager',charsize=1.75
    oplot,[frequencies[45],frequencies[45]],[min((x1)),max((x1))]
    oplot,[frequencies[71],frequencies[71]],[min((x1)),max((x1))]
    loadct,1
    for i=0,n_elements(x1[*,0])-1 do $
      oplot,frequencies,(x1[i,*]),col=220
    
    loadct,53
      oplot,[fx[27],fx[27]],[1e-24,1e-6],col=120,linestyle=2
      oplot,[fx[45],fx[45]],[1e-24,1e-6],col=120,linestyle=2
      oplot,[fx[71],fx[71]],[1e-24,1e-6],col=120,linestyle=2
      oplot,[fx[72],fx[72]],[1e-24,1e-6],col=120,linestyle=2
      al_legend,'LFR-Low',position=[2.5,1e-16*0.8],box=0,textcolors=120
      al_legend,'LFR-High',position=[25,1e-16*0.8],box=0,textcolors=120
      al_legend,'HFR-Low',position=[350,1e-16*0.8],box=0,textcolors=120
      al_legend,'HFR-High',position=[6000,1e-16*0.8],box=0,textcolors=120
    loadct,0
    al_legend,'HFR-Low !cnot calibrated',position=[250,1e-8*0.8],box=0
    oplot,frequencies,(x1med)





    plot_oo,frequencies,x1med,title='Step to calibrate HFR-Low: !cuse the overlap between LFR-High & LFR-Low !cand the Cassini-Voyager 1% level spectra',ytit=TexToIDL('Intensity (V^2/m^2/Hz)'), xtit='Frequency (kHz)',/xstyle,yrange=[min((x1)),(max(x1))],/ysty,/nodata,charsize=1.75
    loadct,1
    for i=0,n_elements(x1[*,0])-1 do oplot,frequencies,x1[i,*],col=220
    loadct,0
    oplot,frequencies,x1med
    loadct,62
    oplot,frequencies[45:71],10^(c1m[45:71]-c1m[44]+alog10(x1med[44])),linestyle=2,color=225
    oplot,frequencies[45:71],10^(c1m[45:71]-c1m[72]+alog10(x1med[72])),linestyle=2,color=225
    oplot,frequencies[45:71],10^c1m_geom,color=225
    loadct,53
      oplot,[fx[27],fx[27]],[1e-24,1e-6],col=120,linestyle=2
      oplot,[fx[45],fx[45]],[1e-24,1e-6],col=120,linestyle=2
      oplot,[fx[71],fx[71]],[1e-24,1e-6],col=120,linestyle=2
      oplot,[fx[72],fx[72]],[1e-24,1e-6],col=120,linestyle=2
      al_legend,'LFR-Low',position=[2.5,1e-16*0.8],box=0,textcolors=120
      al_legend,'LFR-High',position=[25,1e-16*0.8],box=0,textcolors=120
      al_legend,'HFR-Low',position=[350,1e-16*0.8],box=0,textcolors=120
      al_legend,'HFR-High',position=[6000,1e-16*0.8],box=0,textcolors=120
    loadct,0
    al_legend,'HFR-Low !cnot calibrated',position=[250,1e-8*0.8],box=0
    al_legend,'Reference spectrum',position=[4e2,3e-12],box=0,textcolors=120



    plot_oo,frequencies,x1med_tmp,yrange=[min(x1_tmp),max(x1_tmp)],/xsty,/ysty,/nodata,xtit='Frequency(kHz)',ytit=TexToIDL('Intensity (V^2/m^2/Hz)'),tit='Juno/Waves 1% level during PJ!cData after gain application',charsize=1.75
    loadct,1
    for i=0,n_elements(x1_tmp[*,0])-1 do oplot,frequencies,(x1_tmp[i,*]),col=220
    loadct,53
      oplot,[fx[27],fx[27]],[1e-24,1e-6],col=120,linestyle=2
      oplot,[fx[45],fx[45]],[1e-24,1e-6],col=120,linestyle=2
      oplot,[fx[71],fx[71]],[1e-24,1e-6],col=120,linestyle=2
      oplot,[fx[72],fx[72]],[1e-24,1e-6],col=120,linestyle=2
      al_legend,'LFR-Low',position=[2.5,1e-16*0.8],box=0,textcolors=120
      al_legend,'LFR-High',position=[25,1e-16*0.8],box=0,textcolors=120
      al_legend,'HFR-Low',position=[350,1e-16*0.8],box=0,textcolors=120
      al_legend,'HFR-High',position=[6000,1e-16*0.8],box=0,textcolors=120
    loadct,0
    
  
    oplot,frequencies,x1med_tmp

    !p.multi=0




    plot_oo,frequencies,xmedmed,title='Step to calibrate HFR-Low: !cuse the overlap between LFR-High & LFR-Low !cand the Cassini-Voyager 50% level spectra',ytit=TexToIDL('Intensity (V^2/m^2/Hz)'), xtit='Frequency (kHz)',/xstyle,yrange=[1e-19, 1e-5],/nodata,charsize=1.75
    loadct,1
    for i=0,n_elements(xmed[*,0])-1 do oplot,frequencies,xmed[i,*],col=220
    loadct,0
    oplot,frequencies,xmedmed
    loadct,62
    oplot,frequencies[45:71],10^(cmm[45:71]-cmm[44]+alog10(xmedmed[44])),linestyle=2,color=225
    oplot,frequencies[45:71],10^(cmm[45:71]-cmm[72]+alog10(xmedmed[72])),linestyle=2,color=225
    oplot,frequencies[45:71],10^cmm_geom,color=225
    loadct,53
      oplot,[fx[27],fx[27]],[1e-24,1e-5],col=120,linestyle=2
      oplot,[fx[45],fx[45]],[1e-24,1e-5],col=120,linestyle=2
      oplot,[fx[71],fx[71]],[1e-24,1e-5],col=120,linestyle=2
      oplot,[fx[72],fx[72]],[1e-24,1e-5],col=120,linestyle=2
      al_legend,'LFR-Low',position=[2.5,1e-19*0.8],box=0,textcolors=120
      al_legend,'LFR-High',position=[25,1e-19*0.8],box=0,textcolors=120
      al_legend,'HFR-Low',position=[350,1e-19*0.8],box=0,textcolors=120
      al_legend,'HFR-High',position=[6000,1e-19*0.8],box=0,textcolors=120
    loadct,0
    al_legend,'HFR-Low !cnot calibrated',position=[250,1e-6*0.8],box=0
    



    plot_oo,frequencies,x1med,title='Juno/Waves 1% level',ytit=TexToIDL('Flux density (W/m^2/Hz)'), xtit='Frequency (kHz)',/xstyle,yrange=[1e-19, 1e-3],/nodata,charsize=1.75,/ystyle
    loadct,1
    for i=0,n_elements(x1[*,0])-1 do oplot,frequencies,x1[i,*],col=220
    loadct,0
    oplot,frequencies,x1med
;    loadct,62
    oplot,frequencies[45:71],10^(c1m[45:71]-c1m[44]+alog10(x1med[44])),linestyle=2;,color=225
    oplot,frequencies[45:71],10^(c1m[45:71]-c1m[72]+alog10(x1med[72])),linestyle=2;,color=225
    oplot,frequencies[45:71],10^c1m_geom,linestyle=2;,color=225
    loadct,53
      oplot,[fx[27],fx[27]],[1e-24,1e-3],col=120,linestyle=2
      oplot,[fx[45],fx[45]],[1e-24,1e-3],col=120,linestyle=2
      oplot,[fx[71],fx[71]],[1e-24,1e-3],col=120,linestyle=2
      oplot,[fx[72],fx[72]],[1e-24,1e-3],col=120,linestyle=2
      al_legend,'LFR-Low',position=[2.5,1e-17*0.8],box=0,textcolors=120
      al_legend,'LFR-High',position=[25,1e-17*0.8],box=0,textcolors=120
      al_legend,'HFR-Low',position=[350,1e-17*0.8],box=0,textcolors=120
      al_legend,'HFR-High',position=[6000,1e-17*0.8],box=0,textcolors=120
    loadct,0
    al_legend,'HFR-Low !cnot calibrated',position=[250,1e-6*0.8],box=0
    

    plot_oo,frequencies[0:44], g01[0:44],yrange=[0.1,10000],xrange=[min(frequencies),max(frequencies)],xtit='Frequency (kHz)',ytit='Final conversion factor',tit='Juno Waves final conversion factor',/xsty,/ysty,charsize=1.75
    oplot,frequencies[72:-1],g01[72:-1]
    oplot,frequencies[45:71],gain_hfr_low_c1m
    loadct,13
    oplot,frequencies[0:44],g50[0:44],col=255
    oplot,frequencies[72:-1],g50[72:-1],col=255
    oplot,frequencies[45:71],gain_hfr_low_cmm,col=255
    loadct,51
    oplot,frequencies[0:44],gain_final[0:44],thick=5,col=240
    oplot,frequencies[0:14],gain_with_trapped_continuum[0:14],thick=5,col=240,linestyle=2
    oplot,frequencies[45:71],gain_final[45:71],thick=5,col=240
    oplot,frequencies[72:-1],gain_final[72:-1],thick=5,col=240
    loadct,53
    oplot,[fx[27],fx[27]],[0.1,10000],col=120,linestyle=2
    oplot,[fx[45],fx[45]],[0.1,10000],col=120,linestyle=2
    oplot,[fx[71],fx[71]],[0.1,10000],col=120,linestyle=2
    oplot,[fx[72],fx[72]],[0.1,10000],col=120,linestyle=2
    al_legend,'LFR-Low',position=[2.5,0.2],box=0,textcolors=120
    al_legend,'LFR-High',position=[25,0.2],box=0,textcolors=120
    al_legend,'HFR-Low',position=[350,0.2],box=0,textcolors=120
    al_legend,'HFR-High',position=[6000,0.2],box=0,textcolors=120
    loadct,0
    oplot,[frequencies[0],frequencies[-1]],[1.0,1.0],linestyle=1
      



	endps,filename='pdf/intercal_pj_'+data_type_folder
	print, 'Plots complete'
endif
return
END