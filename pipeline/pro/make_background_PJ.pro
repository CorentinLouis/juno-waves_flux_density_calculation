PRO make_background_PJ,PJ=PJ,cut_timeseries_PJ=cut_timeseries_PJ,allPJ=allPJ,LINEAR=LINEAR,dB=dB,version=version

if keyword_set(version) then version_name='_v'+string(format='(I02)',version) else version_name='_v01'

if keyword_Set(linear) then begin
  data_type_treat='lin'
  data_type_folder='linear'
endif else if keyword_Set(dB) then begin
  data_type_treat='db'
  data_type_folder='db'
endif

if ~keyword_set(PJ) and ~keyword_set(allPJ) then begin
	allPJ=1
	print,"The background will be compute using all data"
endif
if keyword_set(PJ) then begin
	if keyword_set(cut_timeseries_PJ) then cut_timeseries_PJ,linear=linear,dB=dB
	peri=[amj_aj(20160704),2016240,2016346,2017033,2017086,2017139,2017192,2017244,2017297,2017350,2018038,2018091,2018144,2018197,2018250,2018302,2018355,2019043,2019096,2019149]
	nperi=n_elements(peri)
	background2=fltarr(nperi,110)
	sigma2=fltarr(nperi,110)
	for i=0,nperi-1 do begin
		restore,'../make_timeseries/ALL_'+data_type_folder+'_noback_subtract/PJ'+strtrim(i,2)+'_timeseries_d60_channels_16-125_z'+data_type_treat+version_name+'.sav',/verb
	
		MAKE_BACKGROUND, 10.*alog10(transpose(timeseries2)),'',b,s,/supmin
		sigma=10.^((b+s)/10.)-10.^(b/10.)
		background=10.^(b/10.)
		background2[i,*]=background
		sigma2[i,*]=sigma
		time=time2
		save, filename='../background/background_PJ'+strtrim(i,2)+'_'+data_type_treat+'.sav',time,frequencies,background,sigma
		                                                                                        
	endfor
	save, filename='../background/background_cut_PJ0-'+strtrim(i,2)+'_'+data_type_treat+'.sav',frequencies,background2,sigma2
	
	plotps,'pdf/background_cut_pj_'+data_type_treat  
	plot_oo,frequencies,background2[0,*],/xstyle,/nodata,tit='Background for each PJ + inbound trajectory',xtit='Frequency (kHz)',ytit=TexToIDL('Intensity (V^2/m^2/Hz)'),yrange=[1e-18, 1e-12],/ystyle
	for i=0,n_elements(background2[*,0])-1 do oplot,frequencies,background2[i,*],col=180
	sigmabackg=fltarr(n_elements(frequencies))  
	meanback=fltarr(n_elements(frequencies)) 
	for i=0,n_elements(frequencies)-1 do meanback[i]=10^mean(alog10(background2[*,i])) 
	for i=0,n_elements(frequencies)-1 do sigmabackg[i]=stddev(alog10(background2[*,i]))
	oplot,frequencies,meanback
	oplot,frequencies,10^(alog10(meanback)+sigmabackg),linestyle=2
	oplot,frequencies,10^(alog10(meanback)-sigmabackg),linestyle=2
	endps,filename='pdf/background_cut_pj_'+data_type_treat  

endif else if keyword_set(allPJ) then begin
	restore,'../make_timeseries/ALL_'+data_type_folder+'_noback_subtract/ALL_timeseries_d60_channels_16-125_z'+data_type_treat+version_name+'.sav',/verb
	wsup0=where(timeseries[0,*] gt 0.)
	MAKE_BACKGROUND, 10.*alog10(transpose(timeseries[*,wsup0])),'',b,s,/supmin
	
	s=10.^((b+s)/10.-10.^(b/10.))
	b=10.^(b/10.)
	background=b
	sigma=s
	save, filename='../background/background_ALLPJ_z'+data_type_treat+version_name+'.sav',time,frequencies,background,sigma
endif
END