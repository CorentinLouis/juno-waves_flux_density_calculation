;----------------------------------------------------------------
  pro INTERCAL, freq,gain, PS=PS,linear=linear,db=db,smoothing=smoothing,version=version
;----------------------------------------------------------------
; freq = Juno Waves frequencies 1 kHz - 40.5 MHz
; gain = Juno to flux densities conversion factor
  

  version_name='_v'+string(format='(I02)',version)

  if keyword_Set(linear) then begin
    data_type_treat='lin'
    data_type_folder='linear'
  endif else if keyword_Set(dB) then begin
    data_type_treat='db'
    data_type_folder='db'
  endif

  if keyword_set(PS) then begin
    plotps,'pdf/intercal_'+data_type_treat
    loadct,13
    ;!p.multi=[0,1,2]
  endif
  restore,'../calibration_data/ALL_calibration_data_freq_'+data_type_folder+'_back_subtract_d60_z'+data_type_treat+'nb_Rmin30_MLatmax15'+version_name+'.sav',/verb
  
  nfx=indgen(110)+16
  nrec=intarr(110)
  nrec(where(nfx le 42))=1
  nrec(where(nfx ge 43 and nfx le 60))=2
  nrec(where(nfx ge 61 and nfx le 87))=3
  nrec(where(nfx ge 88))=4
  
  if not(keyword_set(PS)) then window,0,xs=1000,ys=400
  
  xmm=fltarr(n_elements(fx))
  
  for i=0,n_elements(fx)-1 do begin
    wn0=where(xmed[*,i] gt 0.)
    if wn0[0] ne -1 then xmm[i]=alog10(median(xmed[wn0,i]))
  endfor

  ;xmm=alog10(median(xmed,dimension=1))
  plot_oo,fx,10^xmm,psym=-4,/xsty,yra=[1e-24,1e-16],/ysty,xtit='Frequency (kHz)',ytit=TexToIdl('Flux density @ 1 AU (W/m^2/Hz)'),tit='Juno/Waves 50% and 1% levels',/nodata,charsize=1.75
  oplot,fx,10^xmm,psym=-4
  ;oplot,fx[0:44],10^xmm[0:44],psym=-4
  ;oplot,fx[72:-1],10^xmm[72:-1],psym=-4
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
  loadct,13
  
  xmms=xmm
  for i=1,3 do begin
    w=where(nrec eq i)
    ;oplot,fx(w),10^xmm(w),psym=-4,color=250
    res=poly_fit(alog10(fx(w)),xmm(w),5,yfit=yfit)
    xmms(w)=yfit
    ;oplot,fx(w),10^xmms(w)
  endfor
  ;oplot,fx,10^xmms,col=250
  w=where(nrec eq 4)
  ww=[indgen(17),indgen(17)+19,37]
  ;oplot,fx(w(ww)),xmm(w(ww)),psym=-4,color=250
  
  y=interpol(xmm(w(ww)),fx(w(ww)),fx(w))
  if keyword_set(smoothing) then xmms(w)=smooth(y,3,/edge_trunc)
  ;oplot,fx(w),10^xmms(w),col=250
  oplot,fx,10^xmms,col=250,thick=5
  ;oplot,fx[0:44],10^xmms[0:44],col=250,thick=5
  ;oplot,fx[72:-1],10^xmms[72:-1],col=250,thick=5
  x1m=fltarr(n_elements(fx))
  for i=0,n_elements(fx)-1 do begin
    wn0=where(x1[*,i] gt 0.)
    if wn0[0] ne -1 then x1m[i]=alog10(median(x1[wn0,i]))
  endfor
  ;x1m=alog10(median(x1,dimension=1))
  oplot,fx,10^x1m,psym=-4
  ;oplot,fx[0:44],10^x1m[0:44],psym=-4
  ;oplot,fx[72:-1],10^x1m[72:-1],psym=-4
  x1ms=x1m
  for i=1,3 do begin
    w=where(nrec eq i)
    ;oplot,fx(w),x1m(w),psym=-4,color=250
    res=poly_fit(alog10(fx(w)),x1m(w),5,yfit=yfit)
    x1ms(w)=yfit
    ;oplot,fx(w),x1ms(w)
  endfor
  
  w=where(nrec eq 4)
  ww=[indgen(17),indgen(17)+19,37]
  ;oplot,fx(w(ww)),x1m(w(ww)),psym=-4,color=250
  y=interpol(x1m(w(ww)),fx(w(ww)),fx(w))
  if keyword_set(smoothing) then x1ms(w)=smooth(y,3,/edge_trunc)
  ;oplot,fx(w),x1ms(w)
  oplot,fx,10^x1ms,col=250,thick=5
  ;oplot,fx[0:44],10^x1ms[0:44],col=250,thick=5
  ;oplot,fx[72:-1],10^x1ms[72:-1],col=250,thick=5
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
  
  wfreq=where(fx eq 447.750)
  if not(keyword_set(PS)) then window,1,xs=1000,ys=400
  plot_oo,fx,c1m,yra=[1e-22,1e-17],/xsty,/ysty,psym=-4,xtit='Frequency (kHz)',ytit=TexToIdl('Flux density @ 1 AU (W/m^2/Hz)'),tit='Cassini-RPWS/Voyager-PRA 50% and 1% levels',/nodata,charsize=1.75
  oplot,f01,10^m01,psym=-4
  oplot,f50,10^m50,psym=-4
  oplot,fx[wfreq:-1],10^c92[wfreq:-1],psym=-1
  oplot,fx,10^c1m,col=250,thick=5.
  oplot,fx,10^cmm,col=250,thick=5.
  
  xmms=10.^xmms & x1ms=10.^x1ms
  cmm=10.^cmm   & c1m=10.^c1m
  if not(keyword_set(PS)) then window,2,xs=1000,ys=400
  plot_oo,fx,x1ms,/xsty,yra=[1.e-23,1.e-16],/ysty,xtit='Frequency (kHz)',ytit=TexToIdl('Flux density @ 1 AU (W/m^2/Hz)'),tit='Juno/Waves & Cassini/Voyager 50% and 1% levels',/nodata,charsize=1.75
  oplot,fx,x1ms
  ;oplot,fx[0:44],x1ms[0:44]
  ;oplot,fx[72:-1],x1ms[72:-1]
  oplot,fx,c1m,line=2
  ;oplot,fx[0:44],c1m[0:44],line=2
  ;oplot,fx[72:-1],c1m[72:-1],line=2
  oplot,fx,xmms,color=250
  ;oplot,fx[0:44],xmms[0:44],color=250
  ;oplot,fx[72:-1],xmms[72:-1],color=250
  ;oplot,fx[0:44],cmm[0:44],color=250,line=2
  oplot,fx,cmm,color=250,line=2
  ;oplot,fx[72:-1],cmm[72:-1],color=250,line=2
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
  loadct,13
  
  ; gains 50% et 1%
  g50=cmm/xmms & g01=c1m/x1ms
  if not(keyword_set(PS)) then window,3,xs=1000,ys=400
  g01_tmp=g01
  g50_tmp=g50
  ;g01_tmp[45:71]=0.
  ;g50_tmp[45:71]=0.
  minerror=min([g01_tmp,g50_tmp])
  ;minerror=min([g01_tmp[where(g01_tmp ne 0.)],g50_tmp[where(g50_tmp ne 0.)]])
  maxerror=max([g01_tmp,g50_tmp])
  plot_oo,fx,g01,/xsty,yrange=[minerror,maxerror],xtit='Frequency (kHz)',ytit='Conversion factor',tit='Juno Waves conversion factor',/nodata,charsize=1.75
  oplot,fx,g01
  ;oplot,fx[0:44],g01[0:44]
  ;oplot,fx[72:-1],g01[72:-1]
  oplot,fx,g50,color=250
  ;oplot,fx[0:44],g50[0:44],color=250
  ;oplot,fx[72:-1],g50[72:-1],color=250
  ; comparaisons, moyenne géométrique
  g=sqrt(g50*g01)

  ; # modification April 2021: to taking into account that Juno is obsvering the trapped continuum <5 kHz
  ; # and that Cassini did not, the gain <5.00490 kHz (Channel #30 of Juno/Waves) is equal to the one at 5.00490 kHz 
  wfreq_lt_5=where(fx lt 5)
  gain_with_trapped_continuum=g
  g[wfreq_lt_5]=g[wfreq_lt_5[-1]+1]
  g[wfreq_lt_5]=g[wfreq_lt_5[-1]+1]
  

  loadct,51
  oplot,fx[0:wfreq_lt_5[-1]],gain_with_trapped_continuum[0:wfreq_lt_5[-1]],col=240,linestyle=2,thick=5
  oplot,fx,g,thick=5,col=240
  ;oplot,fx[0:44],g[0:44],thick=5
  ;oplot,fx[72:-1],g[72:-1],thick=5
  loadct,53
      oplot,[fx[27],fx[27]],[1e-4,1e4],col=120,linestyle=2
      oplot,[fx[45],fx[45]],[1e-4,1e4],col=120,linestyle=2
      oplot,[3000,3000],[1e-4,1e4],col=120,linestyle=2
   ;   oplot,[fx[71],fx[71]],[1e-4,1e4],col=120,linestyle=2
   ;   oplot,[fx[72],fx[72]],[1e-4,1e4],col=120,linestyle=2
      al_legend,'LFR-Low',position=[2.5,1e4*0.8],box=0,textcolors=120
      al_legend,'LFR-High',position=[25,1e4*0.8],box=0,textcolors=120
      al_legend,'HFR-Low',position=[350,1e4*0.8],box=0,textcolors=120
      al_legend,'HFR-High',position=[6000,1e4*0.8],box=0,textcolors=120
  loadct,13
  oplot,[1.,1.e5],[1,1],line=1
  
  
  
  
  
 ; plot_oo,fx,g01/g50,line=2.,yrange=[min([g/g01,g/g50]),max([g/g01,g/g50])],/xstyle,ytit='Estimated error on conversion factor',xtit='Frequency (kHz)',tit='Mean error = '+strtrim(mean(g/g50),2)+' %, Max error = x'+strtrim(max(g/g50),2),/nodata,charsize=1.75
 ; oplot,[1.,1.e5],[1,1],line=1
 ; oplot,fx,g/g01,line=1,thick=5
 ; oplot,fx,g/g50,line=1,thick=5,color=250
 ; loadct,53
 ;     oplot,[fx[27],fx[27]],[0.1,10],col=120,linestyle=2
 ;     oplot,[fx[45],fx[45]],[0.1,10],col=120,linestyle=2
 ;     oplot,[fx[71],fx[71]],[0.1,10],col=120,linestyle=2
 ;     oplot,[fx[72],fx[72]],[1e-24,1e-16],col=120,linestyle=2
 ;     al_legend,'LFR-Low',position=[2.5,10*0.9],box=0,textcolors=120
 ;     al_legend,'LFR-High',position=[25,10*0.9],box=0,textcolors=120
 ;     al_legend,'HFR-Low',position=[350,10*0.9],box=0,textcolors=120
 ;     al_legend,'HFR-High',position=[6000,10*0.9],box=0,textcolors=120
 ; loadct,13
  print,'Mean error =',mean(g/g50),mean(g01/g)
  print,'Stdev error =',stddev(g/g50),stddev(g01/g)
  print,'Min error =',min(g/g50),min(g01/g)
  print,'Max error =',max(g/g50),max(g01/g)
 
  freq=fx & gain=g
  save,filename='../gain/gain_juno_casvg_'+data_type_treat+version_name+'.sav',freq,g01,g50,gain,gain_with_trapped_continuum
  if keyword_set(PS) then begin
    endps,filename='pdf/intercal_'+data_type_treat
  endif
  


return
end
