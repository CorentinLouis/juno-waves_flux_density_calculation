;---------------------------------------------------------------------------------------
  pro MAKE_BACKGROUND, data,save_set, b,s,n, SUPMIN=SUPMIN, POSITIVE=POSITIVE, NSIG=NSIG
;---------------------------------------------------------------------------------------
; Calculation of Background and Sigma arrays of a 2-D distribution of 
; intensities, frequency by frequency. A save set can be put in 'save_set'

; data =  (INPUT)  2-D array of intensities, function of (time,freq)
; save_set =  (INPUT)  name of save_set containing b & s
; b,s,n =  (OUTPUT)  arrays of background, 1 sigma fluctuations, 
;                    and number of pixels used per frequency
; SUPMIN = only values >min(data) are used
; POSITIVE = retains only data values > 0
; NSIG = number of sigmas

  nf=n_elements(data(0,*))
  b=fltarr(nf) & s=b & n=b
  if keyword_set(SUPMIN) then mindata=min(data)
  for i=0,nf-1 do begin
    if keyword_set(SUPMIN) then begin
      test=where(data(*,i) gt mindata)
      if test(0) ne -1 then begin
	       BACKGROUND,data(test,i),bb,ss,nn,POSITIVE=POSITIVE,NSIG=NSIG
	       b(i)=bb & s(i)=ss & n(i)=nn
      endif
    endif else begin
      BACKGROUND,data(*,i),bb,ss,nn,POSITIVE=POSITIVE,NSIG=NSIG
      b(i)=bb & s(i)=ss & n(i)=nn
    endelse
  endfor
  if save_set ne '' then save,b,s,n,filename=save_set
return
end
