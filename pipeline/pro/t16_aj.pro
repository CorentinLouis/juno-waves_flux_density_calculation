; -------------------------------------------------
  Function T16_AJ, t16
; -------------------------------------------------
; date conversion T2016.0 -> AAAAJJJ
;		  T2016.0 -> YYYYDDD
; data type MUST BE double precision, scalar or 1D array
; call : aj = T16_AJ(t16)
;	 yd = T16_AJ(t16)

  aa = dindgen(61)+2016.
  deb= double([0,reform(rebin(reform([366,365,365,365],4,1),4,15),60)])
  for i=1L,60 do deb(i)=deb(i)+deb(i-1)

  aj=double(t16)
  for i=0L, n_elements(t16)-1 do begin
    j=t16(i)-deb
    test=where(j ge 1.0)
    aj(i)=aa(test(n_elements(test)-1))*1000.d0+j(test(n_elements(test)-1))
  endfor

return, aj
end
