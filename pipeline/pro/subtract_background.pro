;--------------------------------------------------------------------
  pro SUBTRACT_BACKGROUND, b,s,threshold, data
;--------------------------------------------------------------------
; Subtraction of (background + threshold * sigmas) from a 2-D array of 
; data, and setting of resulting negative values to 0.

; b,s =  (INPUT)  arrays of background and 1 sigma fluctuations
; threshold =  (INPUT)  threshold value for the number of sigmas above background
; data =  (INPUT/OUTPUT)  2-D array of raw data (time,freq) -> data-(b+threshold*s)

  for i=0,n_elements(data(*,0))-1 do data(i,*)=(data(i,*)-b-threshold*s) >0
return
end
