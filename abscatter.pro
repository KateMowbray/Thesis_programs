; Program to create scatter plot of the final looptop/initial energy ratio 
; against the final looptop/initial magnetic field strength ratio 

FUNCTION MAGNITUDE3, arr
  mag = (arr[0]^2+arr[1]^2+arr[2]^2)^0.5
  RETURN, mag
END

;FUNCTION LOOPTOP1, arr
;  n=N_ELEMENTS(arr)
;  q=2
;  m=1
;  WHILE(q GT 1) DO BEGIN
;    IF (arr[m]*arr[m-1] LE 0) THEN BEGIN
;      q=0
;    ENDIF ELSE BEGIN
;      m=m+1
;    ENDELSE
;  ENDWHILE
;  RETURN, m
;END

FUNCTION LOOPTOP2, arr
  n=N_ELEMENTS(arr)
  q=2
  m=1
  WHILE(q GT 1) DO BEGIN
;    IF (arr[n-m]*arr[n-m-1] LT 0) THEN BEGIN ;only works in looptop centred at x=0
    IF ((arr[n-m]-arr[n-m-1])*(arr[n-m-1]-arr[n-m-2]) LT 0) $
       AND (arr[n-m-1]-arr[n-m-2] GT 0) THEN BEGIN
      q=0
    ENDIF ELSE BEGIN
      m=m+1
    ENDELSE
    IF (m GT n) THEN BEGIN
      print, "error for particle ",k+1,", particle makes no loop top passes"
      BREAK
    ENDIF
  ENDWHILE
  RETURN, n-m-1
END


xval=3
zval=3

ee=fltarr(1)
 
bb=fltarr(1)
yy=3.5+FINDGEN(1001)/400.0
earray=fltarr(1001)
barray=fltarr(1001)

ds=getrdata(500*xval+zval+1,/ke);,wkdir="Data/DataR/jetpaperdata/yvar/")
a2=LOOPTOP2(ds.y)
ee[0]=ds.ke[a2]/ds.ke[0]
earray[0]=ds.ke[a2]/ds.ke[0]
ds=getrdata(500*xval+zval+1,/fields);,wkdir="Data/DataR/jetpaperdata/yvar/")
bb[0] = MAGNITUDE3(ds.B[*,a2])/MAGNITUDE3(ds.B[*,0])
barray[0] = MAGNITUDE3(ds.B[*,a2])/MAGNITUDE3(ds.B[*,0])

print, bb[0], ee[0]

timecolor=(1.0-max(ds.t)/100.0)*[255,255,255]

myplot=SCATTERPLOT(xr=[0,14],yr=[0,14],bb,ee, SYMBOL='*', $
        xtit='final/initial field strength ratio', $
        ytit='final/initial energy ratio', SYM_COLOR=timecolor)

for k=1,99 do begin 

  ds=getrdata(500*xval+5*k+zval+1,/ke);,wkdir="Data/DataR/jetpaperdata/yvar/")

  a2=LOOPTOP2(ds.y)

  print, "k = ", k
  
  ee[0] = ds.ke[a2]/ds.ke[0]

  earray[k] = ds.ke[a2]/ds.ke[0]

  ds=getrdata(500*xval+5*k+zval+1,/fields);,wkdir="Data/DataR/jetpaperdata/yvar/")

  bb[0] = MAGNITUDE3(ds.B[*,a2])/MAGNITUDE3(ds.B[*,0])

  barray[k] = MAGNITUDE3(ds.B[*,a2])/MAGNITUDE3(ds.B[*,0])

  timecolor=(1.0-max(ds.t)/100.0)*[255,255,255]

  myplot=SCATTERPLOT(bb,ee, SYMBOL='*',SYM_COLOR=timecolor,OVERPLOT=1)

  print, "bb = ", bb, "ee = ", ee

endfor

;;myplot=SCATTERPLOT(xr=[0,max(bb1)*1.05],yr=[0,max(ee1)*1.05],bb1,ee1, SYMBOL='*', $
;;        xtit='final/initial field strength ratio', $
;;        ytit='final/initial energy ratio', SYM_COLOR='k')
;;myplot=SCATTERPLOT(bb2,ee2,SYMBOL='*',SYM_COLOR='r',OVERPLOT=1)
;;myplot=SCATTERPLOT(bb3,ee3,SYMBOL='*',SYM_COLOR='r',OVERPLOT=1) ;m
;;myplot=SCATTERPLOT(bb4,ee4,SYMBOL='*',SYM_COLOR='c',OVERPLOT=1)
;myplot=SCATTERPLOT(bb5,ee5,SYMBOL='*',SYM_COLOR='r',OVERPLOT=1)

;myplot.Save, "abpfbigt0scatter.png"

;;bbb=FLTARR(11,11)
;;eee=FLTARR(11,11)
;;for i = 0,10 do begin
;;  for j = 0,10 do begin
;;    bbb[i,j]=bb1[11*i+j]
;;    eee[i,j]=ee1[11*i+j]
;;  endfor
;;endfor

END
