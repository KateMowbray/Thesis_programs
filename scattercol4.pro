; km: Program to create scatter plot of the final looptop/initial energy ratio 
;     against the final looptop/initial magnetic field strength ratio 

FUNCTION MAGNITUDE3, arr
  mag = (arr[0]^2+arr[1]^2+arr[2]^2)^0.5
  RETURN, mag
END

;FUNCTION LOOPTOP1, arr ; Finds information about the first loop top pass (mostly irrelevant)
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

FUNCTION LOOPTOP2, arr ; km: Finds index in array of the final loop top pass
  n=long(N_ELEMENTS(arr))
  q=long(2)
  m=long(1)
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

;;c2=SCATTERPLOT(fltarr(2),fltarr(2),SYMBOL='*',AXIS_STYLE=0,LAYOUT=[2,5,1])
; km: some of these commented lines only relevant if you want to show multiple scatter plots in
;     one image

;for i = 0, 1 do begin
;  for j = 0, 4 do begin
;    xval=i
;    zval=j


ee1=fltarr(121) ; arrays storing the energies, currently set up for at most 5 data sets
ee2=fltarr(121)
ee3=fltarr(121)
ee4=fltarr(121)
ee5=fltarr(121)
 
bb1=fltarr(121)
bb2=fltarr(121)
bb3=fltarr(121)
bb4=fltarr(121)
bb5=fltarr(121)

dir='Data/DataR/fermipaperdata/2D_data/2dscatteriny/'

for k=0,120 do begin ; set kmax equal to size of the data set minus 1
  print, k
; need to set a label for each data set as energy values are read in
  ds=getrdata(k+1,/ke,wkdir=dir)
;  dt=getrdata(k+1,/ke,wkdir="Data/DataR/fermipaperdata/2dscatterbz0_005iny/")
;  du=getrdata(k+1,/ke,wkdir="Data/DataR/fermipaperdata/2dscatterbz0_01iny/")
;;  dt=getrdata(k+1,/ke,wkdir="Data/DataR/fermipaperdata/3dscatteryh1/")
;;  du=getrdata(k+1,/ke,wkdir="Data/DataR/fermipaperdata/2dscatterbz0_01iny/")
;;  dv=getrdata(k+1,/ke,wkdir="Data/DataR/ay_one_temp/")
;;  dw=getrdata(k+1,/ke,wkdir="Data/DataR/ay_0_5_yh_1/")
;  dw=getrdata(k+1,/ke,wkdir="Data/minusguidedata/")

;  a1=LOOPTOP1(ds.x)
;  b1=LOOPTOP1(dt.x)
;  c1=LOOPTOP1(du.x)
;  d1=LOOPTOP1(dv.x)
;  e1=LOOPTOP1(dw.x)

; The index for the final loop top pass is found for each data set (can lead to errors if particle
; escapes without ever crossing the loop top, good point to check if an error is returned)
  a2=LOOPTOP2(ds.y)
;  b2=LOOPTOP2(dt.y)
;  c2=LOOPTOP2(du.y)
;;  d2=LOOPTOP2(dv.y)
;  e2=LOOPTOP2(dw.y)

;;print, ds.x[0], ds.y[0], ds.z[0]; 1+zval+5*k+500*xval; "k = ", k, max(ds.ke)/5500
;;  ee1[k] = ds.y[a2]/10000000  
; ee vectors are populated with final to initial energy ratios
  ee1[k] = ds.ke[a2]/ds.ke[0]
;  ee2[k] = dt.ke[b2]/dt.ke[0]
;  ee3[k] = du.ke[c2]/du.ke[0]
;;  ee4[k] = dv.ke[d2]/dv.ke[0]
;  ee5[k] = dw.ke[e2]/ds.ke[0]

; need to take the data sets again to find magnetic field information
; when switching data set the folder in wkdir needs to be changed both here and when information
; on the particle energies is read in. This lower line can be easy to forget about!
  ds=getrdata(k+1,/fields,wkdir=dir)
;  dt=getrdata(k+1,/fields,wkdir="Data/DataR/fermipaperdata/2dscatterbz0_005iny/")
;  du=getrdata(k+1,/fields,wkdir="Data/DataR/fermipaperdata/2dscatterbz0_01iny/")
;;  dt=getrdata(k+1,/fields,wkdir="Data/DataR/fermipaperdata/3dscatteryh1/")
;;  du=getrdata(k+1,/fields,wkdir="Data/DataR/fermipaperdata/2dscatterbz0_01iny/")
;;  dv=getrdata(k+1,/fields,wkdir="Data/DataR/ay_one_temp/")
;;  dw=getrdata(k+1,/fields,wkdir="Data/DataR/ay_0_5_yh_1/")
;  dw=getrdata(k+1,/fields,wkdir="Data/minusguidedata/")

; Ratio of final to initial field is taken
  bb1[k] = MAGNITUDE3(ds.B[*,a2])/MAGNITUDE3(ds.B[*,0])
;  bb2[k] = MAGNITUDE3(dt.B[*,b2])/MAGNITUDE3(dt.B[*,0])
;  bb3[k] = MAGNITUDE3(du.B[*,c2])/MAGNITUDE3(du.B[*,0])
;;  bb4[k] = MAGNITUDE3(dv.B[*,d2])/MAGNITUDE3(dv.B[*,0])
;  bb5[k] = MAGNITUDE3(dw.B[*,e2])/MAGNITUDE3(dw.B[*,0])


endfor

; We want to make sure all the data fits on the chart as we overplot each set of points
xbound=max([max(bb1),max(bb2),max(bb3),max(bb4),max(bb5)])
ybound=max([max(ee1),max(ee2),max(ee3),max(ee4),max(ee5)])



; The plots are constructed, comment lines in/out depending on the number of data sets compared
myplot=SCATTERPLOT(xr=[0,xbound*1.05],yr=[0,ybound*1.05], $
bb1,ee1, SYMBOL='*', font_size=11, $
        xtit='final/initial field strength ratio', $
        ytit='final/initial energy ratio', SYM_COLOR='k',NAME='2D ($B_z = 0$)')
;myplot2=SCATTERPLOT(bb2,ee2,SYMBOL='*',SYM_COLOR='r',OVERPLOT=1,NAME='$B_z = 0.005$T') ; r
;myplot3=SCATTERPLOT(bb3,ee3,SYMBOL='*',SYM_COLOR='b',OVERPLOT=1,NAME='$B_z=0.01$T') ;b
;;myplot=SCATTERPLOT(bb4,ee4,SYMBOL='*',SYM_COLOR='c',OVERPLOT=1)
;myplot=SCATTERPLOT(bb5,ee5,SYMBOL='*',SYM_COLOR='r',OVERPLOT=1)
;leg=LEGEND(TARGET=[myplot,myplot2,myplot3],POSITION=[50,3],/DATA)
;;myplot.Save, "cmtscatter.png"

;;bbb=FLTARR(11,11)
;;eee=FLTARR(11,11)
;;for i = 0,10 do begin
;;  for j = 0,10 do begin
;;    bbb[i,j]=bb1[11*i+j]
;;    eee[i,j]=ee1[11*i+j]
;;  endfor
;;endfor

;;c2=SCATTERPLOT(bb1,ee1,SYMBOL='*',AXIS_STYLE=0,LAYOUT=[2,5,1+(1-i)+2*(4-j)],/CURRENT)

;endfor
;endfor

;c2.Save, 'jetcentremultiscatter.png'

END
