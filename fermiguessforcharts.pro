FUNCTION MAGNITUDE3, arr
  ; gives magnitude of 3D vector, useful for stuff like modB
  mag = (arr[0]^2+arr[1]^2+arr[2]^2)^0.5
  RETURN, mag
END

FUNCTION LOOPTOP2, arr

  ; Finds final loop top pass for a particle orbit
  n=N_ELEMENTS(arr)
  q=2
  m=1
  WHILE(q GT 1) DO BEGIN
    IF (arr[n-m]*arr[n-m-1] LT 0) THEN BEGIN
      q=0
    ENDIF ELSE BEGIN
      m=m+1
    ENDELSE
  ENDWHILE
  RETURN, n-m
END

kmax=91
dataarray=dblarr(kmax)

openr,1,'newestimatordata.dat',/f77
readu,1,dataarray

dataarray=dataarray/5500.0

bb=fltarr(kmax)
ee=fltarr(kmax)

dir ='Data/DataR/';estimatortestdata/centre_85deg/'
;dir= 'Data/DataR/fermipaperdata/3D_data/3dscatter/'

for k = 1, kmax do begin
  print, k
  ds=getrdata(k,/ke,wkdir=dir)
 a2=looptop2(ds.x)
 n= N_ELEMENTS(ds.x)
 ee[k-1] = ds.ke[n-1]/ds.ke[0]
  dt=getrdata(k,/fields,wkdir=dir)
 bb[k-1] = MAGNITUDE3(dt.B[*,a2])/MAGNITUDE3(dt.B[*,0])
endfor

;myplot=SCATTERPLOT(xtit='final B / initial B', ytit = $
;'final energy / initial energy', bb,ee,SYMBOL='*',SYM_COLOR='k')
y1=max(dataarray)
y2=max(ee)
yscl= max([y1,y2])
myplot=SCATTERPLOT(xr=[1.4,1.05*max(bb)],yr=[0,1.05*yscl],xtit= $
'final field strength/ initial field strength', ytit = $
;;;;'initial y position (10Mm)', ytit= $
'final energy / initial energy', bb,ee,SYMBOL='*',SYM_COLOR='k')

;myplot=SCATTERPLOT(findgen(60),findgen(60)*max(ee)/max(bb), SYM_COLOR='r', OVERPLOT=1)
myplot=SCATTERPLOT(bb,dataarray,SYMBOL='*',SYM_COLOR='r',OVERPLOT=1)

;myplot=PLOT(findgen(kmax),ee-dataarray)
close,1

end
