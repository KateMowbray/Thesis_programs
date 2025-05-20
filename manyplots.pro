ds=getrdata(1,/ke);,wkdir='Data/DataR/jet_thesis_data/standard_case/');,wkdir='Data/DataR/newabdata/5500ev70deg/')

;for k=0,12 do begin
;  for k2=0,12 do begin
;    ds=getrdata(13*k+k2+1,/ke);,wkdir='Data/DataR/newabdata/5500ev70deg/')
;    c=PLOT(ds.t, ds.ke,xr=[0:100],AXIS_STYLE=0,LAYOUT=[13,13,13*(12-k)+k2+1],/CURRENT)
;;  dt=getrdata(8*k+6+240,/ke)
;;  c=PLOT(dt.t, dt.ke,xr=[0:100],AXIS_STYLE=0,LAYOUT=[10,9,k+1+30],/CURRENT)
;;  du=getrdata(8*k+6+480,/ke)
;;  c=PLOT(du.t, du.ke,xr=[0:100],AXIS_STYLE=0,LAYOUT=[10,9,k+1+60],/CURRENT)
;;  print, ds.x[0], ds.y[0], ds.z[0]
;   print, 'k=', 13*k+k2+1, max(ds.t), max(ds.ke)/5500
; endfor
;endfor

xval=2
zval=0

for k=0,7 do begin
    r=FLOOR(k/13)
    ds=getrdata(0*81+k+1,/ke);,wkdir='Data/DataR/jet_thesis_data/standard_case/')
;    dt=getrdata(k+1,/ke,wkdir='Data/DataR/jet_thesis_data/t0_5_5/')
    c=PLOT(ds.t, ds.ke, xr=[0,100],AXIS_STYLE=0,LAYOUT=[10,10,k+1],/CURRENT)
;    c=PLOT(dt.t, dt.ke,color='r',OVERPLOT=1,LAYOUT=[13,13,13*(12-k+13*r)+r+1],/CURRENT)
;;  dt=getrdata(8*k+6+240,/ke)
;;  c=PLOT(dt.t, dt.ke,xr=[0:100],AXIS_STYLE=0,LAYOUT=[10,9,k+1+30],/CURRENT)
;;  du=getrdata(8*k+6+480,/ke)
;;  c=PLOT(du.t, du.ke,xr=[0:100],AXIS_STYLE=0,LAYOUT=[10,9,k+1+60],/CURRENT)
;;  c=PLOT(ds.x, ds.y,xr=[-10000000,10000000],yr=[0,67500000],AXIS_STYLE=0,LAYOUT=[10,10,k+1],/CURRENT)
;;  print, max(ds.x)
endfor

END
