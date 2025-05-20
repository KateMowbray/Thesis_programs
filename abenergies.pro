xx=FINDGEN(13)-12
zz=FINDGEN(13)-12

xval=0
yval=3
zval=0
ds=getrdata(141*4*5*xval+141*4*yval+141*zval+1,/ke,wkdir='Data/DataR/bigjetxyzpitchdata/')

c=PLOT(ds.t,ds.ke,xr=[0,100],yr=[0,25000],AXIS_STYLE=0,LAYOUT=[13,13,1]);12*13+1])
;c=PLOT(ds.x,ds.y,yr=[0,6e7],AXIS_STYLE=0,LAYOUT=[13,13,12*13+1])

;;;for i=0,12 do begin
;;;  for j=0,12 do begin
;;;    ds=getrdata(13*i+j+1,/ke);,wkdir='Data/DataR/newabdata/5500ev70deg/')
;;;    c=PLOT(ds.t,ds.ke,xr=[0,100],AXIS_STYLE=0,LAYOUT=[13,13,13*(12-i)+j+1],/CURRENT)
;    c=PLOT(ds.x,ds.y,yr=[0,6e7],AXIS_STYLE=0,LAYOUT=[13,13,13*(12-i)+j+1],/CURRENT)
;;    ttt2(i,12-j)=max(dt.t)
;;    print, 13*i+j+1
;;;   endfor
;;;endfor

for k=0,140 do begin
  ds=getrdata(141*4*5*xval+141*4*yval+141*zval+1+k,/ke,wkdir='Data/DataR/bigjetxyzpitchdata/')
  if (k EQ 0) or (k EQ 140) then begin
     print, ds.x[0], ds.y[0], ds.z[0]
  endif
  print, "k = ", k,  ds.ke[max(size(ds.ke))-1]
  c=PLOT(ds.t,ds.ke,xr=[0,100],yr=[0,25000],AXIS_STYLE=0,LAYOUT=[13,13,k+1],/CURRENT)
;  c=PLOT(ds.x,ds.y,AXIS_STYLE=0,LAYOUT=[13,13,k+1],/CURRENT)
endfor
;c.Save, 'abpfbigt0energies.png'

END
