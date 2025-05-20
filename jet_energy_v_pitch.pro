dir='Data/DataR/jetpaperdata/datasweep/'

y=2
x=0
z=0
ee=fltarr(141)
pp=10+FINDGEN(141)/2.0

for i = 1, 141 do begin
  print, i
  ds=getrdata(x*5*4*141 + y*4*141 + z*141 + i,/ke,wkdir=dir)
  n=SIZE(ds.t,/N_ELEMENTS)
  ee(i-1)=ds.ke(n-1)/1e3 
endfor

c=PLOT(xr=[10,80],pp,ee,xtit='pitch angle (degrees)',ytit='final energy (keV)')
END
