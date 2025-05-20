; Program to bin energies and give approx energy spectrum

kmax=294011

dataarray=DBLARR(kmax)
dataarray2=DBLARR(kmax)
dataarray3=DBLARR(kmax)

openr,1, 'fermiguesscopy1.dat',/f77
readu,1,dataarray

openr,2, 'fermiguesscopy2.dat',/f77
readu,2,dataarray2

openr,3, 'fermiguesscopy3.dat',/f77
readu,3,dataarray3

cc = FLTARR(121)
cc2=FLTARR(121)
cc3=FLTARR(121)
ee = 1.1+FINDGEN(121)*2.2
;dir = 'Data/DataR/jetpaperdata/datasweep/'
dir = 'Data/DataR/jet_thesis_data/standard_case/'
nmax=294011L; 169 ; 11280

FOR i = 1L, nmax DO BEGIN
  print, i
;  ds=getrdata(i,/ke,wkdir=dir)
;  n = SIZE(ds.t, /N_ELEMENTS)
  ke = dataarray(i-1L) ;ds.ke(n-1)
;;  IF (ds.t(n-1) gt 10) then begin
    m = ke/2200
    IF (m gt 120) THEN begin
      print, 'energy exceeds max'
      print, 'Ek = ', ke
    ENDIF
    cc(m) = cc(m) + 1
    m =dataarray2(i-1L)/2200
    cc2(m) = cc2(m) + 1
    m =dataarray3(i-1L)/2200
    cc3(m) = cc3(m) + 1
;;  ENDIF
ENDFOR

myplot=PLOT(ee,cc,xtit='energy (keV)',ytit='count')
myplot=PLOT(ee,cc2,color='r',overplot=1)
;;myplot=PLOT(ee,cc3,color='b',overplot=1)


close,1
close,2
close,3

END
