; Program which tests the accuracy of the estimation method by comparing
; against full orbit data and returning the number of points with a given error

needdata = 0 ; set to 1 if you need to gather the final data points from full orbit data

npts = 41L*71L*101L

a1 = 110.0
a2 = 550.0
a3 = 1100.0
a4 = 5500.0
a4_5 = 10000.0
a5 = (-110.0)
a6 = (-550.0)
a7 = (-1100.0)
a8 = (-5500.0)
a8_5 = (-10000.0)

r1 = 0.02
r2 = 0.05
r3 = 0.1
r3_5 = 0.25
r4 = 0.5
r5 = (-0.02)
r6 = (-0.05)
r7 = (-0.1)
r8 = (-0.25)

a1count = 0L
a2count = 0L
a3count = 0L
a4count = 0L
a4_5count = 0L
a5count = 0L
a6count = 0L
a7count = 0L
a8count = 0L
a8_5count = 0L

r1count = 0L
r2count = 0L
r3count = 0L
r3_5count = 0L
r4count = 0L
r5count = 0L
r6count = 0L
r7count = 0L
r8count = 0L

a1arr = fltarr(1)
a2arr = fltarr(1)
a3arr = fltarr(1)
a4arr = fltarr(1)
a4_5arr = fltarr(1)
a5arr = fltarr(1)
a6arr = fltarr(1)
a7arr = fltarr(1)
a8arr = fltarr(1)
a8_5arr = fltarr(1)

r1arr = fltarr(1)
r2arr = fltarr(1)
r3arr = fltarr(1)
r4arr = fltarr(1)
r5arr = fltarr(1)
r6arr = fltarr(1)
r7arr = fltarr(1)
r8arr = fltarr(1)

if (needdata eq 1) then begin
 fullorbitdata = fltarr(npts)
 for i = 0L, npts-1L do begin
   print, i
   ds=getrdata(i+1L,/ke, wkdir='Data/DataR/hugedataset/')
   N = size(ds.ke, /N_ELEMENTS)
   fullorbitdata(i) = ds.ke(N-1L)
 endfor
endif

dataarray = dblarr(41L*71L*101L) ;dblarr(npts)
openr,1,'fermiguess.dat',/f77
readu,1,dataarray

for j = 0L, npts-1L do begin
  abserr = dataarray(j)-fullorbitdata(j)
  relerr = abserr/fullorbitdata(j)

  if (abserr gt a1) then begin
    a1count = a1count + 1L
    a1arr = [a1arr,[j]]
    if (abserr gt a2) then begin
      a2count = a2count + 1L
      a2arr = [a2arr,[j]]
      if (abserr gt a3) then begin
        a3count = a3count + 1L
        a3arr = [a3arr,[j]]
        if (abserr gt a4) then begin
          a4count = a4count + 1L
          a4arr = [a4arr,[j]] 
          if (abserr gt a4_5) then begin
            a4_5count = a4_5count + 1L 
            a4_5arr = [a4_5arr,[j]] 
          endif 
        endif
      endif
    endif
  endif
  
  if (abserr lt a5) then begin
    a5count = a5count + 1L
    a5arr = [a5arr,[j]]
    if (abserr lt a6) then begin
      a6count = a6count + 1L
      a6arr = [a6arr,[j]]
      if (abserr lt a7) then begin
        a7count = a7count + 1L
        a7arr = [a7arr,[j]]
        if (abserr lt a8) then begin
          a8count = a8count + 1L
          a8arr = [a8arr,[j]]
          if (abserr lt a8_5) then begin
            a8_5count = a8_5count + 1L
            a8_5arr = [a8_5arr,[j]]
          endif
        endif
      endif
    endif
  endif

  if (relerr gt r1) then begin
    r1count = r1count + 1L
    r1arr = [r1arr,[j]]
    if (relerr gt r2) then begin
      r2count = r2count + 1L
      r2arr = [r2arr,[j]]
      if (relerr gt r3) then begin
        r3count = r3count + 1L
        r3arr = [r3arr,[j]]
          if (relerr gt r3_5) then begin
             r3_5count = r3_5count + 1L
             if (relerr gt r4) then begin
               r4count = r4count + 1L
               r4arr = [r4arr,[j]]
             endif
        endif
      endif
    endif
  endif
  
  if (relerr lt r5) then begin
    r5count = r5count + 1L
    r5arr = [r5arr,[j]]
    if (relerr lt r6) then begin
      r6count = r6count + 1L
      r6arr = [r6arr,[j]]
      if (relerr lt r7) then begin
        r7count = r7count + 1L
        r7arr = [r7arr,[j]]
        if (relerr lt r8) then begin
          r8count = r8count + 1L
          r8arr = [r8arr,[j]]
        endif
      endif
    endif
  endif
  
endfor

print, 'errors for ', npts, ' particles'
print, '******************************************'
print, a1count, 'have raw errors greater than ', a1
print, a2count, 'have raw errors greater than ', a2
print, a3count, 'have raw errors greater than ', a3
print, a4count, 'have raw errors greater than ', a4
print, a4_5count, 'have raw errors greater than ', a4_5
print, a5count, 'have raw errors less than ', a5
print, a6count, 'have raw errors less than ', a6
print, a7count, 'have raw errors less than ', a7
print, a8count, 'have raw errors less than ', a8
print, a8_5count, 'have raw errors less than ', a8_5
print, '******************************************'
print, r1count, 'have relative errors greater than ', r1
print, r2count, 'have relative errors greater than ', r2
print, r3count, 'have relative errors greater than ', r3
print, r3_5count, 'have relative errors greater than ', r3_5
print, r4count, 'have relative errors greater than ', r4
print, r5count, 'have relative errors less than ', r5
print, r6count, 'have relative errors less than ', r6
print, r7count, 'have relative errors less than ', r7
print, r8count, 'have relative errors less than ', r8

myplot=SCATTERPLOT(findgen(npts), dataarray-fullorbitdata)

close,1

END
