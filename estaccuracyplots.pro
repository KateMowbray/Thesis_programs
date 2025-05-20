a1pitch10x = fltarr(1)
a1pitch10y = fltarr(1)
a1pitch30x = fltarr(1)
a1pitch30y = fltarr(1)
a1pitch50x = fltarr(1)
a1pitch50y = fltarr(1)
a1pitch70x = fltarr(1)
a1pitch70y = fltarr(1)

a1size = size(a1arr, /N_ELEMENTS)

if (a1size ne 0) then begin
  for p = 0L, a1size-1L do begin
    element = a1arr(p)
    xval = -5.0e6 + floor(element/(41.0*71.0)) * 5.0e4
    yval = 1.0e7 + 5.0e5 * floor((element MOD 2911L)/(41.0))
    pitch = floor(element MOD 41L)

    if (0 le pitch) and (pitch lt 10) then begin
      a1pitch10x = [a1pitch10x, [xval]]
      a1pitch10y = [a1pitch10y, [yval]]
    endif

    if (10 le pitch) and (pitch lt 20) then begin
      a1pitch30x = [a1pitch30x, [xval]]
      a1pitch30y = [a1pitch30y, [yval]]
    endif

    if (20 le pitch) and (pitch lt 30) then begin
      a1pitch50x = [a1pitch50x, [xval]]
      a1pitch50y = [a1pitch50y, [yval]]
    endif

    if (30 le pitch) and (pitch le 40) then begin
      a1pitch70x = [a1pitch70x, [xval]]
      a1pitch70y = [a1pitch70y, [yval]]
    endif
  
  endfor
endif

a2pitch10x = fltarr(1)
a2pitch10y = fltarr(1)
a2pitch30x = fltarr(1)
a2pitch30y = fltarr(1)
a2pitch50x = fltarr(1)
a2pitch50y = fltarr(1)
a2pitch70x = fltarr(1)
a2pitch70y = fltarr(1)

a2size = size(a2arr, /N_ELEMENTS)

if (a2size ne 0) then begin
  for p = 0L, a2size-1L do begin
    element = a2arr(p)
    xval = -5.0e6 + floor(element/(41.0*71.0)) * 5.0e4
    yval = 1.0e7 + 5.0e5 * floor((element MOD 2911L)/(41.0))
    pitch = floor(element MOD 41L)

    if (0 le pitch) and (pitch lt 10) then begin
      a2pitch10x = [a2pitch10x, [xval]]
      a2pitch10y = [a2pitch10y, [yval]]
    endif

    if (10 le pitch) and (pitch lt 20) then begin
      a2pitch30x = [a2pitch30x, [xval]]
      a2pitch30y = [a2pitch30y, [yval]]
    endif

    if (20 le pitch) and (pitch lt 30) then begin
      a2pitch50x = [a2pitch50x, [xval]]
      a2pitch50y = [a2pitch50y, [yval]]
    endif

    if (30 le pitch) and (pitch le 40) then begin
      a2pitch70x = [a2pitch70x, [xval]]
      a2pitch70y = [a2pitch70y, [yval]]
    endif
  
  endfor
endif

a3pitch10x = fltarr(1)
a3pitch10y = fltarr(1)
a3pitch30x = fltarr(1)
a3pitch30y = fltarr(1)
a3pitch50x = fltarr(1)
a3pitch50y = fltarr(1)
a3pitch70x = fltarr(1)
a3pitch70y = fltarr(1)

a3size = size(a3arr, /N_ELEMENTS)

if (a3size ne 0) then begin
  for p = 0L, a3size-1L do begin
    element = a3arr(p)
    xval = -5.0e6 + floor(element/(41.0*71.0)) * 5.0e4
    yval = 1.0e7 + 5.0e5 * floor((element MOD 2911L)/(41.0))
    pitch = floor(element MOD 41L)

    if (0 le pitch) and (pitch lt 10) then begin
      a3pitch10x = [a3pitch10x, [xval]]
      a3pitch10y = [a3pitch10y, [yval]]
    endif

    if (10 le pitch) and (pitch lt 20) then begin
      a3pitch30x = [a3pitch30x, [xval]]
      a3pitch30y = [a3pitch30y, [yval]]
    endif

    if (20 le pitch) and (pitch lt 30) then begin
      a3pitch50x = [a3pitch50x, [xval]]
      a3pitch50y = [a3pitch50y, [yval]]
    endif

    if (30 le pitch) and (pitch le 40) then begin
      a3pitch70x = [a3pitch70x, [xval]]
      a3pitch70y = [a3pitch70y, [yval]]
    endif
  
  endfor
endif

a4pitch10x = fltarr(1)
a4pitch10y = fltarr(1)
a4pitch30x = fltarr(1)
a4pitch30y = fltarr(1)
a4pitch50x = fltarr(1)
a4pitch50y = fltarr(1)
a4pitch70x = fltarr(1)
a4pitch70y = fltarr(1)

a4size = size(a4arr, /N_ELEMENTS)

if (a4size ne 0) then begin
  for p = 0L, a4size-1L do begin
    element = a4arr(p)
    xval = -5.0e6 + floor(element/(41.0*71.0)) * 5.0e4
    yval = 1.0e7 + 5.0e5 * floor((element MOD 2911L)/(41.0))
    pitch = floor(element MOD 41L)

    if (0 le pitch) and (pitch lt 10) then begin
      a4pitch10x = [a4pitch10x, [xval]]
      a4pitch10y = [a4pitch10y, [yval]]
    endif

    if (10 le pitch) and (pitch lt 20) then begin
      a4pitch30x = [a4pitch30x, [xval]]
      a4pitch30y = [a4pitch30y, [yval]]
    endif

    if (20 le pitch) and (pitch lt 30) then begin
      a4pitch50x = [a4pitch50x, [xval]]
      a4pitch50y = [a4pitch50y, [yval]]
    endif

    if (30 le pitch) and (pitch le 40) then begin
      a4pitch70x = [a4pitch70x, [xval]]
      a4pitch70y = [a4pitch70y, [yval]]
    endif
  
  endfor
endif

a4_5pitch10x = fltarr(1)
a4_5pitch10y = fltarr(1)
a4_5pitch30x = fltarr(1)
a4_5pitch30y = fltarr(1)
a4_5pitch50x = fltarr(1)
a4_5pitch50y = fltarr(1)
a4_5pitch70x = fltarr(1)
a4_5pitch70y = fltarr(1)

a4_5size = size(a4_5arr, /N_ELEMENTS)

if (a4_5size ne 0) then begin
  for p = 0L, a4_5size-1L do begin
    element = a4_5arr(p)
    xval = -5.0e6 + floor(element/(41.0*71.0)) * 5.0e4
    yval = 1.0e7 + 5.0e5 * floor((element MOD 2911L)/(41.0))
    pitch = floor(element MOD 41L)

    if (0 le pitch) and (pitch lt 10) then begin
      a4_5pitch10x = [a4_5pitch10x, [xval]]
      a4_5pitch10y = [a4_5pitch10y, [yval]]
    endif

    if (10 le pitch) and (pitch lt 20) then begin
      a4_5pitch30x = [a4_5pitch30x, [xval]]
      a4_5pitch30y = [a4_5pitch30y, [yval]]
    endif

    if (20 le pitch) and (pitch lt 30) then begin
      a4_5pitch50x = [a4_5pitch50x, [xval]]
      a4_5pitch50y = [a4_5pitch50y, [yval]]
    endif

    if (30 le pitch) and (pitch le 40) then begin
      a4_5pitch70x = [a4_5pitch70x, [xval]]
      a4_5pitch70y = [a4_5pitch70y, [yval]]
    endif
  
  endfor
endif

a5pitch10x = fltarr(1)
a5pitch10y = fltarr(1)
a5pitch30x = fltarr(1)
a5pitch30y = fltarr(1)
a5pitch50x = fltarr(1)
a5pitch50y = fltarr(1)
a5pitch70x = fltarr(1)
a5pitch70y = fltarr(1)

a5size = size(a5arr, /N_ELEMENTS)

if (a5size ne 0) then begin
  for p = 0L, a5size-1L do begin
    element = a5arr(p)
    xval = -5.0e6 + floor(element/(41.0*71.0)) * 5.0e4
    yval = 1.0e7 + 5.0e5 * floor((element MOD 2911L)/(41.0))
    pitch = floor(element MOD 41L)

    if (0 le pitch) and (pitch lt 10) then begin
      a5pitch10x = [a5pitch10x, [xval]]
      a5pitch10y = [a5pitch10y, [yval]]
    endif

    if (10 le pitch) and (pitch lt 20) then begin
      a5pitch30x = [a5pitch30x, [xval]]
      a5pitch30y = [a5pitch30y, [yval]]
    endif

    if (20 le pitch) and (pitch lt 30) then begin
      a5pitch50x = [a5pitch50x, [xval]]
      a5pitch50y = [a5pitch50y, [yval]]
    endif

    if (30 le pitch) and (pitch le 40) then begin
      a5pitch70x = [a5pitch70x, [xval]]
      a5pitch70y = [a5pitch70y, [yval]]
    endif
  
  endfor
endif

a6pitch10x = fltarr(1)
a6pitch10y = fltarr(1)
a6pitch30x = fltarr(1)
a6pitch30y = fltarr(1)
a6pitch50x = fltarr(1)
a6pitch50y = fltarr(1)
a6pitch70x = fltarr(1)
a6pitch70y = fltarr(1)

a6size = size(a6arr, /N_ELEMENTS)

if (a6size ne 0) then begin
  for p = 0L, a6size-1L do begin
    element = a6arr(p)
    xval = -5.0e6 + floor(element/(41.0*71.0)) * 5.0e4
    yval = 1.0e7 + 5.0e5 * floor((element MOD 2911L)/(41.0))
    pitch = floor(element MOD 41L)

    if (0 le pitch) and (pitch lt 10) then begin
      a6pitch10x = [a6pitch10x, [xval]]
      a6pitch10y = [a6pitch10y, [yval]]
    endif

    if (10 le pitch) and (pitch lt 20) then begin
      a6pitch30x = [a6pitch30x, [xval]]
      a6pitch30y = [a6pitch30y, [yval]]
    endif

    if (20 le pitch) and (pitch lt 30) then begin
      a6pitch50x = [a6pitch50x, [xval]]
      a6pitch50y = [a6pitch50y, [yval]]
    endif

    if (30 le pitch) and (pitch le 40) then begin
      a6pitch70x = [a6pitch70x, [xval]]
      a6pitch70y = [a6pitch70y, [yval]]
    endif
  
  endfor
endif

a7pitch10x = fltarr(1)
a7pitch10y = fltarr(1)
a7pitch30x = fltarr(1)
a7pitch30y = fltarr(1)
a7pitch50x = fltarr(1)
a7pitch50y = fltarr(1)
a7pitch70x = fltarr(1)
a7pitch70y = fltarr(1)

a7size = size(a7arr, /N_ELEMENTS)

if (a7size ne 0) then begin
  for p = 0L, a7size-1L do begin
    element = a7arr(p)
    xval = -5.0e6 + floor(element/(41.0*71.0)) * 5.0e4
    yval = 1.0e7 + 5.0e5 * floor((element MOD 2911L)/(41.0))
    pitch = floor(element MOD 41L)

    if (0 le pitch) and (pitch lt 10) then begin
      a7pitch10x = [a7pitch10x, [xval]]
      a7pitch10y = [a7pitch10y, [yval]]
    endif

    if (10 le pitch) and (pitch lt 20) then begin
      a7pitch30x = [a7pitch30x, [xval]]
      a7pitch30y = [a7pitch30y, [yval]]
    endif

    if (20 le pitch) and (pitch lt 30) then begin
      a7pitch50x = [a7pitch50x, [xval]]
      a7pitch50y = [a7pitch50y, [yval]]
    endif

    if (30 le pitch) and (pitch le 40) then begin
      a7pitch70x = [a7pitch70x, [xval]]
      a7pitch70y = [a7pitch70y, [yval]]
    endif
  
  endfor
endif

a8pitch10x = fltarr(1)
a8pitch10y = fltarr(1)
a8pitch30x = fltarr(1)
a8pitch30y = fltarr(1)
a8pitch50x = fltarr(1)
a8pitch50y = fltarr(1)
a8pitch70x = fltarr(1)
a8pitch70y = fltarr(1)

a8size = size(a8arr, /N_ELEMENTS)

if (a8size ne 0) then begin
  for p = 0L, a8size-1L do begin
    element = a8arr(p)
    xval = -5.0e6 + floor(element/(41.0*71.0)) * 5.0e4
    yval = 1.0e7 + 5.0e5 * floor((element MOD 2911L)/(41.0))
    pitch = floor(element MOD 41L)

    if (0 le pitch) and (pitch lt 10) then begin
      a8pitch10x = [a8pitch10x, [xval]]
      a8pitch10y = [a8pitch10y, [yval]]
    endif

    if (10 le pitch) and (pitch lt 20) then begin
      a8pitch30x = [a8pitch30x, [xval]]
      a8pitch30y = [a8pitch30y, [yval]]
    endif

    if (20 le pitch) and (pitch lt 30) then begin
      a8pitch50x = [a8pitch50x, [xval]]
      a8pitch50y = [a8pitch50y, [yval]]
    endif

    if (30 le pitch) and (pitch le 40) then begin
      a8pitch70x = [a8pitch70x, [xval]]
      a8pitch70y = [a8pitch70y, [yval]]
    endif
  
  endfor
endif

a8_5pitch10x = fltarr(1)
a8_5pitch10y = fltarr(1)
a8_5pitch30x = fltarr(1)
a8_5pitch30y = fltarr(1)
a8_5pitch50x = fltarr(1)
a8_5pitch50y = fltarr(1)
a8_5pitch70x = fltarr(1)
a8_5pitch70y = fltarr(1)

a8_5size = size(a8_5arr, /N_ELEMENTS)

if (a8_5size ne 0) then begin
  for p = 0L, a8_5size-1L do begin
    element = a8_5arr(p)
    xval = -5.0e6 + floor(element/(41.0*71.0)) * 5.0e4
    yval = 1.0e7 + 5.0e5 * floor((element MOD 2911L)/(41.0))
    pitch = floor(element MOD 41L)

    if (0 le pitch) and (pitch lt 10) then begin
      a8_5pitch10x = [a8_5pitch10x, [xval]]
      a8_5pitch10y = [a8_5pitch10y, [yval]]
    endif

    if (10 le pitch) and (pitch lt 20) then begin
      a8_5pitch30x = [a8_5pitch30x, [xval]]
      a8_5pitch30y = [a8_5pitch30y, [yval]]
    endif

    if (20 le pitch) and (pitch lt 30) then begin
      a8_5pitch50x = [a8_5pitch50x, [xval]]
      a8_5pitch50y = [a8_5pitch50y, [yval]]
    endif

    if (30 le pitch) and (pitch le 40) then begin
      a8_5pitch70x = [a8_5pitch70x, [xval]]
      a8_5pitch70y = [a8_5pitch70y, [yval]]
    endif
  
  endfor
endif









r1pitch10x = fltarr(1)
r1pitch10y = fltarr(1)
r1pitch30x = fltarr(1)
r1pitch30y = fltarr(1)
r1pitch50x = fltarr(1)
r1pitch50y = fltarr(1)
r1pitch70x = fltarr(1)
r1pitch70y = fltarr(1)

r1size = size(r1arr, /N_ELEMENTS)

if (r1size ne 0) then begin
  for p = 0L, r1size-1L do begin
    element = r1arr(p)
    xval = -5.0e6 + floor(element/(41.0*71.0)) * 5.0e4
    yval = 1.0e7 + 5.0e5 * floor((element MOD 2911L)/(41.0))
    pitch = floor(element MOD 41L)

    if (0 le pitch) and (pitch lt 10) then begin
      r1pitch10x = [r1pitch10x, [xval]]
      r1pitch10y = [r1pitch10y, [yval]]
    endif

    if (10 le pitch) and (pitch lt 20) then begin
      r1pitch30x = [r1pitch30x, [xval]]
      r1pitch30y = [r1pitch30y, [yval]]
    endif

    if (20 le pitch) and (pitch lt 30) then begin
      r1pitch50x = [r1pitch50x, [xval]]
      r1pitch50y = [r1pitch50y, [yval]]
    endif

    if (30 le pitch) and (pitch le 40) then begin
      r1pitch70x = [r1pitch70x, [xval]]
      r1pitch70y = [r1pitch70y, [yval]]
    endif
  
  endfor
endif

r2pitch10x = fltarr(1)
r2pitch10y = fltarr(1)
r2pitch30x = fltarr(1)
r2pitch30y = fltarr(1)
r2pitch50x = fltarr(1)
r2pitch50y = fltarr(1)
r2pitch70x = fltarr(1)
r2pitch70y = fltarr(1)

r2size = size(r2arr, /N_ELEMENTS)

if (r2size ne 0) then begin
  for p = 0L, r2size-1L do begin
    element = r2arr(p)
    xval = -5.0e6 + floor(element/(41.0*71.0)) * 5.0e4
    yval = 1.0e7 + 5.0e5 * floor((element MOD 2911L)/(41.0))
    pitch = floor(element MOD 41L)

    if (0 le pitch) and (pitch lt 10) then begin
      r2pitch10x = [r2pitch10x, [xval]]
      r2pitch10y = [r2pitch10y, [yval]]
    endif

    if (10 le pitch) and (pitch lt 20) then begin
      r2pitch30x = [r2pitch30x, [xval]]
      r2pitch30y = [r2pitch30y, [yval]]
    endif

    if (20 le pitch) and (pitch lt 30) then begin
      r2pitch50x = [r2pitch50x, [xval]]
      r2pitch50y = [r2pitch50y, [yval]]
    endif

    if (30 le pitch) and (pitch le 40) then begin
      r2pitch70x = [r2pitch70x, [xval]]
      r2pitch70y = [r2pitch70y, [yval]]
    endif
  
  endfor
endif

r3pitch10x = fltarr(1)
r3pitch10y = fltarr(1)
r3pitch30x = fltarr(1)
r3pitch30y = fltarr(1)
r3pitch50x = fltarr(1)
r3pitch50y = fltarr(1)
r3pitch70x = fltarr(1)
r3pitch70y = fltarr(1)

r3size = size(r3arr, /N_ELEMENTS)

if (r3size ne 0) then begin
  for p = 0L, r3size-1L do begin
    element = r3arr(p)
    xval = -5.0e6 + floor(element/(41.0*71.0)) * 5.0e4
    yval = 1.0e7 + 5.0e5 * floor((element MOD 2911L)/(41.0))
    pitch = floor(element MOD 41L)

    if (0 le pitch) and (pitch lt 10) then begin
      r3pitch10x = [r3pitch10x, [xval]]
      r3pitch10y = [r3pitch10y, [yval]]
    endif

    if (10 le pitch) and (pitch lt 20) then begin
      r3pitch30x = [r3pitch30x, [xval]]
      r3pitch30y = [r3pitch30y, [yval]]
    endif

    if (20 le pitch) and (pitch lt 30) then begin
      r3pitch50x = [r3pitch50x, [xval]]
      r3pitch50y = [r3pitch50y, [yval]]
    endif

    if (30 le pitch) and (pitch le 40) then begin
      r3pitch70x = [r3pitch70x, [xval]]
      r3pitch70y = [r3pitch70y, [yval]]
    endif
  
  endfor
endif

r4pitch10x = fltarr(1)
r4pitch10y = fltarr(1)
r4pitch30x = fltarr(1)
r4pitch30y = fltarr(1)
r4pitch50x = fltarr(1)
r4pitch50y = fltarr(1)
r4pitch70x = fltarr(1)
r4pitch70y = fltarr(1)

r4size = size(r4arr, /N_ELEMENTS)

if (r4size ne 0) then begin
  for p = 0L, r4size-1L do begin
    element = r4arr(p)
    xval = -5.0e6 + floor(element/(41.0*71.0)) * 5.0e4
    yval = 1.0e7 + 5.0e5 * floor((element MOD 2911L)/(41.0))
    pitch = floor(element MOD 41L)

    if (0 le pitch) and (pitch lt 10) then begin
      r4pitch10x = [r4pitch10x, [xval]]
      r4pitch10y = [r4pitch10y, [yval]]
    endif

    if (10 le pitch) and (pitch lt 20) then begin
      r4pitch30x = [r4pitch30x, [xval]]
      r4pitch30y = [r4pitch30y, [yval]]
    endif

    if (20 le pitch) and (pitch lt 30) then begin
      r4pitch50x = [r4pitch50x, [xval]]
      r4pitch50y = [r4pitch50y, [yval]]
    endif

    if (30 le pitch) and (pitch le 40) then begin
      r4pitch70x = [r4pitch70x, [xval]]
      r4pitch70y = [r4pitch70y, [yval]]
    endif
  
  endfor
endif

r5pitch10x = fltarr(1)
r5pitch10y = fltarr(1)
r5pitch30x = fltarr(1)
r5pitch30y = fltarr(1)
r5pitch50x = fltarr(1)
r5pitch50y = fltarr(1)
r5pitch70x = fltarr(1)
r5pitch70y = fltarr(1)

r5size = size(r5arr, /N_ELEMENTS)

if (r5size ne 0) then begin
  for p = 0L, r5size-1L do begin
    element = r5arr(p)
    xval = -5.0e6 + floor(element/(41.0*71.0)) * 5.0e4
    yval = 1.0e7 + 5.0e5 * floor((element MOD 2911L)/(41.0))
    pitch = floor(element MOD 41L)

    if (0 le pitch) and (pitch lt 10) then begin
      r5pitch10x = [r5pitch10x, [xval]]
      r5pitch10y = [r5pitch10y, [yval]]
    endif

    if (10 le pitch) and (pitch lt 20) then begin
      r5pitch30x = [r5pitch30x, [xval]]
      r5pitch30y = [r5pitch30y, [yval]]
    endif

    if (20 le pitch) and (pitch lt 30) then begin
      r5pitch50x = [r5pitch50x, [xval]]
      r5pitch50y = [r5pitch50y, [yval]]
    endif

    if (30 le pitch) and (pitch le 40) then begin
      r5pitch70x = [r5pitch70x, [xval]]
      r5pitch70y = [r5pitch70y, [yval]]
    endif
  
  endfor
endif

r6pitch10x = fltarr(1)
r6pitch10y = fltarr(1)
r6pitch30x = fltarr(1)
r6pitch30y = fltarr(1)
r6pitch50x = fltarr(1)
r6pitch50y = fltarr(1)
r6pitch70x = fltarr(1)
r6pitch70y = fltarr(1)

r6size = size(r6arr, /N_ELEMENTS)

if (r6size ne 0) then begin
  for p = 0L, r6size-1L do begin
    element = r6arr(p)
    xval = -5.0e6 + floor(element/(41.0*71.0)) * 5.0e4
    yval = 1.0e7 + 5.0e5 * floor((element MOD 2911L)/(41.0))
    pitch = floor(element MOD 41L)

    if (0 le pitch) and (pitch lt 10) then begin
      r6pitch10x = [r6pitch10x, [xval]]
      r6pitch10y = [r6pitch10y, [yval]]
    endif

    if (10 le pitch) and (pitch lt 20) then begin
      r6pitch30x = [r6pitch30x, [xval]]
      r6pitch30y = [r6pitch30y, [yval]]
    endif

    if (20 le pitch) and (pitch lt 30) then begin
      r6pitch50x = [r6pitch50x, [xval]]
      r6pitch50y = [r6pitch50y, [yval]]
    endif

    if (30 le pitch) and (pitch le 40) then begin
      r6pitch70x = [r6pitch70x, [xval]]
      r6pitch70y = [r6pitch70y, [yval]]
    endif
  
  endfor
endif

r7pitch10x = fltarr(1)
r7pitch10y = fltarr(1)
r7pitch30x = fltarr(1)
r7pitch30y = fltarr(1)
r7pitch50x = fltarr(1)
r7pitch50y = fltarr(1)
r7pitch70x = fltarr(1)
r7pitch70y = fltarr(1)

r7size = size(r7arr, /N_ELEMENTS)

if (r7size ne 0) then begin
  for p = 0L, r7size-1L do begin
    element = r7arr(p)
    xval = -5.0e6 + floor(element/(41.0*71.0)) * 5.0e4
    yval = 1.0e7 + 5.0e5 * floor((element MOD 2911L)/(41.0))
    pitch = floor(element MOD 41L)

    if (0 le pitch) and (pitch lt 10) then begin
      r7pitch10x = [r7pitch10x, [xval]]
      r7pitch10y = [r7pitch10y, [yval]]
    endif

    if (10 le pitch) and (pitch lt 20) then begin
      r7pitch30x = [r7pitch30x, [xval]]
      r7pitch30y = [r7pitch30y, [yval]]
    endif

    if (20 le pitch) and (pitch lt 30) then begin
      r7pitch50x = [r7pitch50x, [xval]]
      r7pitch50y = [r7pitch50y, [yval]]
    endif

    if (30 le pitch) and (pitch le 40) then begin
      r7pitch70x = [r7pitch70x, [xval]]
      r7pitch70y = [r7pitch70y, [yval]]
    endif
  
  endfor
endif

r8pitch10x = fltarr(1)
r8pitch10y = fltarr(1)
r8pitch30x = fltarr(1)
r8pitch30y = fltarr(1)
r8pitch50x = fltarr(1)
r8pitch50y = fltarr(1)
r8pitch70x = fltarr(1)
r8pitch70y = fltarr(1)

r8size = size(r8arr, /N_ELEMENTS)

if (r8size ne 0) then begin
  for p = 0L, r8size-1L do begin
    element = r8arr(p)
    xval = -5.0e6 + floor(element/(41.0*71.0)) * 5.0e4
    yval = 1.0e7 + 5.0e5 * floor((element MOD 2911L)/(41.0))
    pitch = floor(element MOD 41L)

    if (0 le pitch) and (pitch lt 10) then begin
      r8pitch10x = [r8pitch10x, [xval]]
      r8pitch10y = [r8pitch10y, [yval]]
    endif

    if (10 le pitch) and (pitch lt 20) then begin
      r8pitch30x = [r8pitch30x, [xval]]
      r8pitch30y = [r8pitch30y, [yval]]
    endif

    if (20 le pitch) and (pitch lt 30) then begin
      r8pitch50x = [r8pitch50x, [xval]]
      r8pitch50y = [r8pitch50y, [yval]]
    endif

    if (30 le pitch) and (pitch le 40) then begin
      r8pitch70x = [r8pitch70x, [xval]]
      r8pitch70y = [r8pitch70y, [yval]]
    endif
  
  endfor
endif




myplot1 = SCATTERPLOT(xr=[-5.2e6,0.2e6], yr=[0.9e7,4.6e7], a3pitch10x, a3pitch10y, SYMBOL='*', $
                      SYM_COLOR = [255,192,192], SYM_SIZE = 0.7, tit = 'absolute 10-28 degrees')
myplot1 = SCATTERPLOT(a4pitch10x, a4pitch10y, SYMBOL='*', $
                      SYM_COLOR = [255,128,128], SYM_SIZE = 0.9, OVERPLOT = 1)
;myplot1 = SCATTERPLOT(a3pitch10x, a3pitch10y, SYMBOL='*', $
;                      SYM_COLOR = [255,64,64], SYM_SIZE = 1.1, OVERPLOT = 1)
myplot1 = SCATTERPLOT(a4_5pitch10x, a4_5pitch10y, SYMBOL='*', $
                      SYM_COLOR = [255,0,0], SYM_SIZE = 1.3, OVERPLOT = 1)

myplot2 = SCATTERPLOT(xr=[-5.2e6,0.2e6], yr=[0.9e7,4.6e7], a3pitch30x, a3pitch30y, SYMBOL='*', $
                      SYM_COLOR = [255,192,192], SYM_SIZE = 0.7, tit = 'absolute 30-48 degrees')
myplot2 = SCATTERPLOT(a4pitch30x, a4pitch30y, SYMBOL='*', $
                      SYM_COLOR = [255,128,128], SYM_SIZE = 0.9, OVERPLOT = 1)
;myplot2 = SCATTERPLOT(a3pitch30x, a3pitch30y, SYMBOL='*', $
;                      SYM_COLOR = [255,64,64], SYM_SIZE = 1.1, OVERPLOT = 1)
myplot2 = SCATTERPLOT(a4_5pitch30x, a4_5pitch30y, SYMBOL='*', $
                      SYM_COLOR = [255,0,0], SYM_SIZE = 1.3, OVERPLOT = 1)

myplot3 = SCATTERPLOT(xr=[-5.2e6,0.2e6], yr=[0.9e7,4.6e7], a3pitch50x, a3pitch50y, SYMBOL='*', $
                      SYM_COLOR = [255,192,192], SYM_SIZE = 0.7, tit = 'absolute 50-68 degrees')
myplot3 = SCATTERPLOT(a4pitch50x, a4pitch50y, SYMBOL='*', $
                      SYM_COLOR = [255,128,128], SYM_SIZE = 0.9, OVERPLOT = 1)
;myplot3 = SCATTERPLOT(a3pitch50x, a3pitch50y, SYMBOL='*', $
;                      SYM_COLOR = [255,64,64], SYM_SIZE = 1.1, OVERPLOT = 1)
myplot3 = SCATTERPLOT(a4_5pitch50x, a4_5pitch50y, SYMBOL='*', $
                      SYM_COLOR = [255,0,0], SYM_SIZE = 1.3, OVERPLOT = 1)

myplot4 = SCATTERPLOT(xr=[-5.2e6,0.2e6], yr=[0.9e7,4.6e7], a3pitch70x, a3pitch70y, SYMBOL='*', $
                      SYM_COLOR = [255,192,192], SYM_SIZE = 0.7, tit = 'absolute 70-90 degrees')
myplot4 = SCATTERPLOT(a4pitch70x, a4pitch70y, SYMBOL='*', $
                      SYM_COLOR = [255,128,128], SYM_SIZE = 0.9, OVERPLOT = 1)
;myplot4 = SCATTERPLOT(a3pitch70x, a3pitch70y, SYMBOL='*', $
;                      SYM_COLOR = [255,64,64], SYM_SIZE = 1.1, OVERPLOT = 1)
myplot4 = SCATTERPLOT(a4_5pitch70x, a4_5pitch70y, SYMBOL='*', $
                      SYM_COLOR = [255,0,0], SYM_SIZE = 1.3, OVERPLOT = 1)


myplot1.Save, 'newcase6esterrabs10deg.png'
myplot2.Save, 'newcase6esterrabs30deg.png'
myplot3.Save, 'newcase6esterrabs50deg.png'
myplot4.Save, 'newcase6esterrabs70deg.png'


myplot5 = SCATTERPLOT(xr=[-5.2e6,0.2e6], yr=[0.9e7,4.6e7], a7pitch10x, a7pitch10y, SYMBOL='*', $
                      SYM_COLOR = [192,192,255], SYM_SIZE = 0.7, tit = 'absolute 10-28 degrees')
myplot5 = SCATTERPLOT(a8pitch10x, a8pitch10y, SYMBOL='*', $
                      SYM_COLOR = [128,128,255], SYM_SIZE = 0.9, OVERPLOT = 1)
;myplot5 = SCATTERPLOT(a7pitch10x, a7pitch10y, SYMBOL='*', $
;                      SYM_COLOR = [64,64,255], SYM_SIZE = 1.1, OVERPLOT = 1)
myplot5 = SCATTERPLOT(a8_5pitch10x, a8_5pitch10y, SYMBOL='*', $
                      SYM_COLOR = [0,0,255], SYM_SIZE = 1.3, OVERPLOT = 1)

myplot6 = SCATTERPLOT(xr=[-5.2e6,0.2e6], yr=[0.9e7,4.6e7], a7pitch30x, a7pitch30y, SYMBOL='*', $
                      SYM_COLOR = [192,192,255], SYM_SIZE = 0.7, tit = 'absolute 30-48 degrees')
myplot6 = SCATTERPLOT(a8pitch30x, a8pitch30y, SYMBOL='*', $
                      SYM_COLOR = [128,128,255], SYM_SIZE = 0.9, OVERPLOT = 1)
;myplot6 = SCATTERPLOT(a7pitch30x, a7pitch30y, SYMBOL='*', $
;                      SYM_COLOR = [64,64,255], SYM_SIZE = 1.1, OVERPLOT = 1)
myplot6 = SCATTERPLOT(a8_5pitch30x, a8_5pitch30y, SYMBOL='*', $
                      SYM_COLOR = [0,0,255], SYM_SIZE = 1.3, OVERPLOT = 1)

myplot7 = SCATTERPLOT(xr=[-5.2e6,0.2e6], yr=[0.9e7,4.6e7], a7pitch50x, a7pitch50y, SYMBOL='*', $
                      SYM_COLOR = [192,192,255], SYM_SIZE = 0.7, tit = 'absolute 50-68 degrees')
myplot7 = SCATTERPLOT(a8pitch50x, a8pitch50y, SYMBOL='*', $
                      SYM_COLOR = [128,128,255], SYM_SIZE = 0.9, OVERPLOT = 1)
;myplot7 = SCATTERPLOT(a7pitch50x, a7pitch50y, SYMBOL='*', $
;                      SYM_COLOR = [64,64,255], SYM_SIZE = 1.1, OVERPLOT = 1)
myplot7 = SCATTERPLOT(a8_5pitch50x, a8_5pitch50y, SYMBOL='*', $
                      SYM_COLOR = [0,0,255], SYM_SIZE = 1.3, OVERPLOT = 1)

myplot8 = SCATTERPLOT(xr=[-5.2e6,0.2e6], yr=[0.9e7,4.6e7], a7pitch70x, a7pitch70y, SYMBOL='*', $
                      SYM_COLOR = [192,192,255], SYM_SIZE = 0.7, tit = 'absolute 70-90 degrees')
myplot8 = SCATTERPLOT(a8pitch70x, a8pitch70y, SYMBOL='*', $
                      SYM_COLOR = [128,128,255], SYM_SIZE = 0.9, OVERPLOT = 1)
;myplot8 = SCATTERPLOT(a7pitch70x, a7pitch70y, SYMBOL='*', $
;                      SYM_COLOR = [64,64,255], SYM_SIZE = 1.1, OVERPLOT = 1)
myplot8 = SCATTERPLOT(a8_5pitch70x, a8_5pitch70y, SYMBOL='*', $
                      SYM_COLOR = [0,0,255], SYM_SIZE = 1.3, OVERPLOT = 1)

myplot5.Save, 'newcase6esterrabsunder10deg.png'
myplot6.Save, 'newcase6esterrabsunder30deg.png'
myplot7.Save, 'newcase6esterrabsunder50deg.png'
myplot8.Save, 'newcase6esterrabsunder70deg.png'



myplot9 = SCATTERPLOT(xr=[-5.2e6,0.2e6], yr=[0.9e7,4.6e7], r1pitch10x, r1pitch10y, SYMBOL='*', $
                      SYM_COLOR = [255,192,192], SYM_SIZE = 0.7, tit = 'relative 10-28 degrees')
myplot9 = SCATTERPLOT(r2pitch10x, r2pitch10y, SYMBOL='*', $
                      SYM_COLOR = [255,128,128], SYM_SIZE = 0.9, OVERPLOT = 1)
myplot9 = SCATTERPLOT(r3pitch10x, r3pitch10y, SYMBOL='*', $
                      SYM_COLOR = [255,64,64], SYM_SIZE = 1.1, OVERPLOT = 1)
myplot9 = SCATTERPLOT(r4pitch10x, r4pitch10y, SYMBOL='*', $
                      SYM_COLOR = [255,0,0], SYM_SIZE = 1.3, OVERPLOT = 1)

myplot10 = SCATTERPLOT(xr=[-5.2e6,0.2e6], yr=[0.9e7,4.6e7], r1pitch30x, r1pitch30y, SYMBOL='*', $
                      SYM_COLOR = [255,192,192], SYM_SIZE = 0.7, tit = 'relative 30-48 degrees')
myplot10 = SCATTERPLOT(r2pitch30x, r2pitch30y, SYMBOL='*', $
                      SYM_COLOR = [255,128,128], SYM_SIZE = 0.9, OVERPLOT = 1)
myplot10 = SCATTERPLOT(r3pitch30x, r3pitch30y, SYMBOL='*', $
                      SYM_COLOR = [255,64,64], SYM_SIZE = 1.1, OVERPLOT = 1)
myplot10 = SCATTERPLOT(r4pitch30x, r4pitch30y, SYMBOL='*', $
                      SYM_COLOR = [255,0,0], SYM_SIZE = 1.3, OVERPLOT = 1)

myplot11 = SCATTERPLOT(xr=[-5.2e6,0.2e6], yr=[0.9e7,4.6e7], r1pitch50x, r1pitch50y, SYMBOL='*', $
                      SYM_COLOR = [255,192,192], SYM_SIZE = 0.7, tit = 'relative 50-68 degrees')
myplot11 = SCATTERPLOT(r2pitch50x, r2pitch50y, SYMBOL='*', $
                      SYM_COLOR = [255,128,128], SYM_SIZE = 0.9, OVERPLOT = 1)
myplot11 = SCATTERPLOT(r3pitch50x, r3pitch50y, SYMBOL='*', $
                      SYM_COLOR = [255,64,64], SYM_SIZE = 1.1, OVERPLOT = 1)
myplot11 = SCATTERPLOT(r4pitch50x, r4pitch50y, SYMBOL='*', $
                      SYM_COLOR = [255,0,0], SYM_SIZE = 1.3, OVERPLOT = 1)

myplot12 = SCATTERPLOT(xr=[-5.2e6,0.2e6], yr=[0.9e7,4.6e7], r1pitch70x, r1pitch70y, SYMBOL='*', $
                      SYM_COLOR = [255,192,192], SYM_SIZE = 0.7, tit = 'relative 70-90 degrees')
myplot12 = SCATTERPLOT(r2pitch70x, r2pitch70y, SYMBOL='*', $
                      SYM_COLOR = [255,128,128], SYM_SIZE = 0.9, OVERPLOT = 1)
myplot12 = SCATTERPLOT(r3pitch70x, r3pitch70y, SYMBOL='*', $
                      SYM_COLOR = [255,64,64], SYM_SIZE = 1.1, OVERPLOT = 1)
myplot12 = SCATTERPLOT(r4pitch70x, r4pitch70y, SYMBOL='*', $
                      SYM_COLOR = [255,0,0], SYM_SIZE = 1.3, OVERPLOT = 1)

;myplot9.Save, 'esterrrel10deg6.png'
;myplot10.Save, 'esterrrel30deg6.png'
;myplot11.Save, 'esterrrel50deg6.png'
;myplot12.Save, 'esterrrel70deg6.png'



myplot13 = SCATTERPLOT(xr=[-5.2e6,0.2e6], yr=[0.9e7,4.6e7], r5pitch10x, r5pitch10y, SYMBOL='*', $
                      SYM_COLOR = [192,192,255], SYM_SIZE = 0.7, tit = 'relative 10-28 degrees')
myplot13 = SCATTERPLOT(r6pitch10x, r6pitch10y, SYMBOL='*', $
                      SYM_COLOR = [128,128,255], SYM_SIZE = 0.9, OVERPLOT = 1)
myplot13 = SCATTERPLOT(r7pitch10x, r7pitch10y, SYMBOL='*', $
                      SYM_COLOR = [64,64,255], SYM_SIZE = 1.1, OVERPLOT = 1)
myplot13 = SCATTERPLOT(r8pitch10x, r8pitch10y, SYMBOL='*', $
                      SYM_COLOR = [0,0,255], SYM_SIZE = 1.3, OVERPLOT = 1)

myplot14 = SCATTERPLOT(xr=[-5.2e6,0.2e6], yr=[0.9e7,4.6e7], r5pitch30x, r5pitch30y, SYMBOL='*', $
                      SYM_COLOR = [192,192,255], SYM_SIZE = 0.7, tit = 'relative 30-48 degrees')
myplot14 = SCATTERPLOT(r6pitch30x, r6pitch30y, SYMBOL='*', $
                      SYM_COLOR = [128,128,255], SYM_SIZE = 0.9, OVERPLOT = 1)
myplot14 = SCATTERPLOT(r7pitch30x, r7pitch30y, SYMBOL='*', $
                      SYM_COLOR = [64,64,255], SYM_SIZE = 1.1, OVERPLOT = 1)
myplot14 = SCATTERPLOT(r8pitch30x, r8pitch30y, SYMBOL='*', $
                      SYM_COLOR = [0,0,255], SYM_SIZE = 1.3, OVERPLOT = 1)

myplot15 = SCATTERPLOT(xr=[-5.2e6,0.2e6], yr=[0.9e7,4.6e7], r5pitch50x, r5pitch50y, SYMBOL='*', $
                      SYM_COLOR = [192,192,255], SYM_SIZE = 0.7, tit = 'relative 50-68 degrees')
myplot15 = SCATTERPLOT(r6pitch50x, r6pitch50y, SYMBOL='*', $
                      SYM_COLOR = [128,128,255], SYM_SIZE = 0.9, OVERPLOT = 1)
myplot15 = SCATTERPLOT(r7pitch50x, r7pitch50y, SYMBOL='*', $
                      SYM_COLOR = [64,64,255], SYM_SIZE = 1.1, OVERPLOT = 1)
myplot15 = SCATTERPLOT(r8pitch50x, r8pitch50y, SYMBOL='*', $
                      SYM_COLOR = [0,0,255], SYM_SIZE = 1.3, OVERPLOT = 1)

myplot16 = SCATTERPLOT(xr=[-5.2e6,0.2e6], yr=[0.9e7,4.6e7], r5pitch70x, r5pitch70y, SYMBOL='*', $
                      SYM_COLOR = [192,192,255], SYM_SIZE = 0.7, tit = 'relative 70-90 degrees')
myplot16 = SCATTERPLOT(r6pitch70x, r6pitch70y, SYMBOL='*', $
                      SYM_COLOR = [128,128,255], SYM_SIZE = 0.9, OVERPLOT = 1)
myplot16 = SCATTERPLOT(r7pitch70x, r7pitch70y, SYMBOL='*', $
                      SYM_COLOR = [64,64,255], SYM_SIZE = 1.1, OVERPLOT = 1)
myplot16 = SCATTERPLOT(r8pitch70x, r8pitch70y, SYMBOL='*', $
                      SYM_COLOR = [0,0,255], SYM_SIZE = 1.3, OVERPLOT = 1)


;myplot13.Save, 'esterrrelunder10deg6.png'
;myplot14.Save, 'esterrrelunder30deg6.png'
;myplot15.Save, 'esterrrelunder50deg6.png'
;myplot16.Save, 'esterrrelunder70deg6.png'



myplot17 = SCATTERPLOT(xr=[-5.2e6,0.2e6], yr=[0.9e7,4.6e7], a4pitch10x, a4pitch10y, SYMBOL='*', $
                      SYM_COLOR = [255,128,128], SYM_SIZE = 1.1, tit = 'absolute 10-28 degrees (large errors)')
myplot17 = SCATTERPLOT(a4_5pitch10x, a4_5pitch10y, SYMBOL='*', $
                      SYM_COLOR = [255,0,0], SYM_SIZE = 1.3, OVERPLOT = 1)

myplot18 = SCATTERPLOT(xr=[-5.2e6,0.2e6], yr=[0.9e7,4.6e7], a4pitch30x, a4pitch30y, SYMBOL='*', $
                      SYM_COLOR = [255,128,128], SYM_SIZE = 1.1, tit = 'absolute 30-48 degrees (large errors)')
myplot18 = SCATTERPLOT(a4_5pitch30x, a4_5pitch30y, SYMBOL='*', $
                      SYM_COLOR = [255,0,0], SYM_SIZE = 1.3, OVERPLOT = 1)

myplot19 = SCATTERPLOT(xr=[-5.2e6,0.2e6], yr=[0.9e7,4.6e7], a4pitch50x, a4pitch50y, SYMBOL='*', $
                      SYM_COLOR = [255,128,128], SYM_SIZE = 1.1, tit = 'absolute 50-68 degrees (large errors)')
myplot19 = SCATTERPLOT(a4_5pitch50x, a4_5pitch50y, SYMBOL='*', $
                      SYM_COLOR = [255,0,0], SYM_SIZE = 1.3, OVERPLOT = 1)

myplot20 = SCATTERPLOT(xr=[-5.2e6,0.2e6], yr=[0.9e7,4.6e7], a4pitch70x, a4pitch70y, SYMBOL='*', $
                      SYM_COLOR = [255,128,128], SYM_SIZE = 1.1, tit = 'absolute 70-90 degrees (large errors)')
myplot20 = SCATTERPLOT(a4_5pitch70x, a4_5pitch70y, SYMBOL='*', $
                      SYM_COLOR = [255,0,0], SYM_SIZE = 1.3, OVERPLOT = 1)


myplot17.Save, 'newcase6esterrabskev10deg.png'
myplot18.Save, 'newcase6esterrabskev30deg.png'
myplot19.Save, 'newcase6esterrabskev50deg.png'
myplot20.Save, 'newcase6esterrabskev70deg.png'


myplot21 = SCATTERPLOT(xr=[-5.2e6,0.2e6], yr=[0.9e7,4.6e7], a8pitch10x, a8pitch10y, SYMBOL='*', $
                      SYM_COLOR = [128,128,255], SYM_SIZE = 1.1, tit = 'absolute 10-28 degrees (large errors)')
myplot21 = SCATTERPLOT(a8_5pitch10x, a8_5pitch10y, SYMBOL='*', $
                      SYM_COLOR = [0,0,255], SYM_SIZE = 1.3, OVERPLOT = 1)

myplot22 = SCATTERPLOT(xr=[-5.2e6,0.2e6], yr=[0.9e7,4.6e7], a8pitch30x, a8pitch30y, SYMBOL='*', $
                      SYM_COLOR = [128,128,255], SYM_SIZE = 1.1, tit = 'absolute 30-48 degrees (large errors)')
myplot22 = SCATTERPLOT(a8_5pitch30x, a8_5pitch30y, SYMBOL='*', $
                      SYM_COLOR = [0,0,255], SYM_SIZE = 1.3, OVERPLOT = 1)

myplot23 = SCATTERPLOT(xr=[-5.2e6,0.2e6], yr=[0.9e7,4.6e7], a8pitch50x, a8pitch50y, SYMBOL='*', $
                      SYM_COLOR = [128,128,255], SYM_SIZE = 1.1, tit = 'absolute 50-68 degrees (large errors)')
myplot23 = SCATTERPLOT(a8_5pitch50x, a8_5pitch50y, SYMBOL='*', $
                      SYM_COLOR = [0,0,255], SYM_SIZE = 1.3, OVERPLOT = 1)

myplot24 = SCATTERPLOT(xr=[-5.2e6,0.2e6], yr=[0.9e7,4.6e7], a8pitch70x, a8pitch70y, SYMBOL='*', $
                      SYM_COLOR = [128,128,255], SYM_SIZE = 1.1, tit = 'absolute 70-90 degrees (large errors)')
myplot24 = SCATTERPLOT(a8_5pitch70x, a8_5pitch70y, SYMBOL='*', $
                      SYM_COLOR = [0,0,255], SYM_SIZE = 1.3, OVERPLOT = 1)

myplot21.Save, 'newcase6esterrabsunderkev10deg.png'
myplot22.Save, 'newcase6esterrabsunderkev30deg.png'
myplot23.Save, 'newcase6esterrabsunderkev50deg.png'
myplot24.Save, 'newcase6esterrabsunderkev70deg.png'



myplot25 = SCATTERPLOT(xr=[-5.2e6,0.2e6], yr=[0.9e7,4.6e7], r3pitch10x, r3pitch10y, SYMBOL='*', $
                      SYM_COLOR = [255,64,64], SYM_SIZE = 1.1, tit = 'relative 10-28 degrees (keV errors)')
myplot25 = SCATTERPLOT(r4pitch10x, r4pitch10y, SYMBOL='*', $
                      SYM_COLOR = [255,0,0], SYM_SIZE = 1.3, OVERPLOT = 1)

myplot26 = SCATTERPLOT(xr=[-5.2e6,0.2e6], yr=[0.9e7,4.6e7], r3pitch30x, r3pitch30y, SYMBOL='*', $
                      SYM_COLOR = [255,64,64], SYM_SIZE = 1.1, tit = 'relative 30-48 degrees (keV errors)')
myplot26 = SCATTERPLOT(r4pitch30x, r4pitch30y, SYMBOL='*', $
                      SYM_COLOR = [255,0,0], SYM_SIZE = 1.3, OVERPLOT = 1)

myplot27 = SCATTERPLOT(xr=[-5.2e6,0.2e6], yr=[0.9e7,4.6e7], r3pitch50x, r3pitch50y, SYMBOL='*', $
                      SYM_COLOR = [255,64,64], SYM_SIZE = 1.1, tit = 'relative 50-68 degrees (keV errors)')
myplot27 = SCATTERPLOT(r4pitch50x, r4pitch50y, SYMBOL='*', $
                      SYM_COLOR = [255,0,0], SYM_SIZE = 1.3, OVERPLOT = 1)

myplot28 = SCATTERPLOT(xr=[-5.2e6,0.2e6], yr=[0.9e7,4.6e7], r3pitch70x, r3pitch70y, SYMBOL='*', $
                      SYM_COLOR = [255,64,64], SYM_SIZE = 1.1, tit = 'relative 70-90 degrees (keV errors)')
myplot28 = SCATTERPLOT(r4pitch70x, r4pitch70y, SYMBOL='*', $
                      SYM_COLOR = [255,0,0], SYM_SIZE = 1.3, OVERPLOT = 1)


;myplot25.Save, 'esterrrelkev10deg6.png'
;myplot26.Save, 'esterrrelkev30deg6.png'
;myplot27.Save, 'esterrrelkev50deg6.png'
;myplot28.Save, 'esterrrelkev70deg6.png'



myplot29 = SCATTERPLOT(xr=[-5.2e6,0.2e6], yr=[0.9e7,4.6e7], r7pitch10x, r7pitch10y, SYMBOL='*', $
                      SYM_COLOR = [64,64,255], SYM_SIZE = 1.1, tit = 'relative 10-28 degrees (keV errors)')
myplot29 = SCATTERPLOT(r8pitch10x, r8pitch10y, SYMBOL='*', $
                      SYM_COLOR = [0,0,255], SYM_SIZE = 1.3, OVERPLOT = 1)

myplot30 = SCATTERPLOT(xr=[-5.2e6,0.2e6], yr=[0.9e7,4.6e7], r7pitch30x, r7pitch30y, SYMBOL='*', $
                      SYM_COLOR = [64,64,255], SYM_SIZE = 1.1, tit = 'relative 30-48 degrees (keV errors)')
myplot30 = SCATTERPLOT(r8pitch30x, r8pitch30y, SYMBOL='*', $
                      SYM_COLOR = [0,0,255], SYM_SIZE = 1.3, OVERPLOT = 1)

myplot31 = SCATTERPLOT(xr=[-5.2e6,0.2e6], yr=[0.9e7,4.6e7], r7pitch50x, r7pitch50y, SYMBOL='*', $
                      SYM_COLOR = [64,64,255], SYM_SIZE = 1.1, tit = 'relative 50-68 degrees (keV errors)')
myplot31 = SCATTERPLOT(r8pitch50x, r8pitch50y, SYMBOL='*', $
                      SYM_COLOR = [0,0,255], SYM_SIZE = 1.3, OVERPLOT = 1)

myplot32 = SCATTERPLOT(xr=[-5.2e6,0.2e6], yr=[0.9e7,4.6e7], r7pitch70x, r7pitch70y, SYMBOL='*', $
                      SYM_COLOR = [64,64,255], SYM_SIZE = 1.1, tit = 'relative 70-90 degrees (keV errors)')
myplot32 = SCATTERPLOT(r8pitch70x, r8pitch70y, SYMBOL='*', $
                      SYM_COLOR = [0,0,255], SYM_SIZE = 1.3, OVERPLOT = 1)

;myplot29.Save, 'esterrrelunderkev10deg6.png'
;myplot30.Save, 'esterrrelunderkev30deg6.png'
;myplot31.Save, 'esterrrelunderkev50deg6.png'
;myplot32.Save, 'esterrrelunderkev70deg6.png'


END
