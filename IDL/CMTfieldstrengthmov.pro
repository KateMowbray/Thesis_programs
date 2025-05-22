;PRO CMTfieldstrengthmov
; CMTfieldstrengthmov
;  numerically integrates field lines according to field specified in derivs
;  then displays
;  km: adjusted from CMTFLINE to generate frames for a movie showing evolution
;      in time of multiple field lines
;  REQUIRES - ODEINT, derivs, RKQC

@CMTODEINT
common someothername, length, c1, Bzf, ttt

 lscl=1e7
 zscl=0.1
 dims=[1.,1.,1.]
 ;mysymsize=[lscl,lscl,lscl]*0.1*dims
 mysymsize=[1.0,1.0,1.0]*7.*dims  


t0=00.0
t1=7.0
noframes=70
;;timearr=t0+FINDGEN(noframes+1)*(t1-t0)/noframes
;;fieldarr=FLTARR(noframes+1)
steps=t0+FINDGEN(noframes+1)/(noframes)*(t1-t0)

particleson=1 ; Set to 1 to overplot particles on field lines
colorson=0    ; Set to 1 to colour particles according to their energies


;;----- define variables ----------------
eps=0.0000001
Bzf=0.00
h1=0.01     	    ; step size?
string1='derivsab'

;;----- pick initial conditions ------------
ystart=dblarr(3)
r1=dblarr(3)
r2=dblarr(3)
rsteps=intarr(3)


;r1[0]=-0.9129	;xstart
;r2[0]=0.9129 	;xend
r1[0]=-1.200;-1.1343459;-1.200
r2[0]=-1.200;-1.1343459;-1.200
rsteps[0]=1 	;nx
r1[1]=1.25;0.00;1.25	;ystart
r2[1]=1.25;0.00;1.25 	;yend
rsteps[1]=1	;ny
;r1[2]=0.231564  	;zstart
;r2[2]=0.231564   	;zend
r1[2]=-0.50
r2[2]=-0.50
rsteps[2]=1 	;nz


IF (rsteps[0] EQ 1) THEN  gdx=1.0d0 ELSE gdx=1.0d0/FLOAT(rsteps[0]-1)
IF (rsteps[1] EQ 1) THEN  gdy=1.0d0 ELSE gdy=1.0d0/FLOAT(rsteps[1]-1)
IF (rsteps[2] EQ 1) THEN  gdz=1.0d0 ELSE gdz=1.0d0/FLOAT(rsteps[2]-1)
gds=[gdx,gdy,gdz]
lbox=[r2(0)-r1(0),r2(1)-r1(1),r2(2)-r1(2)]

np=rsteps[0]*rsteps[1]*rsteps[2]

flinecol=[0,0,255]

; Set up an array to hold all of the initial foot point values for the given 
; positions. These will then be fed into the field line generator.

fps=dblarr(np,3)

;;;ttt=1.05
ttt=0.00
nparticles=3
particles=FLTARR(nparticles,noframes+1,4)
mu=fltarr(nparticles)


FOR ix=0,rsteps[0]-1 DO BEGIN
 FOR iy=0,rsteps[1]-1 DO BEGIN
  FOR iz=0,rsteps[2]-1 DO BEGIN
   ;gds=[gdx,gdy,gdz]
   i=[0,0,0]
   ;yorig= R1+lbox*(i*1.0D0)*gds
   ystart= R1+lbox*(i*1.0D0)*[gdx,gdy,gdz]
   CMTODEINT, ystart, eps, h1, string1, yp
   fps[ix*rsteps[2]*rsteps[1]+iy*rsteps[2]+iz,*]=[INTERPOL(yp(0,*),yp(1,*),0), $
      0.0,INTERPOL(yp(2,*),yp(1,*),0)]
   if particleson eq 1 then begin
    FOR ip=0,nparticles-1 DO BEGIN
     ds=getrdata(ip+1,/ke,/fields,wkdir="../Data/DataR/")
     particles(ip,*,0)=INTERPOL(ds.x,ds.t,steps,/SPLINE)
     particles(ip,*,1)=INTERPOL(ds.y,ds.t,steps,/SPLINE)
     particles(ip,*,2)=INTERPOL(ds.z,ds.t,steps,/SPLINE)
     particles(ip,*,3)=INTERPOL(ds.ke,ds.t,steps,/SPLINE)
     mu(ip)=ds.e2[0]/(ds.B[0,0]^2+ds.B[1,0]^2+ds.B[2,0]^2)^0.5
;;     fieldarr=SQRT(INTERPOL(ds.B[0,*],ds.t,steps,/SPLINE)^2+INTERPOL(ds.B[1,*],ds.t,steps,/SPLINE)^2+$
;;                   INTERPOL(ds.B[2,*],ds.t,steps,/SPLINE)^2)
    ENDFOR
   endif
  ENDFOR
 ENDFOR
ENDFOR

if colorson eq 1 then begin
  mincol=min(particles(*,*,3))
  maxcol=max(particles(*,*,3))
endif

for number=0,noframes do begin
ttt=(000.0+t0+(t1-t0)*number/noframes)/10.0 
print, ttt

pno=0

FOR ix=0,np-1 DO BEGIN

; remove comments for these lines if you want color tables for each particle's energy
;;   if colorson eq 1 then begin
;;     mincol=min(particles(ix,*,3))
;;     maxcol=max(particles(ix,*,3))
;;   endif

   
   ystart=reform(fps(ix,*))
   yp=ystart
   CMTODEINT, ystart, eps, -h1, string1, yp, ll, bb
   
   n_data=size(yp)

   ; km: forward field lines
   xf=reform(yp(0,*))
   yf=reform(yp(1,*))
   zf=reform(yp(2,*))
      
   if particleson eq 1 then begin  
     xdata=particles(ix,number,0)
     ydata=particles(ix,number,1)
     zdata=particles(ix,number,2) 
   endif
   

   if colorson eq 1 then begin
     partcolor = (maxcol-particles(ix,number,3))/(maxcol-mincol) * [0,0,255] + $
                 (particles(ix,number,3)-mincol)/(maxcol-mincol) * [255,0,0] + $
                 (particles(ix,number,3)-mincol)*(maxcol-particles(ix,number,3))/$
                 (maxcol-mincol)/(maxcol-mincol) * 4 * [0,255,0] 
   endif else begin
     partcolor=[255,0,0]
   endelse

    IF pno eq 0 THEN BEGIN
     myplot=plot(ll, bb, xr=[0,max(ll)],yr=[0,0.8], xtit='arc length',ytit='field strength',tit='time ='+STRING((ttt-0.00)*10.0,FORMAT='(F4.1)'))
     if particleson eq 1 then begin
      FOR ip=0,nparticles-1 DO BEGIN
       xdata=particles(ip,number,0)
       ydata=particles(ip,number,1)
       zdata=particles(ip,number,2)
       ypp=[xdata,ydata,zdata]/lscl
       CMTODEINT, [xdata,ydata,zdata]/lscl, eps, h1, string1, ypp, llp
       xpoint=max(abs(llp))
       ypoint=particles(ip,number,3)/mu(ip)/0.01
       print, xpoint, ypoint
       myplot=scatterplot([xpoint,xpoint],[ypoint,ypoint],SYMBOL='star',SYM_COLOR='r',OVERPLOT=1)
      endfor
     endif
    ;; text1=TEXT(0.45,0.9,0.9,'time ='+STRING((ttt-1.05)*100.0,FORMAT='(F4.1)'),/NORMAL,FONT_SIZE=14,ONGLASS=1)
    ;; text1=TEXT(0.45,0.9,0.9,'time ='+STRING((ttt-0.00)*10.0,FORMAT='(F5.1)'),/NORMAL,FONT_SIZE=14,ONGLASS=1)
    ENDIF ELSE BEGIN
     myplot=plot(ll, bb, OVERPLOT=1)
     if particleson eq 1 then begin
       ypp=[xdata,ydata,zdata]/lscl
       CMTODEINT, [xdata,ydata,zdata]/lscl, eps, h1, string1, ypp, llp
       xpoint=max(abs(llp))
       ypoint=particles(ix,number,3)/mu(ix)/0.01
       myplot=scatterplot([xpoint,xpoint],[ypoint,ypoint],SYMBOL='star',SYM_COLOR='r',OVERPLOT=1)
     endif
    ENDELSE

   ; XPLOT3D, ts[*,0], ts[*,1], ts[*,2], COLOR=[0,0,0], NAME='start', SYMBOL=oSymbol, THICK=2, /OVERPLOT
   
   ; XPLOT3D, te[*,0], te[*,1], te[*,2], COLOR=[0,0,0], NAME='end', SYMBOL=oSymbol2, THICK=2, /OVERPLOT
   ;  xplot3D, te[*,0]/lscl, te[*,1]/lscl, te[*,2]/zscl, COLOR=[0,0,0], NAME='atest', SYMBOL=oAR3, THICK=2, /OVERPLOT
    pno=pno+1
ENDFOR


;;;myplot.Save, "frames/fsframe"+STRING(number+1,FORMAT='(I4.4)')+".png"


;;;myplot.Close


;;----2D Plot of Field Strength v arclength------------------------------------


endfor

end
