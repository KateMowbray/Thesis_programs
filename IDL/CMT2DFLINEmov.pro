;PRO CMTFLINEmov
; CMTFLINE
;  numerically integrates field lines according to field specified in derivs
;  then displays
;  km: adjusted from CMTFLINE to generate frames for a movie showing evolution
;      in time of multiple field lines
;  REQUIRES - ODEINT, derivs, RKQC

; km: Creates a movie to show evolution of a 2d CMT field with time

@CMTODEINT
common someothername, length, c1, Bzf, ttt

 lscl=1e7
 zscl=0.1
 dims=[1.,1.,1.]
 ;mysymsize=[lscl,lscl,lscl]*0.1*dims
 mysymsize=[1.0,1.0,1.0]*7.*dims  


t0=0.0
t1=100.0
noframes=10;801
;;timearr=t0+FINDGEN(noframes+1)*(t1-t0)/noframes
;;fieldarr=FLTARR(noframes+1)
rotations=0.0
steps=t0+FINDGEN(noframes+1)/(noframes)*(t1-t0)

particleson=0 ; km: Set to 1 to overplot particles on field lines
colorson=0    ; km: Set to 1 to colour particles according to their energies
colorflines=0 ; km: Set to 1 to colour field lines according to field strength
              ;     Not yet optimised but can comment out/in two lines relating to fline color

;;----- define variables ----------------
eps=0.0000001
Bzf=0.00
h1=0.01     	    ; step size?
string1='derivs' ; km: This needs to match with string1 in CMTODEINT

;;----- pick initial conditions ------------
ystart=dblarr(3)
r1=dblarr(3)
r2=dblarr(3)
rsteps=intarr(3)


;r1[0]=-0.9129	;xstart
;r2[0]=0.9129 	;xend
r1[0]=-0.50
r2[0]=-0.001
rsteps[0]=11 	;nx
r1[1]=1.25	;ystart
r2[1]=1.25 	;yend
rsteps[1]=1	;ny
;r1[2]=0.231564  	;zstart
;r2[2]=0.231564   	;zend
r1[2]=-0.001;-1.200
r2[2]=-0.001;0.0001
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

ttt=1.05 ; km: initial normalised of 1.05 necessary for Giuliani et al CMT model, 0.0 for braking jet
;;;ttt=0.00

particles=FLTARR(np,noframes+1,4)

FOR ix=0,rsteps[0]-1 DO BEGIN
 FOR iy=0,rsteps[1]-1 DO BEGIN
  FOR iz=0,rsteps[2]-1 DO BEGIN
   
   ;gds=[gdx,gdy,gdz]
   i=[ix,iy,iz]
   ;yorig= R1+lbox*(i*1.0D0)*gds
   ystart= R1+lbox*(i*1.0D0)*[gdx,gdy,gdz]
   CMTODEINT, ystart, eps, h1, string1, yp
   fps[ix*rsteps[2]*rsteps[1]+iy*rsteps[2]+iz,*]=[INTERPOL(yp(0,*),yp(1,*),0), $
      0.0,INTERPOL(yp(2,*),yp(1,*),0)]
   if particleson eq 1 then begin
     ds=getrdata(ix*rsteps[2]*rsteps[1]+iy*rsteps[2]+iz+1,/ke,/fields,wkdir="../Data/DataR/")
     particles(ix*rsteps[2]*rsteps[1]+iy*rsteps[2]+iz,*,0)=INTERPOL(ds.x,ds.t,steps,/SPLINE)
     particles(ix*rsteps[2]*rsteps[1]+iy*rsteps[2]+iz,*,1)=INTERPOL(ds.y,ds.t,steps,/SPLINE)
     particles(ix*rsteps[2]*rsteps[1]+iy*rsteps[2]+iz,*,2)=INTERPOL(ds.z,ds.t,steps,/SPLINE)
     particles(ix*rsteps[2]*rsteps[1]+iy*rsteps[2]+iz,*,3)=INTERPOL(ds.ke,ds.t,steps,/SPLINE)
;;     fieldarr=SQRT(INTERPOL(ds.B[0,*],ds.t,steps,/SPLINE)^2+INTERPOL(ds.B[1,*],ds.t,steps,/SPLINE)^2+$
;;                   INTERPOL(ds.B[2,*],ds.t,steps,/SPLINE)^2)
   endif
  ENDFOR
 ENDFOR
ENDFOR

if colorson eq 1 then begin
  mincol=min(particles(*,*,3))
  maxcol=max(particles(*,*,3))
endif

for number=0,noframes do begin
ttt=(105.0+t0+(t1-t0)*number/noframes)/100.0 
print, ttt

pno=0

FOR ix=0,np-1 DO BEGIN

; km: remove comments for these lines if you want color tables for each particle's energy
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
     IF number eq 0 THEN BEGIN
       ctable=COLORTABLE(33,INDICES=FINDGEN(256)*max(bb)/256.00)
       maxb=0.8;max(bb)
;       print, maxb
     ENDIF
     myplot=plot(xf, yf, xrange=[-1.5,1.5], yrange=[0,4],$;POSITION=[0.22,0.1,0.95,0.85],$; COLOR=flinecol, $
     ztitle='y (10Mm)', xtitle='x (10Mm)',ytitle='y (10Mm)' $
;     ,RGB_TABLE=ctable, VERT_COLORS=BYTSCL(bb)*max(bb)/maxb $ ;33
     )
;     col=COLORBAR(RANGE=[0,maxb/100],orientation=1, POSITION=[0.12,0.1,0.16,0.85],TITLE='field strength (T)')
     if particleson eq 1 then begin
       myplot=plot([xdata,xdata]/lscl,[ydata,ydata]/lscl,SYM_OBJECT=orb(),SYM_COLOR=partcolor,SYM_SIZE=1.5*0.8,OVERPLOT=1)
     endif
     text1=TEXT(0.45,0.9,'time ='+STRING((ttt-1.05)*100.0,FORMAT='(F5.1)'),/NORMAL,FONT_SIZE=14,ONGLASS=1)
     ;;text1=TEXT(0.45,0.9,'time ='+STRING((ttt-0.00)*10.0,FORMAT='(F5.1)'),/NORMAL,FONT_SIZE=14,ONGLASS=1)
    ENDIF ELSE BEGIN
     myplot=plot(xf, yf, $; COLOR=[0,0,255], OVERPLOT=1)
            ; km: for black field lines instead of coloured ones comment out next line and removed comment midway through previous one
;            RGB_TABLE=ctable,VERT_COLORS=BYTSCL(bb)*max(bb)/maxb, $ ;33
            OVERPLOT=1)
     if particleson eq 1 then begin
       myplot=plot([xdata,xdata]/lscl,[ydata,ydata]/lscl,SYM_OBJECT=orb(),SYM_COLOR=partcolor,SYM_SIZE=1.5*0.8,OVERPLOT=1)
     endif
    ENDELSE
   ; XPLOT3D, ts[*,0], ts[*,1], ts[*,2], COLOR=[0,0,0], NAME='start', SYMBOL=oSymbol, THICK=2, /OVERPLOT
   
   ; XPLOT3D, te[*,0], te[*,1], te[*,2], COLOR=[0,0,0], NAME='end', SYMBOL=oSymbol2, THICK=2, /OVERPLOT
   ;  xplot3D, te[*,0]/lscl, te[*,1]/lscl, te[*,2]/zscl, COLOR=[0,0,0], NAME='atest', SYMBOL=oAR3, THICK=2, /OVERPLOT
    pno=pno+1

ENDFOR


;;myplot=PLOT(timearr,fieldarr[0:number],xrange=[t0,t1],layout=[2,2,4],/CURRENT)

;;myplot=PLOT(timearr,particles[0,[0:number],3],xr=[t0,t1],yr=[min(particles[0,*,3]),max(particles[0,*,3])],layout=[2,2,2],/CURRENT)

myplot.Save, "frames/CMTframe"+STRING(number+1,FORMAT='(I4.4)')+".png"

;;myplot.Close


;;----3D Plot of Field Lines----------------------------------------------


endfor

end
