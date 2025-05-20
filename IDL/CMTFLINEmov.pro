;PRO CMTFLINEmov
; CMTFLINE
;  numerically integrates field lines according to field specified in derivs
;  then displays
;  km: adjusted from CMTFLINE to generate frames for a movie showing evolution
;      in time of multiple field lines in 3D
;  REQUIRES - ODEINT, derivs, RKQC

@CMTODEINT
common someothername, length, c1, Bzf, ttt

 lscl=1e7
 zscl=0.1
 dims=[1.,1.,1.]
 ;mysymsize=[lscl,lscl,lscl]*0.1*dims
 mysymsize=[1.0,1.0,1.0]*7.*dims  


t0=0.0
t1=100.0;50.0
noframes=1;51
;;timearr=t0+FINDGEN(noframes+1)*(t1-t0)/noframes
;;fieldarr=FLTARR(noframes+1)
rotations=0.0
steps=t0+FINDGEN(noframes+1)/(noframes)*(t1-t0)

particleson=0 ; Set to 1 to overplot particles on field lines
colorson=0    ; Set to 1 to colour particles according to their energies
colorflines=0 ; Set to 1 to colour field lines according to field strength
              ; Not yet optimised but can comment out/in two lines relating to fline color

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
r1[0]=-1.200
r2[0]=-0.001
rsteps[0]=3 	;nx
r1[1]=1.25	;ystart
r2[1]=1.25 	;yend
rsteps[1]=1	;ny
;r1[2]=0.231564  	;zstart
;r2[2]=0.231564   	;zend
r1[2]=0.001
r2[2]=1.20
rsteps[2]=3 	;nz


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

particles=FLTARR(np,noframes+1,4)

FOR ix=0,rsteps[0]-1 DO BEGIN
 FOR iy=0,rsteps[1]-1 DO BEGIN
  FOR iz=0,rsteps[2]-1 DO BEGIN
   print, ix, iz
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
;ttt=(105.0+t0+(t1-t0)*number/noframes)/100.0 ;simple model
ttt=(0.0+t0+(t1-t0)*number/noframes)/10.0     ;braking jet
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
     myplot=plot3d(xf, zf, yf, xrange=[-2,2], zrange=[0,6.5], yrange=[-0.2,2.0],$; COLOR=flinecol, $
     ztitle='y (10Mm)', xtitle='x (10Mm)',ytitle='z (10Mm)' $
;     ,RGB_TABLE=33, VERT_COLORS=BYTSCL(bb) $
     )
     if particleson eq 1 then begin
       myplot=plot3d([xdata,xdata]/lscl,[zdata,zdata]/lscl,[ydata,ydata]/lscl,SYM_OBJECT=orb(),SYM_COLOR=partcolor,SYM_SIZE=1.5*0.8,OVERPLOT=1)
     endif
     FOR i1 = 0, 13 do begin
       xline = FINDGEN(2)*(-1.3) + 0.05
       yline = FLTARR(2)+1.25
       zline = FLTARR(2)-0.05 + i1*0.1
       myplot=plot3d(xline,zline,yline,color='r',overplot=1)
       myplot=plot3d((-1)*zline,(-1)*xline,yline,color='r',overplot=1)
     endfor
    ;; text1=TEXT(0.45,0.9,0.9,'time ='+STRING((ttt-1.05)*100.0,FORMAT='(F5.1)'),/NORMAL,FONT_SIZE=14,ONGLASS=1)
     text1=TEXT(0.45,0.9,0.9,'time ='+STRING((ttt-0.00)*10.0,FORMAT='(F5.1)'),/NORMAL,FONT_SIZE=14,ONGLASS=1)
    ENDIF ELSE BEGIN
     myplot=plot3d(xf, zf, yf, $; COLOR=[0,0,255], OVERPLOT=1)
;            RGB_TABLE=33,VERT_COLORS=BYTSCL(bb), $
            OVERPLOT=1)
     if particleson eq 1 then begin
       myplot=plot3d([xdata,xdata]/lscl,[zdata,zdata]/lscl,[ydata,ydata]/lscl,SYM_OBJECT=orb(),SYM_COLOR=partcolor,SYM_SIZE=1.5*0.8,OVERPLOT=1)
     endif
    ENDELSE
   ; XPLOT3D, ts[*,0], ts[*,1], ts[*,2], COLOR=[0,0,0], NAME='start', SYMBOL=oSymbol, THICK=2, /OVERPLOT
   
   ; XPLOT3D, te[*,0], te[*,1], te[*,2], COLOR=[0,0,0], NAME='end', SYMBOL=oSymbol2, THICK=2, /OVERPLOT
   ;  xplot3D, te[*,0]/lscl, te[*,1]/lscl, te[*,2]/zscl, COLOR=[0,0,0], NAME='atest', SYMBOL=oAR3, THICK=2, /OVERPLOT
    pno=pno+1

ENDFOR


myplot.Rotate, 360.0*rotations/noframes*number, /ZAXIS

myplot.Rotate, 60.0, /ZAXIS
;myplot.Rotate, -15.0, /YAXIS
myplot.Rotate, 30.0, /ZAXIS

;;myplot=PLOT(timearr,fieldarr[0:number],xrange=[t0,t1],layout=[2,2,4],/CURRENT)

;;myplot=PLOT(timearr,particles[0,[0:number],3],xr=[t0,t1],yr=[min(particles[0,*,3]),max(particles[0,*,3])],layout=[2,2,2],/CURRENT)

myplot.Save, "frames/CMTframe"+STRING(number+1,FORMAT='(I4.4)')+".png"

;;myplot.Rotate, -360.0*rotations/noframes*number, /ZAXIS

;;myplot.Rotate, -90.0, /ZAXIS

;;myplot.Save, "frames/CMTsideframe"+STRING(number+1,FORMAT='(I4.4)')+".png"

;;myplot.Rotate, 90.0, /ZAXIS

;These lines give a top down view
;;myplot.Rotate, -30.0, /ZAXIS
;;myplot.Rotate, 60.0, /XAXIS
;;myplot.Save, "frames/topframe"+STRING(number+1,FORMAT='(I4.4)')+".png"


;These next 2 lines give an x,y plane projection, useful for 2D model (but erases x axis!!)
;;myplot.Rotate, -30.0, /ZAXIS
;;myplot.Rotate, -30.0, /XAXIS
;;myplot.Save, "frames/CMTsideframe"+STRING(number+1,FORMAT='(I4.4)')+".png"

;;myplot.Close


;;----3D Plot of Field Lines----------------------------------------------


endfor

end
