FUNCTION sech, a

num = 1.0/cosh(a)
return, num
end

FUNCTION sec, a

num = 1.0/cos(a)
return, num
end

FUNCTION derivs, x, y
common someothername, length, c1, Bzf, ttt
;------------------------------------------------------------------
;******************************************************************
; right hand side of field line ODE for 2D CMT field
;******************************************************************
!except=0
;
;JT B-field

pi=3.14159
length=1e7
Lscl=1.0
dydx = dblarr(3)
Lx=Lscl
Lz=Lscl
Ly=Lscl
c1=(-1.5)*1e5

X0=y[0]/Lscl
;Y0=0.4*ttt*alog(1+y[1]/0.4/ttt)*(1+tanh((y[1]-1)*0.9))/2 + $
;   (1-tanh((y[1]-1)*0.9))/2*y[1]
;ab model stuff

tzero=100.0
yprm=1.0
xprm=100.0
d=0.3
w=0.1
k=7.0
chi=1.0
yzero=7.0
zeta=0.3
tprm=5.0
w2=2.3
alpha=1.0
beta=0.5
szero=0.5
vphi=(-2.0)
sigma=4.0
s=0.7

r= y[0] ;SQRT(y[0]^2+y[2]^2)

lilf = (1.0-tanh(ttt-tzero))/2.0
lilg = 1-lilf
F = (1.0-tanh(y[1]-yprm))/2.0
G = (tanh(r+xprm)-tanh(r-xprm))/2.0
J = d*exp((-r*r*y[1])/w)
bigT = k*(1.0+(chi*sin(pi*y[1]/yzero)-1.0)*(1-tanh(zeta*(ttt-tprm)))/2.0)* $
       tan(pi*r*r/2.0/w2)*tanh(y[1])
phi= alpha*(1.0+(chi*sin(pi*y[1]/yzero)-1.0)*(1.0-tanh(zeta*(ttt-tprm)))/2.0)*$
     (y[1]+1.0)^beta*(szero*r*r+1.0)^beta*(y[1]-vphi*sigma*tanh(ttt/sigma)-yzero+J)+bigT

Y0=(s*ALOG(1.0+y[1]/s)*(1-F*G)+y[1]*F*G+y[1]*(1.0+tanh(phi))/2.0)*lilf+y[1]*lilg

dlilfdt = (-sech(ttt-tzero)*sech(ttt-tzero))/2.0
dlilgdt = (-1.0)*dlilfdt
dFdy = (-sech(y[1]-yprm)*sech(y[1]-yprm))/2.0
dGdr = (sech(r-xprm)*sech(r-xprm)-sech(r+xprm)*sech(r+xprm))/2.0
dJdr = (-2*r*y[1])/w*d*exp((-r*r*y[1])/w)
dJdy = (-r*r)/w*d*exp((-r*r*y[1])/w)
dbigTdr = k*(1.0+(chi*sin(pi*y[1]/yzero)-1.0)*(1-tanh(zeta*(ttt-tprm)))/2.0)*$
          r*pi/w2*sec(pi/2.0*r*r/w2)*sec(pi/2.0*r*r/w2)*tanh(y[1])
dbigTdy = k*chi*pi/yzero*cos(pi*y[1]/yzero)*(1-tanh(zeta*(ttt-tprm)))/2.0*$
          tan(pi/2.0*r*r/w2)*tanh(y[1]) + $
          k*(1.0+(chi*sin(pi*y[1]/yzero)-1.0)*(1-tanh(zeta*(ttt-tprm)))/2.0)*$
          tan(pi/2.0*r*r/w2)*sech(y[1])*sech(y[1])
dbigTdt = k*(chi*sin(pi*y[1]/yzero)-1.0)*(-zeta*sech(zeta*(ttt-tprm))*sech(zeta*(ttt-tprm))/2.0)*$
          tan(pi/2.0*r*r/w2)*tanh(y[1])
dphidr = alpha*(1.0+(chi*sin(pi*y[1]/yzero)-1.0)*(1.0-tanh(zeta*(ttt-tprm)))/2.0)*$
         (y[1]+1.0)^beta*beta*2*r*szero*(szero*r*r+1.0)^(beta-1.0)*$
         (y[1]-vphi*sigma*tanh(ttt/sigma)-yzero+J) + $
         alpha*(1.0+(chi*sin(pi*y[1]/yzero)-1.0)*(1.0-tanh(zeta*(ttt-tprm)))/2.0)*$
         (y[1]+1.0)^beta*(szero*r*r+1.0)^beta*dJdr + dbigTdr
dphidy = alpha*(pi/yzero*chi*cos(pi*y[1]/yzero)*(1.0-tanh(zeta*(ttt-tprm)))/2.0)*$
         (y[1]+1.0)^beta*(szero*r*r+1.0)^beta*(y[1]-vphi*sigma*tanh(ttt/sigma)-yzero+J) + $
         alpha*(1.0+(chi*sin(pi*y[1]/yzero)-1.0)*(1.0-tanh(zeta*(ttt-tprm)))/2.0)*$
         beta*(y[1]+1.0)^(beta-1.0)*(szero*r*r+1.0)^beta*(y[1]-vphi*sigma*tanh(ttt/sigma)-yzero+J) + $
         alpha*(1.0+(chi*sin(pi*y[1]/yzero)-1.0)*(1.0-tanh(zeta*(ttt-tprm)))/2.0)*$
         (y[1]+1.0)^beta*(szero*r*r+1.0)^beta*(1.0+dJdy) + dbigTdy 
dphidt = alpha*(chi*sin(pi*y[1]/yzero)-1.0)*(-zeta*sech(zeta*(ttt-tprm))*sech(zeta*(ttt-tprm)))/2.0*$
         (y[1]+1.0)^beta*(szero*r*r+1.0)^beta*(y[1]-vphi*sigma*tanh(ttt/sigma)-yzero+J) + $ 
         alpha*(1.0+(chi*sin(pi*y[1]/yzero)-1.0)*(1.0-tanh(zeta*(ttt-tprm)))/2.0)*$
         (y[1]+1.0)^beta*(szero*r*r+1.0)^beta*(-vphi*sech(ttt/sigma)*sech(ttt/sigma)) + dbigTdt

drdx = 1.0 ;y[0]/SQRT(y[0]^2+y[2]^2)
drdz = 0.0 ;y[2]/SQRT(y[0]^2+y[2]^2)

dY0dr = (s*ALOG(1.0+y[1]/s)*(-F*dGdr)+y[1]*F*dGdr+y[1]*dphidr*sech(phi)*sech(phi)/2.0)*lilf

dY0dX = drdx*dY0dr
dY0dY = (1.0/(1.0+y[1]/s)*(1.0-F*G)+s*ALOG(1.0+y[1]/s)*(-dFdy*G)+F*G+ $
        y[1]*dFdy*G+(1.0+tanh(phi))/2.0+y[1]*dphidy*sech(phi)*sech(phi)/2.0)*lilf+lilg
dY0dZ = drdz*dY0dr

dA0dX0=32*c1*(Y0+1)*X0/((4*X0*X0-4*X0+1+4*Y0*Y0+8*Y0+4)* $ 
                        (4*X0*X0+4*X0+1+4*Y0*Y0+8*Y0+4))
dA0dY0=4*c1*(-4*X0*X0+1+4*Y0*Y0+8*Y0+4)/((4*X0*X0-4*X0+1+4*Y0*Y0+8*Y0+4)* $
                                         (4*X0*X0+4*X0+1+4*Y0*Y0+8*Y0+4))

dX0dX=1.0
dX0dY=0.0
dX0dZ=0.0

;dY0dX = 0.0
;dY0dY = 0.5*(1.0+tanh((y[1]-1.0)*0.9))/(1.0+y[1]/0.4/ttt)+ $
;        0.5*0.4*ttt*alog(1+y[1]/0.4/ttt)*(1.0-tanh((y[1]-1.0)*0.9)^2.0)*0.9 - $
;        0.5*(1.0-tanh((y[1]-1.0)*0.9)^2.0)*0.9*y[1]+0.5-0.5*tanh((y[1]-1.0)*0.9)
;dY0dZ = 0.0


Bx=(dY0dY*dA0dY0 - dX0dY*dA0dX0)/length
By=(-(dY0dX*dA0dY0 + dX0dX*dA0dX0))/length
Bz=(dX0dX*dY0dY-dY0dX*dX0dY)*Bzf
 
 dydx[0]=Bx/0.01
 dydx[1]=By/0.01
 dydx[2]=Bz/0.01

 
return, dydx
end

 derivs2_5, x, y
common someothername, length, c1, Bzf, ttt
;------------------------------------------------------------------
;******************************************************************
; right hand side of field line ODE for 2.5D CMT w. shear
;******************************************************************
!except=0
;
;JT B-field

length=1e7
Lscl=1.0
dydx = dblarr(3)
Lx=Lscl
Lz=Lscl
Ly=Lscl
c1=(-1.5)*1e5

a3d=1.0
delta=1.0

X0=y[0]
Y0=0.4*ttt*alog(1+y[1]/0.4/ttt)*(1+tanh((y[1]-1)*0.9))/2 + $
   (1-tanh((y[1]-1)*0.9))/2*y[1]
Z0=y[2]+delta*(Y0-y[1])*y[0]/(a3d*a3d+y[0]*y[0])

;print, "X0 = ", X0
;print, "Y0 = ", Y0

dA0dX0=32*c1*(Y0+1)*X0/((4*X0*X0-4*X0+1+4*Y0*Y0+8*Y0+4)* $ 
                        (4*X0*X0+4*X0+1+4*Y0*Y0+8*Y0+4))
dA0dY0=4*c1*(-4*X0*X0+1+4*Y0*Y0+8*Y0+4)/((4*X0*X0-4*X0+1+4*Y0*Y0+8*Y0+4)* $
                                         (4*X0*X0+4*X0+1+4*Y0*Y0+8*Y0+4))

;print, "dA0dX0 = ", dA0dX0  
;print, "dA0dY0 = ", dA0dY0

dX0dX=1.0
dX0dY=0.0
dX0dZ=0.0

dY0dX = 0.0
dY0dY = 0.5*(1.0+tanh((y[1]/Lscl-1.0)*0.9))/(1.0+y[1]/Lscl/0.4/ttt)+ $
        0.5*0.4*ttt*alog(1+y[1]/Lscl/0.4/ttt)*(1.0-tanh((y[1]/Lscl-1.0)*0.9)^2.0)*0.9 - $
        0.5*(1.0-tanh((y[1]/Lscl-1.0)*0.9)^2.0)*0.9*y[1]/Lscl+0.5-0.5*tanh((y[1]/Lscl-1.0)*0.9)
dY0dZ = 0.0

dZ0dX = 0.0;delta*(Y0-y[1])*(a3d*a3d-y[0]*y[0])/ $ 
        ;((a3d*a3d+y[0]*y[0])*(a3d*a3d+y[0]*y[0]))
dZ0dY = 0.0;delta*(dY0dY-1.0)*y[0]/(a3d*a3d+y[0]*y[0])
dZ0dZ = 1.0

Bx=(dY0dY*dA0dY0 - dX0dY*dA0dX0)/length
By=(-(dY0dX*dA0dY0 + dX0dX*dA0dX0))/length
Bz=((dY0dX*dZ0dY-dZ0dX*dY0dY)*(dA0dY0) - $
    (dZ0dX*dX0dY-dX0dX*dZ0dY)*(dA0dX0))/length + $
   (dX0dX*dY0dY-dY0dX*dX0dY)*Bzf
 
;print, "Bx = ", Bx
;print, "By = ", By
;print, "Bz = ", Bz

 dydx[0]=Bx/0.01
 dydx[1]=By/0.01
 dydx[2]=Bz/0.01

 
return, dydx
end


FUNCTION derivs3, x, y
common someothername, length, c1, Bzf, ttt
;------------------------------------------------------------------
;******************************************************************
; right hand side of field line ODE for 3D CMT
;******************************************************************
!except=0
;
;JT B-field

length=1e7
Lscl=1.0
dydx = dblarr(3)
Lx=Lscl
Lz=Lscl
Ly=Lscl
c1=(-1.5)*1e5
cc=0.4;(0.4/SQRT(1.05))^(2.0/3.0);SQRT(0.4/1.05)
esp=1.0 ;2.0

a3d=1.0
delta=1.0
ay=1.0
yh=0.0

Y0=(cc*ttt)^esp*alog(1+y[1]/(cc*ttt)^esp)*(1+tanh((y[1]-1)*0.9))/2 + $
   (1-tanh((y[1]-1)*0.9))/2*y[1]
;X0=y[0]-delta*(Y0-y[1])*y[2]/(a3d*a3d+y[0]*y[0]+y[2]*y[2])
;Z0=y[2]+delta*(Y0-y[1])*y[0]/(a3d*a3d+y[0]*y[0]+y[2]*y[2])
X0=y[0]-delta*(Y0-y[1])*y[2]/(a3d*a3d+y[0]*y[0]+ay*(y[1]-yh)*ay*(y[1]-yh)+y[2]*y[2])
Z0=y[2]+delta*(Y0-y[1])*y[0]/(a3d*a3d+y[0]*y[0]+ay*(y[1]-yh)*ay*(y[1]-yh)+y[2]*y[2])

;print, "X0 = ", X0
;print, "Y0 = ", Y0
;print, "Z0 = ", Z0

a1=SQRT((X0+0.5)*(X0+0.5)+(Y0+1)*(Y0+1)+Z0*Z0)
a2=SQRT((X0-0.5)*(X0-0.5)+(Y0+1)*(Y0+1)+Z0*Z0)

B0x = c1*((X0+0.5)/a1/a1/a1 - (X0-0.5)/a2/a2/a2)
B0y = c1*((Y0+1)/a1/a1/a1 - (Y0+1)/a2/a2/a2)
B0z = c1*(Z0/a1/a1/a1 - Z0/a2/a2/a2)

;print, "B0x = ", B0x
;print, "B0y = ", B0y
;print, "B0z = ", B0z

dY0dX = 0.0
;;dY0dY = 0.5*(1.0+tanh((y[1]/Lscl-1.0)*0.9))/(1.0+y[1]/Lscl/(cc*ttt)^esp)+ $
;;        0.5*(cc*ttt)^esp*alog(1+y[1]/Lscl/(cc*ttt)^esp)*(1.0-tanh((y[1]/Lscl-1.0)*0.9)^2.0)*0.9 - $
;;        0.5*(1.0-tanh((y[1]/Lscl-1.0)*0.9)^2.0)*0.9*y[1]/Lscl+0.5-0.5*tanh((y[1]/Lscl-1.0)*0.9)

dY0dY = 0.5*(1.0+tanh((y[1]-1.0)*0.9))/(1.0+y[1]/(cc*ttt)^esp)+ $
        0.5*(cc*ttt)^esp*alog(1+y[1]/(cc*ttt)^esp)*(1.0-tanh((y[1]-1.0)*0.9)^2.0)*0.9 - $
        0.5*(1.0-tanh((y[1]-1.0)*0.9)^2.0)*0.9*y[1]+0.5-0.5*tanh((y[1]-1.0)*0.9)

dY0dZ = 0.0

;dX0dX=1.0+2.0*y[0]*delta*(Y0-y[1])*y[2]/((a3d*a3d+y[0]*y[0]+y[2]*y[2])* $
;                                        (a3d*a3d+y[0]*y[0]+y[2]*y[2]))
;dX0dY=delta*(1.0-dY0dY)*y[2]/(a3d*a3d+y[0]*y[0]+y[2]*y[2])
;dX0dZ=(-1.0)*delta*(Y0-y[1])*(a3d*a3d+y[0]*y[0]-y[2]*y[2])/ $
;      ((a3d*a3d+y[0]*y[0]+y[2]*y[2])*(a3d*a3d+y[0]*y[0]+y[2]*y[2]))

dX0dX=1.0+2.0*y[0]*delta*(Y0-y[1])*y[2]/$
                ((a3d*a3d+y[0]*y[0]+ay*(y[1]-yh)*ay*(y[1]-yh)+y[2]*y[2])*$
                 (a3d*a3d+y[0]*y[0]+ay*(y[1]-yh)*ay*(y[1]-yh)+y[2]*y[2]))+$
      delta*(-dY0dX)*y[2]/(a3d*a3d+y[0]*y[0]+ay*(y[1]-yh)*ay*(y[1]-yh)+y[2]*y[2])
dX0dY=delta*(1.0-dY0dY)*y[2]/(a3d*a3d+y[0]*y[0]+ay*(y[1]-yh)*ay*(y[1]-yh)+y[2]*y[2]) + $
      2.0*delta*(Y0-y[1])*ay*ay*(y[1]-yh)*y[2]/$
      ((a3d*a3d+y[0]*y[0]+ay*(y[1]-yh)*ay*(y[1]-yh)+y[2]*y[2])*$
       (a3d*a3d+y[0]*y[0]+ay*(y[1]-yh)*ay*(y[1]-yh)+y[2]*y[2]))
dX0dZ=(-1.0)*delta*(Y0-y[1])*(a3d*a3d+y[0]*y[0]+ay*(y[1]-yh)*ay*(y[1]-yh)-y[2]*y[2])/ $
      ((a3d*a3d+y[0]*y[0]+ay*(y[1]-yh)*ay*(y[1]-yh)+y[2]*y[2])*$
       (a3d*a3d+y[0]*y[0]+ay*(y[1]-yh)*ay*(y[1]-yh)+y[2]*y[2]))+$
      delta*(-dY0dZ)*y[2]/(a3d*a3d+y[0]*y[0]+ay*(y[1]-yh)*ay*(y[1]-yh)+y[2]*y[2])


;dZ0dX = delta*(Y0-y[1])*(a3d*a3d-y[0]*y[0]+y[2]*y[2])/ $ 
;        ((a3d*a3d+y[0]*y[0]+y[2]*y[2])*(a3d*a3d+y[0]*y[0]+y[2]*y[2]))
;dZ0dY = delta*(dY0dY-1.0)*y[0]/(a3d*a3d+y[0]*y[0]+y[2]*y[2])
;dZ0dZ = 1.0-2.0*delta*(Y0-y[1])*y[0]*y[2]/((a3d*a3d+y[0]*y[0]+y[2]*y[2])* $
;                                          (a3d*a3d+y[0]*y[0]+y[2]*y[2]))

dZ0dX = delta*(Y0-y[1])*(a3d*a3d-y[0]*y[0]+ay*(y[1]-yh)*ay*(y[1]-yh)+y[2]*y[2])/ $ 
        ((a3d*a3d+y[0]*y[0]+ay*(y[1]-yh)*ay*(y[1]-yh)+y[2]*y[2])*$
         (a3d*a3d+y[0]*y[0]+ay*(y[1]-yh)*ay*(y[1]-yh)+y[2]*y[2]))+$
        delta*dY0dX*y[0]/(a3d*a3d+y[0]*y[0]+ay*(y[1]-yh)*ay*(y[1]-yh)+y[2]*y[2])  
dZ0dY = delta*(dY0dY-1.0)*y[0]/(a3d*a3d+y[0]*y[0]+ay*(y[1]-yh)*ay*(y[1]-yh)+y[2]*y[2]) -$
        2.0*delta*(Y0-y[1])*ay*ay*(y[1]-yh)*y[0]/$
            ((a3d*a3d+y[0]*y[0]+ay*(y[1]-yh)*ay*(y[1]-yh)+y[2]*y[2])*$
             (a3d*a3d+y[0]*y[0]+ay*(y[1]-yh)*ay*(y[1]-yh)+y[2]*y[2]))
dZ0dZ = 1.0-2.0*delta*(Y0-y[1])*y[0]*y[2]/$
            ((a3d*a3d+y[0]*y[0]+ay*(y[1]-yh)*ay*(y[1]-yh)+y[2]*y[2])* $
             (a3d*a3d+y[0]*y[0]+ay*(y[1]-yh)*ay*(y[1]-yh)+y[2]*y[2]))+$
        delta*dY0dZ*y[0]/(a3d*a3d+y[0]*y[0]+ay*(y[1]-yh)*ay*(y[1]-yh)+y[2]*y[2])

Bx=((dY0dY*dZ0dZ-dZ0dY*dY0dZ)*B0x+(dZ0dY*dX0dZ-dX0dY*dZ0dZ)*B0y+  $
    (dX0dY*dY0dZ-dY0dY*dX0dZ)*B0z)/length
By=((dY0dZ*dZ0dX-dZ0dZ*dY0dX)*B0x+(dZ0dZ*dX0dX-dX0dZ*dZ0dX)*B0y+  $
    (dX0dZ*dY0dX-dY0dZ*dX0dX)*B0z)/length
Bz=((dY0dX*dZ0dY-dZ0dX*dY0dY)*B0x+(dZ0dX*dX0dY-dX0dX*dZ0dY)*B0y+  $
    (dX0dX*dY0dY-dY0dX*dX0dY)*B0z)/length ; + (dX0dX*dY0dY-dY0dX*dX0dY)*Bzf

;print, "Bx = ", Bx
;print, "By = ", By
;print, "Bz = ", Bz

;print,"************************"
 
 dydx[0]=Bx/0.01
 dydx[1]=By/0.01
 dydx[2]=Bz/0.01

 
return, dydx
end

FUNCTION derivsab, x, y
common someothername, length, c1, Bzf, ttt
;------------------------------------------------------------------
;******************************************************************
; right hand side of field line ODE for 3D CMT
;******************************************************************
!except=0
;
;JT B-field

length=1e7
Lscl=1.0
dydx = dblarr(3)
Lx=Lscl
Lz=Lscl
Ly=Lscl
c1=(-1.5)*1e5
cc=0.4;(0.4/SQRT(1.05))^(2.0/3.0);SQRT(0.4/1.05)
esp=1.0 ;2.0
pi=3.14159

a3d=1.0
delta=0.00
ay=1.0/2.0
yh=0.0

tzero=100.0;3.5;100.0
yprm=1.0
xprm=100.0
d=0.3
w=0.1
k=7.0
chi=1.0
yzero=7.0
zeta=0.3
tprm=5.0
w2=2.3
alpha=1.0
beta=0.5
szero=0.5
vphi=(-2.0)
sigma=4.0
s=0.7

r= SQRT(y[0]^2+y[2]^2)

lilf = (1.0-tanh(ttt-tzero))/2.0

lilg = 1-lilf
F = (1.0-tanh(y[1]-yprm))/2.0
G = (tanh(r+xprm)-tanh(r-xprm))/2.0
J = d*exp((-r*r*y[1])/w)
bigT = k*(1.0+(chi*sin(pi*y[1]/yzero)-1.0)*(1-tanh(zeta*(ttt-tprm)))/2.0)* $
       tan(pi*y[0]*y[0]/2.0/w2)*tanh(y[1])
phi= alpha*(1.0+(chi*sin(pi*y[1]/yzero)-1.0)*(1.0-tanh(zeta*(ttt-tprm)))/2.0)*$
     (y[1]+1.0)^beta*(szero*r*r+1.0)^beta*(y[1]-vphi*sigma*tanh(ttt/sigma)-yzero+J)+bigT

Y0=(s*ALOG(1.0+y[1]/s)*(1-F*G)+y[1]*F*G+y[1]*(1.0+tanh(phi))/2.0)*lilf + $;y[1]*lilg
   lilg*((cc*tttt)^esp*alog(1+y[1]/(cc*tttt)^esp)*(1+tanh((y[1]-1)*0.9))/2 + $
    (1-tanh((y[1]-1)*0.9))/2*y[1])

X0=y[0]-delta*(Y0-y[1])*y[2]/(a3d*a3d+y[0]*y[0]+ay*(y[1]-yh)*ay*(y[1]-yh)+y[2]*y[2])
Z0=y[2]+delta*(Y0-y[1])*y[0]/(a3d*a3d+y[0]*y[0]+ay*(y[1]-yh)*ay*(y[1]-yh)+y[2]*y[2])

;print, "X0 = ", X0
;print, "Y0 = ", Y0
;print, "Z0 = ", Z0

a1=SQRT((X0+0.5)*(X0+0.5)+(Y0+1)*(Y0+1)+Z0*Z0)
a2=SQRT((X0-0.5)*(X0-0.5)+(Y0+1)*(Y0+1)+Z0*Z0)

B0x = c1*((X0+0.5)/a1/a1/a1 - (X0-0.5)/a2/a2/a2)
B0y = c1*((Y0+1)/a1/a1/a1 - (Y0+1)/a2/a2/a2)
B0z = c1*(Z0/a1/a1/a1 - Z0/a2/a2/a2)

;print, "B0x = ", B0x
;print, "B0y = ", B0y
;print, "B0z = ", B0z

dlilfdt = (-sech(ttt-tzero)*sech(ttt-tzero))/2.0


dlilgdt = (-1.0)*dlilfdt
dFdy = (-sech(y[1]-yprm)*sech(y[1]-yprm))/2.0
dGdr = (sech(r+xprm)*sech(r+xprm)-sech(r-xprm)*sech(r-xprm))/2.0
dJdr = (-2*r*y[1])/w*d*exp((-r*r*y[1])/w)
dJdy = (-r*r)/w*d*exp((-r*r*y[1])/w)
dbigTdx = k*(1.0+(chi*sin(pi*y[1]/yzero)-1.0)*(1-tanh(zeta*(ttt-tprm)))/2.0)*$
          y[0]*pi/w2*sec(pi/2.0*y[0]*y[0]/w2)*sec(pi/2.0*y[0]*y[0]/w2)*tanh(y[1])
dbigTdr = 0.0
dbigTdy = k*chi*pi/yzero*cos(pi*y[1]/yzero)*(1-tanh(zeta*(ttt-tprm)))/2.0*$
          tan(pi/2.0*y[0]*y[0]/w2)*tanh(y[1]) + $
          k*(1.0+(chi*sin(pi*y[1]/yzero)-1.0)*(1-tanh(zeta*(ttt-tprm)))/2.0)*$
          tan(pi/2.0*y[0]*y[0]/w2)*sech(y[1])*sech(y[1])
dbigTdt = k*(chi*sin(pi*y[1]/yzero)-1.0)*(-zeta*sech(zeta*(ttt-tprm))*sech(zeta*(ttt-tprm))/2.0)*$
          tan(pi/2.0*y[0]*y[0]/w2)*tanh(y[1])
dphidr = alpha*(1.0+(chi*sin(pi*y[1]/yzero)-1.0)*(1.0-tanh(zeta*(ttt-tprm)))/2.0)*$
         (y[1]+1.0)^beta*beta*2*r*szero*(szero*r*r+1.0)^(beta-1.0)*$
         (y[1]-vphi*sigma*tanh(ttt/sigma)-yzero+J) + $
         alpha*(1.0+(chi*sin(pi*y[1]/yzero)-1.0)*(1.0-tanh(zeta*(ttt-tprm)))/2.0)*$
         (y[1]+1.0)^beta*(szero*r*r+1.0)^beta*dJdr + dbigTdr
dphidx = dbigTdx
dphidy = alpha*(pi/yzero*chi*cos(pi*y[1]/yzero)*(1.0-tanh(zeta*(ttt-tprm)))/2.0)*$
         (y[1]+1.0)^beta*(szero*r*r+1.0)^beta*(y[1]-vphi*sigma*tanh(ttt/sigma)-yzero+J) + $
         alpha*(1.0+(chi*sin(pi*y[1]/yzero)-1.0)*(1.0-tanh(zeta*(ttt-tprm)))/2.0)*$
         beta*(y[1]+1.0)^(beta-1.0)*(szero*r*r+1.0)^beta*(y[1]-vphi*sigma*tanh(ttt/sigma)-yzero+J) + $
         alpha*(1.0+(chi*sin(pi*y[1]/yzero)-1.0)*(1.0-tanh(zeta*(ttt-tprm)))/2.0)*$
         (y[1]+1.0)^beta*(szero*r*r+1.0)^beta*(1.0+dJdy) + dbigTdy 
dphidt = alpha*(chi*sin(pi*y[1]/yzero)-1.0)*(-zeta*sech(zeta*(ttt-tprm))*sech(zeta*(ttt-tprm)))/2.0*$
         (y[1]+1.0)^beta*(szero*r*r+1.0)^beta*(y[1]-vphi*sigma*tanh(ttt/sigma)-yzero+J) + $ 
         alpha*(1.0+(chi*sin(pi*y[1]/yzero)-1.0)*(1.0-tanh(zeta*(ttt-tprm)))/2.0)*$
         (y[1]+1.0)^beta*(szero*r*r+1.0)^beta*(-vphi*sech(ttt/sigma)*sech(ttt/sigma)) + dbigTdt

drdx = y[0]/SQRT(y[0]^2+y[2]^2)
drdz = y[2]/SQRT(y[0]^2+y[2]^2)

dY0dr = (s*ALOG(1.0+y[1]/s)*(-F*dGdr)+y[1]*F*dGdr+y[1]*dphidr*sech(phi)*sech(phi)/2.0)*lilf
partdY0dx = y[1]*dphidx*sech(phi)*sech(phi)/2.0*lilf

dY0dX = drdx*dY0dr + partdY0dx
dY0dY = (1.0/(1.0+y[1]/s)*(1.0-F*G)+s*ALOG(1.0+y[1]/s)*(-dFdy*G)+F*G+ $
        y[1]*dFdy*G+(1.0+tanh(phi))/2.0+y[1]*dphidy*sech(phi)*sech(phi)/2.0)*lilf+lilg

dY0dZ = drdz*dY0dr

;dX0dX=1.0+2.0*y[0]*delta*(Y0-y[1])*y[2]/((a3d*a3d+y[0]*y[0]+y[2]*y[2])* $
;                                        (a3d*a3d+y[0]*y[0]+y[2]*y[2]))
;dX0dY=delta*(1.0-dY0dY)*y[2]/(a3d*a3d+y[0]*y[0]+y[2]*y[2])
;dX0dZ=(-1.0)*delta*(Y0-y[1])*(a3d*a3d+y[0]*y[0]-y[2]*y[2])/ $
;      ((a3d*a3d+y[0]*y[0]+y[2]*y[2])*(a3d*a3d+y[0]*y[0]+y[2]*y[2]))

dX0dX=1.0+2.0*y[0]*delta*(Y0-y[1])*y[2]/$
                ((a3d*a3d+y[0]*y[0]+ay*(y[1]-yh)*ay*(y[1]-yh)+y[2]*y[2])*$
                 (a3d*a3d+y[0]*y[0]+ay*(y[1]-yh)*ay*(y[1]-yh)+y[2]*y[2]))+$
      delta*(-dY0dX)*y[2]/(a3d*a3d+y[0]*y[0]+ay*(y[1]-yh)*ay*(y[1]-yh)+y[2]*y[2])
dX0dY=delta*(1.0-dY0dY)*y[2]/(a3d*a3d+y[0]*y[0]+ay*(y[1]-yh)*ay*(y[1]-yh)+y[2]*y[2]) + $
      2.0*delta*(Y0-y[1])*ay*ay*(y[1]-yh)*y[2]/$
      ((a3d*a3d+y[0]*y[0]+ay*(y[1]-yh)*ay*(y[1]-yh)+y[2]*y[2])*$
       (a3d*a3d+y[0]*y[0]+ay*(y[1]-yh)*ay*(y[1]-yh)+y[2]*y[2]))
dX0dZ=(-1.0)*delta*(Y0-y[1])*(a3d*a3d+y[0]*y[0]+ay*(y[1]-yh)*ay*(y[1]-yh)-y[2]*y[2])/ $
      ((a3d*a3d+y[0]*y[0]+ay*(y[1]-yh)*ay*(y[1]-yh)+y[2]*y[2])*$
       (a3d*a3d+y[0]*y[0]+ay*(y[1]-yh)*ay*(y[1]-yh)+y[2]*y[2]))+$
      delta*(-dY0dZ)*y[2]/(a3d*a3d+y[0]*y[0]+ay*(y[1]-yh)*ay*(y[1]-yh)+y[2]*y[2])


;dZ0dX = delta*(Y0-y[1])*(a3d*a3d-y[0]*y[0]+y[2]*y[2])/ $ 
;        ((a3d*a3d+y[0]*y[0]+y[2]*y[2])*(a3d*a3d+y[0]*y[0]+y[2]*y[2]))
;dZ0dY = delta*(dY0dY-1.0)*y[0]/(a3d*a3d+y[0]*y[0]+y[2]*y[2])
;dZ0dZ = 1.0-2.0*delta*(Y0-y[1])*y[0]*y[2]/((a3d*a3d+y[0]*y[0]+y[2]*y[2])* $
;                                          (a3d*a3d+y[0]*y[0]+y[2]*y[2]))

dZ0dX = delta*(Y0-y[1])*(a3d*a3d-y[0]*y[0]+ay*(y[1]-yh)*ay*(y[1]-yh)+y[2]*y[2])/ $ 
        ((a3d*a3d+y[0]*y[0]+ay*(y[1]-yh)*ay*(y[1]-yh)+y[2]*y[2])*$
         (a3d*a3d+y[0]*y[0]+ay*(y[1]-yh)*ay*(y[1]-yh)+y[2]*y[2]))+$
        delta*dY0dX*y[0]/(a3d*a3d+y[0]*y[0]+ay*(y[1]-yh)*ay*(y[1]-yh)+y[2]*y[2])  
dZ0dY = delta*(dY0dY-1.0)*y[0]/(a3d*a3d+y[0]*y[0]+ay*(y[1]-yh)*ay*(y[1]-yh)+y[2]*y[2]) -$
        2.0*delta*(Y0-y[1])*ay*ay*(y[1]-yh)*y[0]/$
            ((a3d*a3d+y[0]*y[0]+ay*(y[1]-yh)*ay*(y[1]-yh)+y[2]*y[2])*$
             (a3d*a3d+y[0]*y[0]+ay*(y[1]-yh)*ay*(y[1]-yh)+y[2]*y[2]))
dZ0dZ = 1.0-2.0*delta*(Y0-y[1])*y[0]*y[2]/$
            ((a3d*a3d+y[0]*y[0]+ay*(y[1]-yh)*ay*(y[1]-yh)+y[2]*y[2])* $
             (a3d*a3d+y[0]*y[0]+ay*(y[1]-yh)*ay*(y[1]-yh)+y[2]*y[2]))+$
        delta*dY0dZ*y[0]/(a3d*a3d+y[0]*y[0]+ay*(y[1]-yh)*ay*(y[1]-yh)+y[2]*y[2])

Bx=((dY0dY*dZ0dZ-dZ0dY*dY0dZ)*B0x+(dZ0dY*dX0dZ-dX0dY*dZ0dZ)*B0y+  $
    (dX0dY*dY0dZ-dY0dY*dX0dZ)*B0z)/length
By=((dY0dZ*dZ0dX-dZ0dZ*dY0dX)*B0x+(dZ0dZ*dX0dX-dX0dZ*dZ0dX)*B0y+  $
    (dX0dZ*dY0dX-dY0dZ*dX0dX)*B0z)/length
Bz=((dY0dX*dZ0dY-dZ0dX*dY0dY)*B0x+(dZ0dX*dX0dY-dX0dX*dZ0dY)*B0y+  $
    (dX0dX*dY0dY-dY0dX*dX0dY)*B0z)/length ; + (dX0dX*dY0dY-dY0dX*dX0dY)*Bzf

;print, "Bx = ", Bx
;print, "By = ", By
;print, "Bz = ", Bz

;print,"************************"
 
 dydx[0]=Bx/0.01
 dydx[1]=By/0.01
 dydx[2]=Bz/0.01

 
return, dydx
end


;##################################################################
; Runge-Kutta routines with adaptive step size control taken
; from Numerical Recipes
; rk4 from idl lib
;##################################################################
PRO RKQC, y, dydx,x,htry,eps,yscal,hdid,hnext,string1

common someothername, length, c1, Bzf, ttt
;----------------------------------------------------
; *****************************************************************
;  RKQC is the stepper routine which makes one adaptive step 
;  forward
; *****************************************************************
;-----------------------------------------------------------------
; assign some fixed values to variables
;-----------------------------------------------------------------

pgrow = -0.2
pshrnk= -0.25
fcor  = 1.0/15.0
safety = 0.9
errcon = 6.0e-4
ytemp = dblarr(3)
ysav  = ytemp
dysav = ysav

;-----------------------------------------------------------------
; save initial values for step since it is needed twice
;-----------------------------------------------------------------

xsav=x
ysav=y
dysav=dydx

;-----------------------------------------------------------------
; try stepsize htry first and define hh as stepsize for the 
; two half-steps
;-----------------------------------------------------------------

h = htry
jump1: hh= 0.5 *h


;-----------------------------------------------------------------
; take two half-steps first
;-----------------------------------------------------------------

ytemp=rk4(ysav,dysav,xsav,hh,string1)

x=xsav+hh

dydx=derivs(x,ytemp)
y=rk4(ytemp,dydx,x,hh,string1)

;-----------------------------------------------------------------
; one step next
;-----------------------------------------------------------------

x=xsav+h
ytemp=rk4(ysav,dysav,xsav,h,string1)

;-----------------------------------------------------------------
; now evaluate the accuracy
;-----------------------------------------------------------------

errmax = 0.
ytemp = y-ytemp

errmax= max(abs(ytemp/yscal))

;-----------------------------------------------------------------
; scale error to required accuracy eps and reduce stepsize if 
; necessary
;-----------------------------------------------------------------

errmax = errmax/eps
if errmax gt 1.0 then begin
   h = safety*double(h)*(errmax^pshrnk)
   goto, jump1
endif else begin            
;-----------------------------------------------------------------
; if stepsize ok, compute stepsize for next step
;-----------------------------------------------------------------

   hdid = h                 
   if errmax gt errcon then begin
       hnext = safety*h*(errmax^pgrow)
   endif else begin
       hnext = 4.0*h
   endelse
endelse



;-----------------------------------------------------------------
; reduce truncation error to fifth order
;-----------------------------------------------------------------

y = y + ytemp*fcor

hnext=hnext

end ;_____of RKQC
;*****************************************************************
; ODEINT is the driver routine for the ODE solver
;*****************************************************************

PRO CMTODEINT, ystart,eps,h1,string1, yp,ll,bb
!except=0
common someothername, length, c1, Bzf, ttt
;common path, kount;, yp
;maxsteps = 1e4-2
maxsteps= 500000L
tiny = 1.0e-20

;-----------------------------------------------------------------
; we always start at arclength(=x) zero; set initial steplength to
; h1; and step counter to zero
;-----------------------------------------------------------------
x=0
h=h1
kount = 0L
;-----------------------------------------------------------------
; assign initial values and store them in solution vector
;-----------------------------------------------------------------
y = ystart
yp=dblarr(3)
yp(*,0) = ystart(*)

ll=dblarr(1) ;these lines for arc length / field strength plot
ll(0)=0
bb=dblarr(1)

;-----------------------------------------------------------------
; start integrating ODE
;-----------------------------------------------------------------
for iit=1,maxsteps do begin
	dydx = derivs(x,y)
;-----------------------------------------------------------------
; for monitoring accuracy; see Numerical Recipes for explanation
;-----------------------------------------------------------------
  	yscal = abs(y) + abs(h*dydx) + tiny
;-----------------------------------------------------------------
; make one adaptive step forward and write result to yp
;-----------------------------------------------------------------
  RKQC,y, dydx,x,h,eps,yscal,hdid,hnext,string1
  kount = kount + 1
  yp=[[yp],[dblarr(3)]]
  yp(*,kount) = y(*)

  ll=[[ll],[dblarr(1)]] ;these lines for arc length / field strength plot
  ll(kount)=ll(kount-1)+SQRT((y(0)-yp(0,kount-1))^2+(y(1)-yp(1,kount-1))^2+$
                             (y(2)-yp(2,kount-1))^2)
  ;;ll[kount]=y(0)
  bb=[[bb],[dblarr(1)]]
  bb(kount-1)=sqrt(dydx(0)^2+dydx(1)^2+dydx(2)^2)

  IF ((ABS(y(0)) GT 10.0) OR (y(1) LT 0.0)) THEN BEGIN
  ;IF ((ABS(y(0)) GT 10.0) OR (ABS(y(2)) GT 10.0) OR (y(1) GT 20.0) OR (y(1) LT 0.0)) THEN BEGIN
  ;IF ((ABS(y(0)) GT 10.0) OR (ABS(y(1)) GT 10.0) OR (y(2) GT 20.0) OR (y(2) LT 0.0)) THEN BEGIN
  ;IF ((ABS(y(0)) GT 10.0) OR (ABS(y(2)) GT 10.0) OR (y(1) GT 20.0) OR (y(1) LT 0.0)) THEN BEGIN
  ;IF ((ABS(y(0)) GT 10.0) OR (y(2) LT 0.0) OR (y(2) GT 20.0) OR (y(1) GT 20.0) OR (y(1) LT 0.0)) THEN BEGIN

  conearr=ASIN(SQRT(bb/max(bb)))
  conearr=conearr*180/3.14159

;;;  plot, ll, bb, xtit='arc length', ytit='field strength'
;  print, min(bb(0:kount-1))
;  print,ll(7300)
;  print,bb(7300)
;  print,yp(*,7300)
  

   return
  ENDIF
  ;print, y(*)
  ;STOP

  ; limiting step size growth  
  IF abs(hnext) lt 0.1 THEN h = hnext
  
  
endfor

end;____________of ODEINT
