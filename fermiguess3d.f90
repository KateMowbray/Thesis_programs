PROGRAM fermiguessfor

IMPLICIT NONE

INTEGER, PARAMETER :: dp=kind(1.d0)
REAL(dp) :: lscl, bscl, tscl, oneuponl, oneupontl
REAL(dp) :: c1, d, esp, cc, a, a2, Lv, a3d, delta, ay, yh
REAL(dp) :: xmin, xmax, ymin, ymax, zmin, zmax, Emin, Emax, pitchmin, pitchmax
INTEGER :: xsteps, ysteps, zsteps, Esteps, pitchsteps
INTEGER :: xindex, yindex, zindex, Eindex, pitchindex
REAL(dp) :: xinit ,yinit ,zinit, Einit, pitchinit, tinit
REAL(dp) :: looptoptol
INTEGER :: steps, looptopsteps, Bydir, ndata, datapoint
REAL(dp), DIMENSION(3) :: Rinit, fp, BBinit
REAL(dp) :: mu, Binit, eguess
REAL(dp), DIMENSION(:), ALLOCATABLE :: dataarray

 lscl=1e7
 bscl=0.01
 tscl=100.0

 oneuponl=1.0/lscl
 oneupontl=1.0/lscl/tscl

 c1 = 0.15e6
 d=lscl
 esp=1.0
 cc=0.4
 a=0.9
 a2=0.9
 Lv=lscl
 a3d = 1.0
 delta=0.0
 ay=1.0
 yh=0.0

 steps=60!51
 looptopsteps=81!81!41+0*10
 looptoptol=0.2-0*0.05

 xmin = -0.500e7!-0.4916833e7
 xmax = -0.001e7
 xsteps=1001

 ymin = 1.25e7
 ymax = 1.25e7
 ysteps=1 

 zmin = -0.001e7
 zmax = -0.001e7
 zsteps=1

 Emin = 5.5e3
 Emax = 5.5e3
 Esteps=1

 pitchmin = 70.0
 pitchmax = 70.0
 pitchsteps=1

 ndata=xsteps*ysteps*zsteps*Esteps*pitchsteps

 ALLOCATE(dataarray(ndata))
 datapoint=1

 Tinit = 105.00

Tinit = Tinit/tscl

DO xindex = 0, xsteps-1
  xinit = xmin + xindex*(xmax-xmin)/(xsteps-1)
  DO yindex = 0, ysteps-1
    yinit = ymin + yindex*(ymax-ymin)/(ysteps)
    DO zindex = 0, zsteps-1
      zinit = zmin + zindex*(zmax-zmin)/(zsteps)
      Rinit(1) = xinit/lscl
      Rinit(2) = yinit/lscl
      Rinit(3) = zinit/lscl
!      PRINT*, 'R init = ', Rinit
      BBinit = twoDODEint(Rinit,Tinit)
      IF (BBinit(2) .gt. 0.0) then
        Bydir = 1
      else
        Bydir = (-1)
      end if
      Binit = magnitude3(BBinit)*bscl
      fp=footpoint(Rinit, Tinit)
      DO Eindex = 0, Esteps-1
        Einit = Emin + Eindex*(Emax-Emin)/Esteps
!!        PRINT*, 'E init = ', Einit
        DO pitchindex = 0, pitchsteps-1
          pitchinit = pitchmin + pitchindex*(pitchmax-pitchmin)/pitchsteps
!!          PRINT*, 'pitch init = ', pitchinit

          mu = Einit*SIN(pitchinit*3.14159/180.0)*SIN(pitchinit*3.14159/180.0)/Binit

          eguess=heightweighted(fp, pitchinit, Binit, mu, steps, looptopsteps, looptoptol, Rinit(1))
!          PRINT*, 'eguess = ', eguess
          dataarray(datapoint) = eguess
          datapoint = datapoint+1
        END DO
      END DO
    END DO
  END DO
END DO  

OPEN(1,FORM='UNFORMATTED',STATUS='REPLACE',FILE='fermiguess.dat')
WRITE(1) dataarray
 CLOSE(1)

DEALLOCATE(dataarray)

CONTAINS

SUBROUTINE INTERPOL(xData, yData, xVal, yVal)

IMPLICIT NONE

REAL(dp), INTENT(IN) :: xData(:), yData(:), xVal(:)
REAL(dp), INTENT(OUT) :: yVal(:)
INTEGER :: inputIndex, dataIndex
REAL(dp) :: minXdata, maxXdata, minYdata, weight

minXdata=xData(1)
maxXdata=xData(size(xData))
dataIndex=0


yval(1)=yData(1)

DO inputindex=2, size(xVal)-1
 DO WHILE (xData(dataIndex+1) .LT. xVal(inputIndex))
   dataIndex = dataIndex+1
 END DO
 weight = (xVal(inputindex) - xData(dataIndex))/(xData(dataindex+1)-xData(dataIndex))
 yval(inputindex) = (1.0-weight)*yData(dataindex) + weight*yData(dataIndex+1)
END DO

yval(size(xval))=yData(size(yData)-1)

END SUBROUTINE INTERPOL

FUNCTION DOT(aa, bb)
IMPLICIT NONE

REAL(dp), DIMENSION(3) :: aa, bb
REAL(dp) :: DOT

  ! 3D dot product, useful for later vector calculations
  DOT = aa(1)*bb(1) + aa(2)*bb(2) + aa(3)*bb(3) 

END FUNCTION DOT

FUNCTION CROSS(aa, bb)
IMPLICIT NONE

REAL(dp), DIMENSION(3) :: aa, bb, cross
  ! similar to above, 3D cross product, useful for later vector calculations like E
  cross(1)= aa(2)*bb(3) - aa(3)*bb(2) 
  cross(2)= aa(3)*bb(1) - aa(1)*bb(3) 
  cross(3)= aa(1)*bb(2) - aa(2)*bb(1) 

END FUNCTION CROSS

FUNCTION MAGNITUDE3(arr)
IMPLICIT NONE

REAL(dp), DIMENSION(3) :: arr
REAL(dp) :: magnitude3
  ! gives magnitude of 3D vector, useful for stuff like modB
  magnitude3 = SQRT(arr(1)*arr(1)+arr(2)*arr(2)+arr(3)*arr(3))

END FUNCTION MAGNITUDE3


FUNCTION twoDODEint(R, T)

 ! Returns 2D field for given position and time, useful for field line calculations

IMPLICIT NONE

REAL(dp), DIMENSION(3) :: R
REAL(dp) :: T
REAL(dp), DIMENSION(3) :: dydx, twoDODEint, dXdx, dXdy, dXdz, B0
REAL(dp) :: X0, Y0, Z0, dX0dX, dX0dY, dX0dZ, dY0dX, dY0dY, dY0dZ, dZ0dX, dZ0dY, dZ0dZ
REAL(dp) :: dX0dt, dY0dt, dZ0dt
REAL(dp) :: d1, d2, oneuponden
REAL(dp) :: Bx, By, Bz


oneuponden=1.0/(a3d*a3d+R(1)*R(1)+ay*(R(2)-yh)*ay*(R(2)-yh)+R(3)*R(3))

 Y0=(cc*T)**esp*log (1.+(R(2)/((cc*T)**esp)))*  &	
      (1.+tanh((R(2)-Lv/Lscl)*a))*0.5 +           &
       (1.-tanh((R(2)-Lv/Lscl)*a2))*0.5*R(2)

 X0=R(1)-delta*(Y0-R(2))*R(3)*oneuponden

 Z0=R(3)+delta*(Y0-R(2))*R(1)*oneuponden


 dY0dX=0.
 dY0dY=0.5*(1.+tanh((R(2)-Lv/Lscl)*a))/(1.+R(2)/(cc*T)**esp)+        &
      0.5*(cc*T)**esp*log(1.+R(2)/(cc*T)**esp)*               &
      (1.-tanh((R(2)-Lv/Lscl)*a)**2)*a-0.5*                         &
      (1.-tanh((R(2)-Lv/Lscl)*a2)**2)*a2*R(2)+0.5-0.5*tanh((R(2)-Lv/Lscl)*a2) 
 dY0dZ=0.

 
 dX0dX=1.0+2.0*R(1)*delta*(Y0-R(2))*R(3)*oneuponden*oneuponden - &
       delta*dY0dX*R(3)*oneuponden
 dX0dY=delta*(1.0-dY0dY)*R(3)*oneuponden+ &
       2.0*delta*(Y0-R(2))*ay*ay*(R(2)-yh)*R(3)*oneuponden*oneuponden
 dX0dZ=(-1.0)*delta*(Y0-R(2))*(a3d*a3d+R(1)*R(1)+ay*(R(2)-yh)*ay*(R(2)-yh)-R(3)*R(3))* &
                             oneuponden*oneuponden - &
       delta*dY0dZ*R(3)*oneuponden


 dZ0dX=delta*(Y0-R(2))*(a3d*a3d+ay*(R(2)-yh)*ay*(R(2)-yh)+R(3)*R(3)-R(1)*R(1))* &
                             oneuponden*oneuponden + &
       delta*dY0dX*R(1)*oneuponden
 dZ0dY=delta*(dY0dY-1.0)*R(1)*oneuponden-2.0*delta*(Y0-R(2))*ay*ay*(R(2)-yh)*R(1)*oneuponden*oneuponden
 dZ0dZ=1.0-2*delta*(Y0-R(2))*R(1)*R(3)*oneuponden*oneuponden + &
       delta*dY0dZ*R(1)*oneuponden
   
 dY0dt=0.5*(cc*T)**esp*esp*log(1.+R(2)/(cc*T)**esp)*          &
      (1.+tanh((R(2)-Lv/Lscl)*a))/T-0.5*R(2)*esp*                  &
      (1.+tanh((R(2)-Lv/Lscl)*a))/(T*(1.+R(2)/(cc*T)**esp)) 
 dX0dt=(-1.0)*delta*dY0dt*R(3)*oneuponden
 dZ0dt=delta*dY0dt*R(1)*oneuponden         


 d1 = SQRT((X0+0.5)*(X0+0.5)+(Y0+1)*(Y0+1)+Z0*Z0)
 d2 = SQRT((X0-0.5)*(X0-0.5)+(Y0+1)*(Y0+1)+Z0*Z0)

 B0(1) = (X0+0.5)/d1/d1/d1 - (X0-0.5)/d2/d2/d2
 B0(2) = (Y0+1)/d1/d1/d1 - (Y0+1)/d2/d2/d2
 B0(3) = Z0/d1/d1/d1 - Z0/d2/d2/d2

 B0 = B0 * c1

 dXdx(1)=dX0dX
 dXdx(2)=dY0dX
 dXdx(3)=dZ0dX

 dXdy(1)=dX0dY
 dXdy(2)=dY0dY
 dXdy(3)=dZ0dY

 dXdz(1)=dX0dZ
 dXdz(2)=dY0dZ
 dXdz(3)=dZ0dZ

 dydx(1)= DOT(CROSS(dXdy,dXdz),B0)*oneuponL/bscl
 dydx(2)= DOT(CROSS(dXdz,dXdx),B0)*oneuponL/bscl
 dydx(3)= DOT(CROSS(dXdx,dXdy),B0)*oneuponL/bscl

twoDODEint=dydx

END FUNCTION twoDODEint

FUNCTION footpoint(R, T)

 ! Calculate the footpoint position for a field line that passes through the point xx at time t.
 ! h is negative and small, this calculation is only made once, but the footpoint position is used
 ! many times, so a higher degree of accuracy is preferable. When tracing field lines from the 
 ! foot point, h should switch sign.

 IMPLICIT NONE

 REAL(dp), DIMENSION(3) :: R
 REAL(dp) :: T
 REAL(dp) :: h1, hrk
 REAL(dp), DIMENSION(3) :: k1,k2,k3,k4
 REAL(dp), DIMENSION(3) :: fp, footpoint
 REAL(dp), DIMENSION(3) :: rrprev, dydx, rr, rrrk, rrprevrk

 h1=0.002!0.0005 ! set stepsize


 ! Ensure that we integrate down the loop towards the foot point
 h1 = h1 * (-Bydir)
 hrk = 10*h1

 rrrk=R
 rrprevrk=R

  DO WHILE (rrrk(2) .GT. 0.0)
   k1 = twoDODEint(rrrk,T)
   k2 = twoDODEint(rrrk+k1*hrk/2.0,T)
   k3 = twoDODEint(rrrk+k2*hrk/2.0,T)
   k4 = twoDODEint(rrrk+k3*hrk,T)
   rrprevrk=rrrk
   rrrk = rrrk + hrk*(k1+2*k2+2*k3+k4)/6.0
  end do

 rr=rrprevrk
 rrprev=rrprevrk 

 DO WHILE (rr(2) .GT. 0.0)
   dydx = twoDODEint(rr,T)
   rrprev = rr
   rr = rr + h1*dydx
 END DO

 ! linear interpolation of footpoint position (y=0)

 fp(1) = rrprev(1)+(-rrprev(2))*(rr(1)-rrprev(1))/(rr(2)-rrprev(2)) 
 fp(2) = rrprev(2)+(-rrprev(2))*(rr(2)-rrprev(2))/(rr(2)-rrprev(2)) 
 fp(3) = rrprev(3)+(-rrprev(2))*(rr(3)-rrprev(3))/(rr(2)-rrprev(2)) 

 footpoint = fp

END FUNCTION footpoint

FUNCTION looptoppts(R, T, xval, nosteps)

 ! Finds equally spaced set of points on the looptop

 IMPLICIT NONE

 REAL(dp), DIMENSION(3) :: R
 REAL(dp) :: T, xval
 INTEGER :: nosteps
 REAL(dp), DIMENSION(nosteps,3) :: pointarr, looptoppts
 REAL(dp), DIMENSION(nosteps) :: sarray, xarray, yarray, zarray
 REAL(dp) :: h1, s, hrk
 REAL(dp), DIMENSION(3) :: rrprev, dydx, lptp, rr, rrrk, rrprevrk
 REAL(dp), DIMENSION(3) :: k1,k2,k3,k4
 INTEGER :: k, i, q1
 REAL(dp), DIMENSION(:,:), ALLOCATABLE :: yp

 h1=min(0.005,abs(xval)*0.1) ! Choose stepsize

 !ensure we integrate up the field line
 h1 = h1 * Bydir

 rr=R
 rrprev=R
 s=0.0
 k=1

 ! We integrate over the looptop once so we know how to allocate the matrix of positions

 if (abs(xval) .lt. 0.001) then
   DO q1 = 1,nosteps
    pointarr(q1,1)=R(1)
    pointarr(q1,2)=R(2)
    pointarr(q1,3)=R(3)
   END DO

 else


 hrk = 10*h1

 rrrk=R
 rrprevrk=R

  DO WHILE ((abs(rrrk(1))-xval)*(abs(rrprevrk(1))-xval) .GT. 0.0)
   k1 = twoDODEint(rrrk,T)
   k2 = twoDODEint(rrrk+k1*hrk/2.0,T)
   k3 = twoDODEint(rrrk+k2*hrk/2.0,T)
   k4 = twoDODEint(rrrk+k3*hrk,T)
   rrprevrk=rrrk
   rrrk = rrrk + hrk*(k1+2*k2+2*k3+k4)/6.0
  end do


 rr=rrprevrk
 rrprev=rrprevrk 

 DO WHILE (rr(1)*rrprev(1) .GE. 0.0)
   dydx = twoDODEint(rr,T)
   rrprev=rr
   rr = rr + h1*dydx
   if (ABS(rr(1)) .lt. xval) then 
     k = k+1
   end if
 END DO

 ALLOCATE(yp(4,k))
 yp=0.0
 
! We can use the above routine to find when we're nearly at the loop top and use this to speed up
! the next routine (we don't have to integrate up to the loop top twice)
 rr=rrprevrk
 rrprev=rrprevrk
 k=1

 DO WHILE (rr(1)*rrprev(1) .GE. 0.0) 
   dydx = twoDODEint(rr,T)
   rrprev=rr
   rr = rr + h1*dydx
   if (ABS(rr(1)) .lt. xval) then 
   yp(1,k)=rr(1)
   yp(2,k)=rr(2)
   yp(3,k)=rr(3)
   yp(4,k)=s
   k=k+1
   s=s+MAGNITUDE3(dydx)*h1
   end if
 END DO


 DO i = 1, nosteps
   sarray(i) = (i-1.0)/(nosteps-1.0)*s
 END DO


 CALL INTERPOL(yp(4,:),yp(1,:),sarray,xarray) 
 CALL INTERPOL(yp(4,:),yp(2,:),sarray,yarray)
 CALL INTERPOL(yp(4,:),yp(3,:),sarray,zarray)


 DO q1 = 1,nosteps
  pointarr(q1,1)=xarray(q1)
  pointarr(q1,2)=yarray(q1)
  pointarr(q1,3)=zarray(q1)
 END DO

 
 yp=0.0
DEALLOCATE(yp)

end if

 looptoppts = pointarr


END FUNCTION looptoppts

FUNCTION timeinloops(fp, mu, T, looptoptol, btol)

 ! Function to estimate the proportion of time in an orbit that is spent in the Fermi acceleration
 ! region. Often updated to reflect changes in field topology and particle energies.

IMPLICIT NONE

REAL(dp), DIMENSION(3) :: fp
REAL(dp) :: mu, T, looptoptol, btol
REAL(dp) :: timeinloops
REAL(dp), DIMENSION(3) :: rr, rrprev, dydx
REAL(dp) :: h1, ustep, mirdist, lptpdist
REAL(dp), DIMENSION(3) :: k1,k2,k3,k4

  rr=fp
  rrprev=fp
  mirdist=0.0
  lptpdist=0.0
  h1=0.01!0.05!0.001
  ustep=0.0

 !ensure we integrate up the field line
 h1 = h1 * Bydir

  DO WHILE (rr(1)*rrprev(1) .GE. 0.0) 
  ! Make sure we haven't crossed the loop top yet
   k1 = twoDODEint(rr,T)
   k2 = twoDODEint(rr+k1*h1/2.0,T)
   k3 = twoDODEint(rr+k2*h1/2.0,T)
   k4 = twoDODEint(rr+k3*h1,T)
   rrprev=rr
   rr = rr + h1*(k1+2*k2+2*k3+k4)/6.0
   if (magnitude3(k1)*bscl .lt. Btol) then 
   ! Check if we're within the mirror points
     ustep = SQRT(mu*btol-mu*magnitude3(k1)*bscl)
     mirdist = mirdist+magnitude3(k1)/ustep
      if (abs(rr(1)) .lt. looptoptol) then 
     ! Check if we're in the Fermi acceleration region
       lptpdist=lptpdist+magnitude3(k1)/ustep
      end if
    end if
  END DO
 
 timeinloops = 1.0*lptpdist/mirdist
 ! Take the ratio of time spent in acceleration region to time spent between mirror points

 if (timeinloops .gt. 1.0 .OR. mirdist .EQ. 0) then 
   timeinloops = 1.0
 ! Set the time ratio to one if the mirror points are within the acceleration region. Above
 ! method would calculate as proportion of time spent in acceleration region as being over 100%
 end if

!PRINT*, 'timeinloops = ', timeinloops

END FUNCTION timeinloops

FUNCTION t0calc(fp, mu, T, initialB)

 ! Function to find t0, an estimate of when particles first enter the acceleration region

IMPLICIT NONE

REAL(dp), DIMENSION(3) :: fp
REAL(dp) :: mu, T, initialB
REAL(dp) :: t0, t0calc
REAL(dp), DIMENSION(3) :: rr, rrprev, dydx
REAL(dp) :: h1, ustep, btol
REAL(dp), DIMENSION(3) :: k1,k2,k3,k4

  rr=fp
  rrprev=fp
  Btol=Einit/mu
  h1=0.02!0.001
  ustep=0.0
  t0=0.0

 !ensure we integrate up the field line
 h1 = h1 * Bydir

  DO WHILE (rr(1)*rrprev(1) .GE. 0.0)
   k1 = twoDODEint(rr,T)
   k2 = twoDODEint(rr+k1*h1/2.0,T)
   k3 = twoDODEint(rr+k2*h1/2.0,T)
   k4 = twoDODEint(rr+k3*h1,T)
   rrprev=rr
   rr = rr + h1*(k1+2*k2+2*k3+k4)/6.0
   if (magnitude3(k1)*bscl .lt. Btol) then 
     ustep=SQRT(2*(mu*(btol-MAGNITUDE3(k1)*bscl))*1.6e-19/9.1e-31)/lscl*tscl
     if (magnitude3(k1)*bscl .lt. initialB) then
        t0 = t0 + h1*magnitude3(k1)/ustep
     else 
        t0 = t0 + 2*h1*magnitude3(k1)/ustep
     end if
   ! this needs fine tuning but the idea here is that we double the time between the initial position
   ! and mirror point since a particle initially moving towards the footpoint will cross this portion
   ! of the field line twice. If it initially moves towards the loop top it won't cross this portion
   end if
  END DO

 t0calc = t0

!PRINT*, 't0 = ', t0calc

END FUNCTION t0calc

SUBROUTINE ExBdrift(R, tstart, tsteps, trange, ffs)

 ! Returns the position and field strength after the ExB drift has transported a point for a 
 ! specified time

IMPLICIT NONE

REAL(dp), DIMENSION(3), INTENT(INOUT) :: R
REAL(dp), INTENT(IN) :: tstart, trange
INTEGER, INTENT(IN) :: tsteps
REAL(dp), INTENT(OUT) :: ffs
REAL(dp), DIMENSION(3) :: b, E, vf, ue
REAL(dp) :: tstep, T, X0, Y0, Z0, omodb, dX0dX, dX0dY, dX0dZ, dY0dX, dY0dY, dY0dZ, dZ0dX, dZ0dY, dZ0dZ
REAL(dp) :: dX0dt, dY0dt, dZ0dt, determinant
REAL(dp) :: dA0dX0, dA0dY0
INTEGER :: i

 tstep = trange/tsteps

DO i = 0, tsteps-1 +1

 T=tstart+tstep*i

 X0 = R(1)
 Y0 = (cc*T)**esp*log (1.+(R(2)/((cc*T)**esp)))*  &	
      (1.+tanh((R(2)-Lv/Lscl)*a))*0.5 +           &
       (1.-tanh((R(2)-Lv/Lscl)*a2))*0.5*R(2)
 Z0 = R(3)

 dY0dX=0.
 dY0dY=0.5*(1.+tanh((R(2)-Lv/Lscl)*a))/(1.+R(2)/(cc*T)**esp)+        &
      0.5*(cc*T)**esp*log(1.+R(2)/(cc*T)**esp)*               &
      (1.-tanh((R(2)-Lv/Lscl)*a)**2)*a-0.5*                         &
      (1.-tanh((R(2)-Lv/Lscl)*a2)**2)*a2*R(2)+0.5-0.5*tanh((R(2)-Lv/Lscl)*a2) 
 dY0dZ=0.

 dX0dX=1.0
 dX0dY=0.0
 dX0dZ=0.0

 dZ0dX=0.0
 dZ0dY=0.0
 dZ0dZ=1.0

 dY0dt=0.5*(cc*T)**esp*esp*log(1.+R(2)/(cc*T)**esp)*          &
      (1.+tanh((R(2)-Lv/Lscl)*a))/T-0.5*R(2)*esp*                  &
      (1.+tanh((R(2)-Lv/Lscl)*a))/(T*(1.+R(2)/(cc*T)**esp)) 
 dX0dt=0.0
 dZ0dt=0.0 


dA0dX0 = 32.*c1*(Y0*Lscl+d)*Lscl*Lscl*Lscl*X0/  &
 ((4.*Lscl*Lscl*X0*X0-4.*Lscl*Lscl*X0+Lscl*Lscl+4.*Y0*Y0*Lscl*Lscl+8.*Y0*Lscl*d+4.*d*d)* &
 (4.*Lscl*Lscl*X0*X0+4.*Lscl*Lscl*X0+Lscl*Lscl+4.*Y0*Y0*Lscl*Lscl+8.*Y0*Lscl*d+4.*d*d))

dA0dY0 = 4.*c1*Lscl*Lscl*  &
 (-4.*Lscl*Lscl*X0*X0+Lscl*Lscl+4.*Y0*Y0*Lscl*Lscl+8.*Y0*Lscl*d+4.*d*d)/  &
 ((4.*Lscl*Lscl*X0*X0-4.*Lscl*Lscl*X0+Lscl*Lscl+4*Y0*Y0*Lscl*Lscl+8.*Y0*Lscl*d+4.*d*d)*  &
 (4.*Lscl*Lscl*X0*X0+4.*Lscl*Lscl*X0+Lscl*Lscl+4.*Y0*Y0*Lscl*Lscl+8.*Y0*Lscl*d+4.*d*d))



 B(1) = (dY0dY*dA0dY0 - dX0dY*dA0dX0)/lscl/bscl
 B(2) = (-(dY0dX*dA0dY0 + dX0dX*dA0dX0))/lscl/bscl
 B(3) = 0.0


 omodb = 1.0/SQRT(B(1)*B(1)+B(2)*B(2)+B(3)*B(3))


 DETERMINANT = dX0dx*(dY0dy*dZ0dz-dZ0dy*dY0dz)-dX0dy*(dY0dx*dZ0dz-dZ0dx*dY0dz)+&
               dX0dz*(dY0dx*dZ0dy-dZ0dx*dY0dy)
 
 Vf(1)= (-1.0)*(dX0dt*(dY0dy*dZ0dz-dZ0dy*dY0dz) +   &
                 dY0dt*(dX0dz*dZ0dy-dX0dy*dZ0dz) +   &
                 dZ0dt*(dX0dy*dY0dz-dX0dz*dY0dy) )/DETERMINANT
 
 Vf(2)= (-1.0)*(dX0dt*(dZ0dx*dY0dz-dZ0dz*dY0dx) +   &
                 dY0dt*(dX0dx*dZ0dz-dX0dz*dZ0dx) +   &
                 dZ0dt*(dX0dz*dY0dx-dX0dx*dY0dz) )/DETERMINANT

 Vf(3)= (-1.0)*(dX0dt*(dY0dx*dZ0dy-dY0dy*dZ0dx) +   &
                 dY0dt*(dX0dy*dZ0dx-dX0dx*dZ0dy) +   &
                 dZ0dt*(dX0dx*dY0dy-dX0dy*dY0dx) )/DETERMINANT

 E = -CROSS(vf,B)
 
 UE = CROSS(E,B)*oMODB*oMODB

 R=R+ue*tstep

END DO

ffs=MAGNITUDE3(B)*bscl

END SUBROUTINE ExBdrift

SUBROUTINE betatronstep(fp, T, dt, mu, btol, h, betstep, mirpointx)

! Estimates betatron energisation for a particle. No orbit data used.

IMPLICIT NONE

REAL(dp), DIMENSION(3), INTENT(IN) :: fp
REAL(dp), INTENT(IN) :: T, dt, mu, btol, h
REAL(dp), INTENT(OUT) :: betstep, mirpointx
REAL(dp), DIMENSION(3) :: rr, rrprev, dydx, k1, k2, k3, k4
REAL(dp) :: ustep, ffs, ttotal, hrk, orbitlength, h1
REAL(dp), DIMENSION(:,:), ALLOCATABLE :: rrarray
REAL(dp), DIMENSION(:), ALLOCATABLE :: tarray, barray
REAL(dp), DIMENSION(3) :: rrrk, rrprevrk, mirrorpoint, halfstep, halfstepfield
INTEGER :: npts, arraydim, kount, i

  betstep=0.0

  ! Ensure we integrate up the field line
  h1 = h
  h1 = h1 * Bydir

  ustep=0.0
  hrk=10*h1
  rrrk=fp
  rrprevrk=fp
  orbitlength=0.0
  kount=1

! We'll find the length of the orbit first so that we can set up the matrix correctly later on
  DO WHILE (rrrk(1)*rrprevrk(1) .GE. 0.0)
   k1 = twoDODEint(rrrk,T)
   k2 = twoDODEint(rrrk+k1*hrk/2.0,T)
   k3 = twoDODEint(rrrk+k2*hrk/2.0,T)
   k4 = twoDODEint(rrrk+k3*hrk,T)
   rrprevrk=rrrk
   rrrk = rrrk + hrk*(k1+2*k2+2*k3+k4)/6.0
   if (magnitude3(k1)*bscl .lt. Btol) then 
     orbitlength = orbitlength + MAGNITUDE3(rrrk-rrprevrk)
   else
     mirrorpoint = rrrk
   end if
  end do

  rr = mirrorpoint
  rrprev=mirrorpoint
  mirpointx = abs(mirrorpoint(1))
  arraydim = ceiling(abs(1.5*orbitlength/h1)) + 1 ! Allow some leeway for the dimensions of the array

  ALLOCATE(rrarray(3,arraydim))
  ALLOCATE(barray(arraydim))
  ALLOCATE(tarray(arraydim))

  rrarray=0.0
  barray=0.0
  tarray=0.0

  DO WHILE (rr(1)*rrprev(1) .GE. 0.0 .and. orbitlength .gt. 0.0)
   dydx = twoDODEint(rr,T)
   rrprev=rr
   if (magnitude3(dydx)*bscl .le. Btol) then 
     ustep = SQRT(mu*btol-mu*magnitude3(dydx)*bscl)
 !    print*, ustep, mu*btol-mu*magnitude3(dydx)*bscl, magnitude3(dydx)*bscl, btol
!     if (ustep .ne. 0.0) then
!!       if (kount .eq. 1) then
!!         halfstep = rr + h1/2.0*dydx/magnitude3(dydx)
!!         halfstepfield = twoDODEint(halfstep,T)
!!         ustep = (SQRT(mu*btol-mu*magnitude3(dydx)*bscl)+ &
!!                  SQRT(mu*btol-mu*magnitude3(halfstepfield)*bscl))/2.0
!!       end if
       rrarray(:,kount)=rr
       tarray(kount)=magnitude3(dydx)/ustep
       barray(kount)=magnitude3(dydx)*bscl
       kount=kount+1
!     end if
   end if
   rr = rr + h1*dydx/magnitude3(dydx)
   if (rr(2) .lt. rrprev(2) .and. rr(2) .gt. 0.0) then
   print*, rr
!   print*, dydx
   end if

  end do

!!!  tarray(kount) = tarray(kount-1)

 ! Evolve position in time and find ratio of field strength increase

 ! Deal with case where orbit is narrow enough to not be picked up by the routines above
 IF (kount .EQ. 1) then
   tarray(1)    = 1.0
   rrarray(:,1) = rr
   dydx = twoDODEint(rr,T)
   barray(1) = magnitude3(dydx)*bscl
   CALL driftevo(rr,T,1,dt,ffs)
   betstep = betstep+mu/Einit*(ffs - barray(1))
 end if

 npts = kount
 ttotal = SUM(tarray)


 DO i = 0, npts-2 
  CALL driftevo(rrarray(:,i+1),T,2,dt,ffs)
  betstep = betstep+mu/Einit*(ffs - barray(i+1))*tarray(i+1)/ttotal
 ! print*, betstep, ffs, barray(i+1), tarray(i+1), ttotal
 ! Find the change in field strengths across the points on the orbit, average to find betatron
 ! contribution. We weight by time spent at each point.
 END DO

  rrarray=0.0
  barray=0.0
  tarray=0.0

  DEALLOCATE(rrarray)
  DEALLOCATE(barray)
  DEALLOCATE(tarray)

!PRINT*, 'betstep = ', betstep

END SUBROUTINE betatronstep 



FUNCTION heightweighted(fp, pitch, initfs, mu, steps, looptopsteps, looptoptol, xinit)

 ! This has gone through a few iterations, hence the out of date function name. This was previously
 ! used as a more sophisticated estimate than height swept out by loop top as a guess for Fermi
 ! acceleration, height was weighted by curvature for a better guess. Since E x B drift can be
 ! directly calculated from field data alone, there is no longer a need to keep track of the height
 ! swept out (effectively E x B drift integrated over time). Current iteration of this code only
 ! needs initial pitch angle, initial field strength at particle position and the foot point position
 ! that lies on the same field line as the initial position. No particle orbit calculations are 
 ! required for any of this data. Code not optimised atm, assumes 2D with particular parameter values

IMPLICIT NONE

REAL(dp), DIMENSION(3) :: fp
REAL(dp) :: pitch, initfs, mu, looptoptol, xinit
INTEGER :: steps, looptopsteps
REAL(dp) :: weightedh, heightweighted
REAL(dp), DIMENSION(3) :: lptp, R
REAL(dp) :: t0, btol, T, timeweight, betstep, bettotal, U, omodb, gradbt, determinant, curvature
REAL(dp) :: uaverage, curveterm, fenterm, fpfs, beth, gam, mirpointx
INTEGER :: i, j, zcount
REAL(dp), DIMENSION(:), ALLOCATABLE :: carray, uu
REAL(dp), DIMENSION(:,:), ALLOCATABLE :: ptsarray, lptppts
REAL(dp) :: X0, Y0, Z0, y, d1, d2, den, oneuponden
REAL(dp), DIMENSION(3) :: B, E, vf, UE, gradb, dbdx, dbdy, dbdz, dbdt, B0, dXdx, dXdy, dXdz
REAL(dp), DIMENSION(3) :: d2Xdx2, d2Xdxdy, d2Xdxdz, d2Xdxdt, d2Xdydx, d2Xdy2, d2Xdydz, d2Xdydt
REAL(dp), DIMENSION(3) :: d2Xdzdx, d2Xdzdy, d2Xdz2, d2Xdzdt, d2Xdtdx, d2Xdtdy, d2Xdtdz, d2Xdt2
REAL(dp), DIMENSION(3) :: dB0dx, dB0dy, dB0dz, dB0dt
REAL(dp) :: dX0dX, dX0dY,dX0dZ, dY0dX, dY0dY,dY0dZ,dZ0dX,dZ0dY,dZ0dZ, dX0dt, dY0dt,dZ0dt
REAL(dp), DIMENSION(3) :: dB0dX0, dB0dY0, dB0dZ0
REAL(dp) :: d2X0dXdY,d2X0dYdX,d2X0dXdZ,d2X0dZdX,d2X0dYdZ,d2X0dZdY,d2X0dX2,d2X0dY2,d2X0dZ2
REAL(dp) :: d2Y0dXdY,d2Y0dYdX,d2Y0dXdZ,d2Y0dZdX,d2Y0dYdZ,d2Y0dZdY,d2Y0dX2,d2Y0dY2,d2Y0dZ2
REAL(dp) :: d2Z0dXdY,d2Z0dYdX,d2Z0dXdZ,d2Z0dZdX,d2Z0dYdZ,d2Z0dZdY,d2Z0dX2,d2Z0dY2,D2Z0dZ2
REAL(dp) :: d2X0dtdX,d2X0dXdt,d2X0dtdY,d2X0dYdt,d2X0dtdZ,d2X0dZdt
REAL(dp) :: d2Y0dtdX,d2Y0dXdt,d2Y0dtdY,d2Y0dYdt,d2Y0dtdZ,d2Y0dZdt
REAL(dp) :: d2Z0dtdX,d2Z0dXdt,d2Z0dtdY,d2Z0dYdt,d2Z0dtdZ,d2Z0dZdt
REAL(dp) :: d2X0dt2, d2Y0dt2, d2Z0dt2
REAL(dp), DIMENSION(3) :: dlilbdx, dlilbdy, dlilbdz, dlilbdt, dlilbds, UEgradlilb, fpfield


 weightedh=0.0
 
 ! Steps is the number of loop top passes that we count. Different values shouldn't affect the Fermi
 ! contributions relative to each other and should only affect accuracy. Steps are uniformly spaced
 ! in time and use the loop top position at each of these times to estimate Fermi contribution

 ALLOCATE(carray(steps))
 ALLOCATE(ptsarray(looptopsteps,2)) ! Stores Fermi term and position at each point, useful for
                                    ! multiplying by average velocity later on
 ALLOCATE(uu(looptopsteps))
 ALLOCATE(lptppts(looptopsteps,3))

carray=0.0
ptsarray=0.0
uu=0.0
lptppts=0.0
beth=0.002!0.01!0.001
 
 fpfield = twoDODEint(fp,T)
 fpfs = MAGNITUDE3(fpfield)*bscl
! print*, 'fp = ', fp
! print*, 'fpfield = ', fpfield
! print*, 'fpfs = ', fpfs
 
 t0 = t0calc(fp, mu, tinit, initfs)
! print*, 't0 done'

 btol = Einit/mu ! Initial tolerance for when we hit mirror points
 gam = 1.0 + (mu*btol*1.602e-19)/(9.11e-31 * 3.0e8 * 3.0e8)


 IF (btol .gt. fpfs) then
   heightweighted=mu*Btol
   print*, 'particle escape at time =', Tinit*Tscl, 'energy = ', heightweighted
 else

!! betstep=betatronstep(fp, tinit, t0, mu, btol, beth)
 CALL betatronstep(fp, tinit, t0, mu, btol, beth, betstep, mirpointx)
 betstep=betstep/gam
! print*, 'betstep done'
 bettotal=betstep

 IF (abs(xinit) .le. looptoptol) THEN
   t0 = 0.0
   bettotal = 0.0
 END IF

 btol=btol+Einit/mu*bettotal ! Find mirror point and new btol at t=t0
 gam = 1.0 + (mu*btol*1.602e-19)/(9.11e-31 * 3.0e8 * 3.0e8)

 DO i = 0,steps-1-1 
   T=tinit+t0+i*(1.0-t0)/(steps-1.0)
   IF (btol .gt. fpfs) then
     heightweighted=mu*Btol
     print*, 'particle escape at time =', T*Tscl, 'energy = ', heightweighted
     exit
   END IF
   timeweight=timeinloops(fp,mu,T,min(looptoptol,mirpointx),btol) ! Find proportion of time spent in loops
!   print*, i, 'timeweight done'
! Find the points where we'll calculate e5
   CALL betatronstep(fp, T, (1.0-t0)/(steps-1), mu, btol, beth, betstep, mirpointx)
   lptppts=looptoppts(fp, t, min(looptoptol,mirpointx),looptopsteps) 
!   print*, i, 'looptoppts done'
!   betstep=betatronstep(fp, t, (1.0-t0)/(steps-1), mu, btol, beth)
   betstep=betstep/gam
   bettotal=bettotal+betstep
!   print*, i, 'betstep done'
   DO j = 0, looptopsteps-1

 ! find loop top position based on time and foot point position
 
  lptp=lptppts(j+1,:)

 ! find Fermi contribution at position and time related to looptop we've calculated

 R=lptp


 ! Field calculation leading up to finding e5 (curvature term in dgammadt)

 den=(a3d*a3d+R(1)*R(1)+ay*(R(2)-yh)*ay*(R(2)-yh)+R(3)*R(3))
 
oneuponden=1.0/(a3d*a3d+R(1)*R(1)+ay*(R(2)-yh)*ay*(R(2)-yh)+R(3)*R(3))

 Y0=(cc*T)**esp*log (1.+(R(2)/((cc*T)**esp)))*  &	
      (1.+tanh((R(2)-Lv/Lscl)*a))*0.5 +           &
       (1.-tanh((R(2)-Lv/Lscl)*a2))*0.5*R(2)

 X0=R(1)-delta*(Y0-R(2))*R(3)*oneuponden

 Z0=R(3)+delta*(Y0-R(2))*R(1)*oneuponden


 dY0dX=0.
 dY0dY=0.5*(1.+tanh((R(2)-Lv/Lscl)*a))/(1.+R(2)/(cc*T)**esp)+        &
      0.5*(cc*T)**esp*log(1.+R(2)/(cc*T)**esp)*               &
      (1.-tanh((R(2)-Lv/Lscl)*a)**2)*a-0.5*                         &
      (1.-tanh((R(2)-Lv/Lscl)*a2)**2)*a2*R(2)+0.5-0.5*tanh((R(2)-Lv/Lscl)*a2) 
 dY0dZ=0.

 
 dX0dX=1.0+2.0*R(1)*delta*(Y0-R(2))*R(3)*oneuponden*oneuponden - &
       delta*dY0dX*R(3)*oneuponden
 dX0dY=delta*(1.0-dY0dY)*R(3)*oneuponden+ &
       2.0*delta*(Y0-R(2))*ay*ay*(R(2)-yh)*R(3)*oneuponden*oneuponden
 dX0dZ=(-1.0)*delta*(Y0-R(2))*(a3d*a3d+R(1)*R(1)+ay*(R(2)-yh)*ay*(R(2)-yh)-R(3)*R(3))* &
                             oneuponden*oneuponden - &
       delta*dY0dZ*R(3)*oneuponden


 dZ0dX=delta*(Y0-R(2))*(a3d*a3d+ay*(R(2)-yh)*ay*(R(2)-yh)+R(3)*R(3)-R(1)*R(1))* &
                             oneuponden*oneuponden + &
       delta*dY0dX*R(1)*oneuponden
 dZ0dY=delta*(dY0dY-1.0)*R(1)*oneuponden-2.0*delta*(Y0-R(2))*ay*ay*(R(2)-yh)*R(1)*oneuponden*oneuponden
 dZ0dZ=1.0-2*delta*(Y0-R(2))*R(1)*R(3)*oneuponden*oneuponden + &
       delta*dY0dZ*R(1)*oneuponden
   
 dY0dt=0.5*(cc*T)**esp*esp*log(1.+R(2)/(cc*T)**esp)*          &
      (1.+tanh((R(2)-Lv/Lscl)*a))/T-0.5*R(2)*esp*                  &
      (1.+tanh((R(2)-Lv/Lscl)*a))/(T*(1.+R(2)/(cc*T)**esp)) 
 dX0dt=(-1.0)*delta*dY0dt*R(3)*oneuponden
 dZ0dt=delta*dY0dt*R(1)*oneuponden         


 d1 = SQRT((X0+0.5)*(X0+0.5)+(Y0+1)*(Y0+1)+Z0*Z0)
 d2 = SQRT((X0-0.5)*(X0-0.5)+(Y0+1)*(Y0+1)+Z0*Z0)

 B0(1) = (X0+0.5)/d1/d1/d1 - (X0-0.5)/d2/d2/d2
 B0(2) = (Y0+1)/d1/d1/d1 - (Y0+1)/d2/d2/d2
 B0(3) = Z0/d1/d1/d1 - Z0/d2/d2/d2

 B0 = B0 * c1

 dXdx(1)=dX0dX
 dXdx(2)=dY0dX
 dXdx(3)=dZ0dX

 dXdy(1)=dX0dY
 dXdy(2)=dY0dY
 dXdy(3)=dZ0dY

 dXdz(1)=dX0dZ
 dXdz(2)=dY0dZ
 dXdz(3)=dZ0dZ

 dB0dX0(1)=(d1*d1-3*(X0+0.5)*(X0+0.5))/d1/d1/d1/d1/d1 - (d2*d2-3*(X0-0.5)*(X0-0.5))/d2/d2/d2/d2/d2
 dB0dX0(2)=(-Y0-1)*3*(X0+0.5)/d1/d1/d1/d1/d1 - (-Y0-1)*3*(X0-0.5)/d2/d2/d2/d2/d2
 dB0dX0(3)=(-Z0)*3*(X0+0.5)/d1/d1/d1/d1/d1 - (-Z0)*3*(X0-0.5)/d2/d2/d2/d2/d2

 dB0dY0(1)=(-X0-0.5)*3*(Y0+1)/d1/d1/d1/d1/d1 - (-X0+0.5)*3*(Y0+1)/d2/d2/d2/d2/d2
 dB0dY0(2)=(d1*d1-3*(Y0+1)*(Y0+1))/d1/d1/d1/d1/d1 - (d2*d2-3*(Y0+1)*(Y0+1))/d2/d2/d2/d2/d2
 dB0dY0(3)=(-Z0)*3*(Y0+1)/d1/d1/d1/d1/d1 - (-Z0)*3*(Y0+1)/d2/d2/d2/d2/d2

 dB0dZ0(1)=(-X0-0.5)*3*Z0/d1/d1/d1/d1/d1 - (-X0+0.5)*3*Z0/d2/d2/d2/d2/d2
 dB0dZ0(2)=(-Y0-1)*3*Z0/d1/d1/d1/d1/d1 - (-Y0-1)*3*Z0/d2/d2/d2/d2/d2
 dB0dZ0(3)=(d1*d1-3*Z0*Z0)/d1/d1/d1/d1/d1 - (d2*d2-3*Z0*Z0)/d2/d2/d2/d2/d2
 
 !km: Again, multiply by c1 factor in Grady(2009) eq (36)
 dB0dX0 = dB0dX0 * c1
 dB0dY0 = dB0dY0 * c1
 dB0dZ0 = dB0dZ0 * c1


 y=R(2)
 
 d2y0dy2= -0.1e1 / (0.1e1 + y / (cc * t) ** esp) ** 2 * (0.1e1 + tan&
     &h((y - Lv/Lscl) * a)) / (cc * t) ** esp / 0.2e1 + 0.1e1 / (0.1e1 + y /&
     & (cc * t) ** esp) * (0.1e1 - tanh((y - Lv/Lscl) * a) ** 2) * a - (cc& 
     &* t) ** esp * log(0.1e1 + y / (cc * t) ** esp) * tanh((y - Lv/Lscl) * a&
     &) * (0.1e1 - tanh((y - Lv/Lscl) * a) ** 2) * a*a + tanh((y - Lv/Lscl)& 
     &* a2) * (0.1e1 - tanh((y - Lv/Lscl) * a2) ** 2) * a2*a2 * y - (0.1e1& 
     &- tanh((y - Lv/Lscl) * a2) ** 2) * a2
 d2Y0dXdY=0.
 d2Y0dYdX=0.
 d2Y0dXdZ=0.
 d2Y0dZdX=0.
 d2Y0dYdZ=0.
 d2Y0dZdY=0.
 d2Y0dX2=0.
 d2Y0dZ2=0.

 d2X0dXdY = 2*delta*(dY0dY-1.0)*R(1)*R(3)*oneuponden*oneuponden - &
      8*delta*(Y0-R(2))*R(1)*R(3)*ay*ay*(R(2)-yh)*oneuponden*oneuponden*oneuponden + &
      delta*dY0dX*2*ay*ay*(R(2)-yh)*R(3)*oneuponden*oneuponden-delta*d2Y0dXdY*R(3)*oneuponden
 d2X0dYdX = d2X0dXdY
 d2X0dXdZ = 2*delta*(Y0-R(2))*R(1)*(den-4*R(3)*R(3))*oneuponden*oneuponden*oneuponden + &
            delta*dY0dZ*2*R(1)*R(3)*oneuponden*oneuponden-delta*d2Y0dXdZ*R(3)*oneuponden - &
            delta*dY0dX*(den-2*R(3)*R(3))*oneuponden*oneuponden
 d2X0dZdX = d2X0dXdZ
 d2X0dYdZ = (-delta)*(dY0dY-1.0)*(den-2*R(3)*R(3))*oneuponden*oneuponden + &
            2*delta*(Y0-R(2))*ay*ay*(R(2)-yh)*(den-4*R(3)*R(3))*oneuponden*oneuponden*oneuponden + &
            delta*dY0dZ*2*ay*ay*(R(2)-yh)*R(3)*oneuponden*oneuponden-delta*d2Y0dYdZ*R(3)*oneuponden
 d2X0dZdY = d2X0dYdZ
 d2X0dX2 = 2*delta*(Y0-R(2))*R(3)*(den-4*R(1)*R(1))*oneuponden*oneuponden*oneuponden + &
           2*delta*dY0dX*2*R(1)*R(3)*oneuponden*oneuponden-delta*d2Y0dX2*R(3)*oneuponden
 d2X0dY2 = (-delta)*d2Y0dY2*R(3)*oneuponden+4*delta*(dY0dY-1.0)*ay*ay*(R(2)-yh)*R(3)*oneuponden*oneuponden + &
           2*delta*(Y0-R(2))*R(3)*(ay*ay*den-4*ay*ay*(R(2)-yh)*ay*ay*(R(2)-yh))*oneuponden*oneuponden*oneuponden
 d2X0dZ2 = delta*(Y0-R(2))*(4*R(3)*(a3d*a3d+R(1)*R(1)+ay*(R(2)-yh)*ay*(R(2)-yh)-R(3)*R(3))+ &
           2*R(3)*den)*oneuponden*oneuponden*oneuponden - &
           2*delta*dY0dZ*(den-2*R(3)*R(3))*oneuponden*oneuponden-delta*d2Y0dZ2*R(3)*oneuponden

 d2Z0dXdY = delta*(dY0dY-1.0)*(den-2*R(1)*R(1))*oneuponden*oneuponden + &
            delta*(Y0-R(2))*(2*ay*ay*(R(2)-yh)*den-4*ay*ay*(R(2)-yh)* &
            (a3d*a3d+ay*(R(2)-yh)*ay*(R(2)-yh)+R(3)*R(3)-R(1)*R(1)))*oneuponden*oneuponden*oneuponden + &
            delta*d2Y0dXdY*R(1)*oneuponden-delta*dY0dX*2*ay*ay*(R(2)-yh)*R(1)*oneuponden*oneuponden
 d2Z0dYdX = d2Z0dXdY
 d2Z0dXdZ = delta*(Y0-R(2))*R(3)*(2*den-4*(a3d*a3d+ay*(R(2)-yh)*ay*(R(2)-yh)+R(3)*R(3)-R(1)*R(1)))* &
            oneuponden*oneuponden*oneuponden + &
            delta*d2Y0dXdZ*R(1)*oneuponden+delta*dY0dZ*(den-2*R(1)*R(1))*oneuponden*oneuponden - &
            delta*dY0dX*2*R(1)*R(3)*oneuponden*oneuponden
 d2Z0dZdX = d2Z0dXdZ
 d2Z0dYdZ = (-delta)*(dY0dY-1.0)*2*R(3)*R(1)*oneuponden*oneuponden + &
            delta*(Y0-R(2))*8*R(1)*R(3)*ay*ay*(R(2)-yh)*oneuponden*oneuponden*oneuponden + &
            delta*d2Y0dYdZ*R(1)*oneuponden-delta*dY0dZ*2*ay*ay*(R(2)-yh)*R(1)*oneuponden*oneuponden
 d2Z0dZdY = d2Z0dYdZ
 d2Z0dX2 = (-delta)*(Y0-R(2))*(2*R(1)*den+4*R(1)* &
           (a3d*a3d+ay*(R(2)-yh)*ay*(R(2)-yh)+R(3)*R(3)-R(1)*R(1)))*oneuponden*oneuponden*oneuponden + &
           delta*d2Y0dX2*R(1)*oneuponden+2*delta*dY0dX*(den-2*R(1)*R(1))*oneuponden*oneuponden
 d2Z0dY2 = delta*d2Y0dY2*R(1)*oneuponden-4*delta*(dY0dY-1.0)*R(1)*ay*ay*(R(2)-yh)*oneuponden*oneuponden + &
           delta*(Y0-R(2))*R(1)*(8*ay*ay*(R(2)-yh)*ay*ay*(R(2)-yh)-2*ay*ay*den)*oneuponden*oneuponden*oneuponden
 d2Z0dZ2 = 2*(-delta)*(Y0-R(2))*R(1)*(den-4*R(3)*R(3))*oneuponden*oneuponden*oneuponden + &
           delta*d2Y0dZ2*R(1)*oneuponden-2*delta*dY0dZ*2*R(1)*R(3)*oneuponden*oneuponden

d2y0dtdx =0.
d2y0dxdt=d2y0dtdx
d2y0dtdz =0.
d2y0dzdt =d2y0dtdz

d2y0dtdy = 0.1e1 / (0.1e1 + y / (cc * t) ** esp) ** 2 * (0.1e1 + tanh &
     ((y - Lv/Lscl) * a)) * y / (cc * t) ** esp * esp / t / 0.2e1 + (cc * & 
     t) ** esp * esp / t * log(0.1e1 + y / (cc * t) ** esp) * (0.1e1 - &
     tanh((y - Lv/Lscl) * a) ** 2) * a / 0.2e1 - y * esp / t / (0.1e1 + y &
     / (cc * t) ** esp) * (0.1e1 - tanh((y - Lv/Lscl) * a) ** 2) * a / 0.2e1

d2y0dydt=d2y0dtdy


d2y0dt2 = (cc * t) ** esp * esp ** 2 / t ** 2 * log(0.1e1 + y / (cc &
      * t) ** esp) * (0.1e1 + tanh((y - Lv/Lscl) * a)) / 0.2e1 - (cc * t) ** &
      esp * esp / t ** 2 * log(0.1e1 + y / (cc * t) ** esp) * (0.1e1 +  &
      tanh((y - Lv/Lscl) * a)) / 0.2e1 - esp ** 2 / t ** 2 * y / (0.1e1 + y  &
      / (cc * t) ** esp) * (0.1e1 + tanh((y - Lv/Lscl) * a)) / 0.2e1 + y *   &
     esp / t ** 2 / (0.1e1 + y / (cc * t) ** esp) * (0.1e1 + tanh((y -  &
     Lv/Lscl) * a)) / 0.2e1 - y ** 2 * esp ** 2 / t ** 2 / (0.1e1 + y / (cc  &
      * t) ** esp) ** 2 * (0.1e1 + tanh((y - Lv/Lscl) * a)) / (cc * t) **    &
     esp / 0.2e1

d2X0dtdx = 2*delta*dY0dt*R(1)*R(3)*oneuponden*oneuponden-delta*d2Y0dXdt*R(3)*oneuponden
d2X0dxdt = d2X0dtdx
d2X0dydt = (-delta)*d2Y0dtdy*R(3)*oneuponden+2*delta*dY0dt*ay*ay*(R(2)-yh)*R(3)*oneuponden*oneuponden
d2X0dtdy = d2X0dydt
d2X0dtdz = (-delta)*dY0dt*(den-2*R(3)*R(3))*oneuponden*oneuponden-delta*d2Y0dzdt*R(3)*oneuponden
d2X0dzdt = d2X0dtdz
d2X0dt2 = (-delta)*d2y0dt2*R(3)*oneuponden


d2Z0dtdx = delta*dY0dt*(den-2*R(1)*R(1))*oneuponden*oneuponden+delta*d2Y0dXdt*R(1)*oneuponden
d2Z0dxdt = d2Z0dtdx
d2Z0dtdy = delta*d2Y0dydt*R(1)*oneuponden-2*delta*dY0dt*R(1)*ay*ay*(R(2)-yh)*oneuponden*oneuponden
d2Z0dydt = d2Z0dtdy
d2Z0dtdz = (-delta)*dY0dt*2*R(1)*R(3)*oneuponden*oneuponden+delta*d2Y0dZdt*R(1)*oneuponden
d2Z0dzdt = d2Z0dtdz
d2Z0dt2 = delta*d2Y0dt2*R(1)*oneuponden

  d2Xdxdt(1)=d2X0dxdt
  d2Xdxdt(2)=d2Y0dxdt
  d2Xdxdt(3)=d2Z0dxdt

  d2Xdydt(1)=d2X0dydt
  d2Xdydt(2)=d2Y0dydt
  d2Xdydt(3)=d2Z0dydt

  d2Xdzdt(1)=d2X0dzdt
  d2Xdzdt(2)=d2Y0dzdt
  d2Xdzdt(3)=d2Z0dzdt

  d2Xdx2(1)=d2X0dx2
  d2Xdx2(2)=d2Y0dx2
  d2Xdx2(3)=d2Z0dx2

  d2Xdydx(1)=d2X0dydx
  d2Xdydx(2)=d2Y0dydx
  d2Xdydx(3)=d2Z0dydx

  d2Xdzdx(1)=d2X0dzdx
  d2Xdzdx(2)=d2Y0dzdx
  d2Xdzdx(3)=d2Z0dzdx

  d2Xdxdy(1)=d2X0dxdy
  d2Xdxdy(2)=d2Y0dxdy
  d2Xdxdy(3)=d2Z0dxdy

  d2Xdy2(1)=d2X0dy2
  d2Xdy2(2)=d2Y0dy2
  d2Xdy2(3)=d2Z0dy2

  d2Xdzdy(1)=d2X0dzdy
  d2Xdzdy(2)=d2Y0dzdy
  d2Xdzdy(3)=d2Z0dzdy

  d2Xdxdz(1)=d2X0dxdz
  d2Xdxdz(2)=d2Y0dxdz
  d2Xdxdz(3)=d2Z0dxdz

  d2Xdydz(1)=d2X0dydz
  d2Xdydz(2)=d2Y0dydz
  d2Xdydz(3)=d2Z0dydz

  d2Xdz2(1)=d2X0dz2
  d2Xdz2(2)=d2Y0dz2
  d2Xdz2(3)=d2Z0dz2

  dB0dx = dX0dx*dB0dx0 + dY0dx*dB0dy0 + dZ0dx*dB0dz0
  dB0dy = dX0dy*dB0dx0 + dY0dy*dB0dy0 + dZ0dy*dB0dz0
  dB0dz = dX0dz*dB0dx0 + dY0dz*dB0dy0 + dZ0dz*dB0dz0
  dB0dt = dX0dt*dB0dx0 + dY0dt*dB0dy0 + dZ0dt*dB0dz0

 !dBxdx
 dBdx(1)=oneuponL*oneuponL*(DOT(CROSS(d2Xdydx,dXdz)+CROSS(dXdy,d2Xdzdx),B0)+&
         DOT(CROSS(dXdy,dXdz),dB0dx))
 !dBxdy
 dBdy(1)=oneuponL*oneuponL*(DOT(CROSS(d2Xdy2,dXdz)+CROSS(dXdy,d2Xdzdy),B0)+&
         DOT(CROSS(dXdy,dXdz),dB0dy))
 !dBxdz
 dBdz(1)=oneuponL*oneuponL*(DOT(CROSS(d2Xdydz,dXdz)+CROSS(dXdy,d2Xdz2),B0)+&
         DOT(CROSS(dXdy,dXdz),dB0dz))
 !dBydx
 dBdx(2)=oneuponL*oneuponL*(DOT(CROSS(d2Xdzdx,dXdx)+CROSS(dXdz,d2Xdx2),B0)+&
         DOT(CROSS(dXdz,dXdx),dB0dx))
 !dBydy
 dBdy(2)=oneuponL*oneuponL*(DOT(CROSS(d2Xdzdy,dXdx)+CROSS(dXdz,d2Xdxdy),B0)+&
         DOT(CROSS(dXdz,dXdx),dB0dy))
 !dBydz
 dBdz(2)=oneuponL*oneuponL*(DOT(CROSS(d2Xdz2,dXdx)+CROSS(dXdz,d2Xdxdz),B0)+&
         DOT(CROSS(dXdz,dXdx),dB0dz))
 !dBzdx
 dBdx(3)=oneuponL*oneuponL*(DOT(CROSS(d2Xdx2,dXdy)+CROSS(dXdx,d2Xdydx),B0)+&
         DOT(CROSS(dXdx,dXdy),dB0dx)) 
 !dBzdy
 dBdy(3)=oneuponL*oneuponL*(DOT(CROSS(d2Xdxdy,dXdy)+CROSS(dXdx,d2Xdy2),B0)+&
         DOT(CROSS(dXdx,dXdy),dB0dy)) 

 !dBzdz
 dBdz(3)=oneuponL*oneuponL*(DOT(CROSS(d2Xdxdz,dXdy)+CROSS(dXdx,d2Xdydz),B0)+&
         DOT(CROSS(dXdx,dXdy),dB0dz)) 

 !dBxdt
 dBdt(1)=oneuponTL*(DOT(CROSS(d2Xdydt,dXdz)+CROSS(dXdy,d2Xdzdt),B0)+&
         DOT(CROSS(dXdy,dXdz),dB0dt))
 !dBydt
 dBdt(2)=oneuponTL*(DOT(CROSS(d2Xdzdt,dXdx)+CROSS(dXdz,d2Xdxdt),B0)+&
         DOT(CROSS(dXdz,dXdx),dB0dt))
 !dBzdt
 dBdt(3)=oneuponTL*(DOT(CROSS(d2Xdxdt,dXdy)+CROSS(dXdx,d2Xdydt),B0)+&
         DOT(CROSS(dXdx,dXdy),dB0dt)) 


 ! Make derivatives dimensionless
 dBdx = lscl/bscl*dBdx
 dBdy = lscl/bscl*dBdy
 dBdz = lscl/bscl*dBdz
 dBdt = tscl/bscl*dBdt



 B(1)= DOT(CROSS(dXdy,dXdz),B0)*oneuponL
 B(2)= DOT(CROSS(dXdz,dXdx),B0)*oneuponL
 B(3)= DOT(CROSS(dXdx,dXdy),B0)*oneuponL

 B=B/bscl

 ! Estimate the parallel velocity, we take the velocity at the first theoretical loop top pass
 ! if energy doesn't change. Not necessarily accurate , but captures how particles will have 
 ! different parallel velocities at the loop top relating to initial pitch angle and position.
 ! Updated to account for energy changes


 U=SQRT(2*(mu*(btol-MAGNITUDE3(B)*bscl))*1.6e-19/9.1e-31)/lscl*tscl*gam

 if (btol-MAGNITUDE3(B)*bscl .LE. 0.0) then 
   u = 0.0
 end if

! if (i .eq. 0) then
!!  print*, MAGNITUDE3(B)*bscl, btol-MAGNITUDE3(B)*bscl, mu*(btol-MAGNITUDE3(B)*bscl)
!  print*, u
! end if

 if (j .eq. looptopsteps-1) then
!! print*, 'btol is ', btol, 'R is ', R
!! print*, 'U is', U
 u = ptsarray(looptopsteps-1,2)
 end if

 ! Make sure velocity isn't square root of negative number, necessary is mirror points are in 
 ! acceleration region

 ptsarray(j+1,2)=U


 omodb = 1.0/SQRT(B(1)*B(1)+B(2)*B(2)+B(3)*B(3))

 GRADB(1) = DOT(B,dBdX)*omodb
 GRADB(2) = DOT(B,dBdY)*omodb
 GRADB(3) = DOT(B,dBdZ)*omodb
 GRADBT = DOT(B,dBdT)*omodb

 dlilbdx = dBdX*omodb - B*GRADB(1)*omodb*omodb
 dlilbdy = dBdY*omodb - B*GRADB(2)*omodb*omodb
 dlilbdz = dBdZ*omodb - B*GRADB(3)*omodb*omodb
 dlilbdt = dBdT*omodb - B*GRADBT*omodb*omodb


 DETERMINANT = dX0dx*(dY0dy*dZ0dz-dZ0dy*dY0dz)-dX0dy*(dY0dx*dZ0dz-dZ0dx*dY0dz)+&
               dX0dz*(dY0dx*dZ0dy-dZ0dx*dY0dy)
 
 Vf(1)= (-1.0)*(dX0dt*(dY0dy*dZ0dz-dZ0dy*dY0dz) +   &
                 dY0dt*(dX0dz*dZ0dy-dX0dy*dZ0dz) +   &
                 dZ0dt*(dX0dy*dY0dz-dX0dz*dY0dy) )/DETERMINANT
 
 Vf(2)= (-1.0)*(dX0dt*(dZ0dx*dY0dz-dZ0dz*dY0dx) +   &
                 dY0dt*(dX0dx*dZ0dz-dX0dz*dZ0dx) +   &
                 dZ0dt*(dX0dz*dY0dx-dX0dx*dY0dz) )/DETERMINANT

 Vf(3)= (-1.0)*(dX0dt*(dY0dx*dZ0dy-dY0dy*dZ0dx) +   &
                 dY0dt*(dX0dy*dZ0dx-dX0dx*dZ0dy) +   &
                 dZ0dt*(dX0dx*dY0dy-dX0dy*dY0dx) )/DETERMINANT

 E = -CROSS(vf,B)
 
 UE = CROSS(E,B)*oMODB*oMODB

 dlilbds = (B(1)*dlilbdx + B(2)*dlilbdy + B(3)*dlilbdz)*omodb

 UEgradlilb = UE(1)*dlilbdx + UE(2)*dlilbdy + UE(3)*dlilbdz


 dlilbdt = U*dlilbds
 ! Only have 1 term here because the other two cancel due to loop symmetry

 curvature = U*DOT(dlilbdt,UE)
! print, 'step = ', j, 'n=',n, 'curvature = ', curvature

 if (j .EQ. looptopsteps-1) then 
   curvature=curvature/2.0
 end if
 ! We halve the contribution at the very loop top, we've used the symmetry of the loop so we
 ! want to double count everywhere but here

 ptsarray(j+1,1)=curvature*timeweight

 end do

uaverage = (SUM(ptsarray(:,2))-0.5*ptsarray(looptopsteps,2))/(SIZE(ptsarray(:,2))-0.5)
!print, uaverage
 zcount=0
 curveterm=0.0
! We need to multiply the contirbutions by the average speed and divide by the speed at that point.
! The points are spaced equally but particles will spend less time in areas where they move faster.
! We need to weight contributions by the time spent in certain areas in the loop top.

DO j=0,looptopsteps-1
  if (ptsarray(j+1,2) .ne. 0.0) then 
    curveterm=curveterm+ptsarray(j+1,1)*uaverage/ptsarray(j+1,2)
  else 
   zcount=zcount+1
   ! Count the points where velocity is zero, we don't want these to count towards average speed
  end if
end do

 if (i .eq. 0 ) then

!! PRINT*, 'fermi terms = ', ptsarray(:,1)
! PRINT*, 'zcount = ', zcount
! PRINT*, SUM(ptsarray(:,2))-0.5*ptsarray(looptopsteps,2), SUM(ptsarray(:,2)), ptsarray(looptopsteps,2)
! print*, PTSARRAY(:,2)
!! PRINT*, 'velocities = ', ptsarray(:,2)
! PRINT*, 'uaverage = ', uaverage*(SIZE(ptsarray(:,2))-0.5)/(SIZE(ptsarray(:,2))-0.5-zcount)
! PRINT*, 'correcting term = ', (SIZE(ptsarray(:,2))-0.5)/(SIZE(ptsarray(:,2))-0.5-zcount)
! PRINT*, 'curveterm = ', curveterm*(SIZE(ptsarray(:,2))-0.5)/(SIZE(ptsarray(:,2))-0.5-zcount)*&
!                                   (SIZE(ptsarray(:,2))-0.5)/(SIZE(ptsarray(:,2))-0.5-zcount)/&
!                                   (looptopsteps-0.5)*1e10/1.6e-19*9.1e-31/steps*(1.0-t0)
endif

 curveterm=curveterm*&
                  (SIZE(ptsarray(:,2))-0.5)/(SIZE(ptsarray(:,2))-0.5-zcount)*&
                  (SIZE(ptsarray(:,2))-0.5)/(SIZE(ptsarray(:,2))-0.5-zcount)
! We multiply by the square of the fraction of points with non zero velocities. The first 
! multiplication accounts for u average using points with zero velocity for its calculation.
! The second multiplication accounts for Fermi contributions where u=0 being calculated abouve,
! effectively creating an average that includes points that the particle never crosses.

!print, uaverage*FLOAT(SIZE(ptsarray[*,1])/(SIZE(ptsarray[*,1])-zcount))

carray(steps-i)=carray(steps-i)+curveterm/gam
! Add the Fermi contribution from this loop

if (abs(mirpointx) .lt. 0.1) then
   carray(steps- i) = carray(steps - i)*abs(mirpointx)/0.1!*abs(mirpointx)/0.15
end if 

 fenterm = SUM(carray)/(looptopsteps-0.5)*1e10/1.6e-19*9.1e-31/steps*(1.0-t0)

!!IF (i .eq. 1) then
!1  print*, fenterm
!!end if

 ! This term calculates Fermi acceleration so far, used for updated velocity calculation

!! print*, 'fermi = ', fenterm
!! print*, 'betatron - ', Einit*bettotal

 btol=fenterm/mu+bettotal*Einit/mu+Einit/mu
 gam = 1.0 + (mu*btol*1.602e-19)/(9.11e-31 * 3.0e8 * 3.0e8)

!print*, 't = ', t, 'energy = ', btol*mu

!print, 't =', t, ' energy =', mu*btol, ' zeros =', zcount
! print, 'Bterm = ', exbvector[3]
! print, 'fenterm = ', fenterm/mu
! print, 't = ',t,'energy = ', mu*btol
! Find a new Btol from the movement of the mirror point between this loop calculation and the next

end do

! Add up contributions for all loop calculations

!DO i = 0, steps-1
!  weightedh = weightedh + carray[steps-i-1]
!end do

!print,'fermi term =', weightedh/(looptopsteps-0.5)*1e10/1.6e-19*9.1e-31*1000/steps

!print, weightedh/(looptopsteps-0.5)*100*1e10/9e16*1000/steps

! print, 'est = ', mu*btol

carray=0.0
ptsarray=0.0
uu=0.0
lptppts=0.0

 DEALLOCATE(carray)
 DEALLOCATE(ptsarray) 
 DEALLOCATE(uu)
 DEALLOCATE(lptppts)

 heightweighted = mu*btol

 endif
! Multiply the dgamma/dt term by vscl^2 and mass of electron, divide by charge for eV, x1000 to
! correct earlier division by 1000.

! return,weightedh/(looptopsteps-0.5)*100*1e10/9e16/steps

! print*, 'fenterm = ', fenterm

END FUNCTION heightweighted

FUNCTION ExB(R,T)

REAL(dp), DIMENSION(3) :: R, ExB
REAL(dp), DIMENSION(3) :: b, E, vf, dXdx, dXdy, dXdz, B0
REAL(dp) :: tstep, T, X0, Y0, Z0, omodb, dX0dX, dX0dY, dX0dZ, dY0dX, dY0dY, dY0dZ, dZ0dX, dZ0dY, dZ0dZ
REAL(dp) :: dX0dt, dY0dt, dZ0dt, determinant, d1, d2, oneuponden

oneuponden=1.0/(a3d*a3d+R(1)*R(1)+ay*(R(2)-yh)*ay*(R(2)-yh)+R(3)*R(3))

 Y0=(cc*T)**esp*log (1.+(R(2)/((cc*T)**esp)))*  &	
      (1.+tanh((R(2)-Lv/Lscl)*a))*0.5 +           &
       (1.-tanh((R(2)-Lv/Lscl)*a2))*0.5*R(2)

 X0=R(1)-delta*(Y0-R(2))*R(3)*oneuponden

 Z0=R(3)+delta*(Y0-R(2))*R(1)*oneuponden


 dY0dX=0.
 dY0dY=0.5*(1.+tanh((R(2)-Lv/Lscl)*a))/(1.+R(2)/(cc*T)**esp)+        &
      0.5*(cc*T)**esp*log(1.+R(2)/(cc*T)**esp)*               &
      (1.-tanh((R(2)-Lv/Lscl)*a)**2)*a-0.5*                         &
      (1.-tanh((R(2)-Lv/Lscl)*a2)**2)*a2*R(2)+0.5-0.5*tanh((R(2)-Lv/Lscl)*a2) 
 dY0dZ=0.

 
 dX0dX=1.0+2.0*R(1)*delta*(Y0-R(2))*R(3)*oneuponden*oneuponden - &
       delta*dY0dX*R(3)*oneuponden
 dX0dY=delta*(1.0-dY0dY)*R(3)*oneuponden+ &
       2.0*delta*(Y0-R(2))*ay*ay*(R(2)-yh)*R(3)*oneuponden*oneuponden
 dX0dZ=(-1.0)*delta*(Y0-R(2))*(a3d*a3d+R(1)*R(1)+ay*(R(2)-yh)*ay*(R(2)-yh)-R(3)*R(3))* &
                             oneuponden*oneuponden - &
       delta*dY0dZ*R(3)*oneuponden


 dZ0dX=delta*(Y0-R(2))*(a3d*a3d+ay*(R(2)-yh)*ay*(R(2)-yh)+R(3)*R(3)-R(1)*R(1))* &
                             oneuponden*oneuponden + &
       delta*dY0dX*R(1)*oneuponden
 dZ0dY=delta*(dY0dY-1.0)*R(1)*oneuponden-2.0*delta*(Y0-R(2))*ay*ay*(R(2)-yh)*R(1)*oneuponden*oneuponden
 dZ0dZ=1.0-2*delta*(Y0-R(2))*R(1)*R(3)*oneuponden*oneuponden + &
       delta*dY0dZ*R(1)*oneuponden
   
 dY0dt=0.5*(cc*T)**esp*esp*log(1.+R(2)/(cc*T)**esp)*          &
      (1.+tanh((R(2)-Lv/Lscl)*a))/T-0.5*R(2)*esp*                  &
      (1.+tanh((R(2)-Lv/Lscl)*a))/(T*(1.+R(2)/(cc*T)**esp)) 
 dX0dt=(-1.0)*delta*dY0dt*R(3)*oneuponden
 dZ0dt=delta*dY0dt*R(1)*oneuponden         


 d1 = SQRT((X0+0.5)*(X0+0.5)+(Y0+1)*(Y0+1)+Z0*Z0)
 d2 = SQRT((X0-0.5)*(X0-0.5)+(Y0+1)*(Y0+1)+Z0*Z0)

 B0(1) = (X0+0.5)/d1/d1/d1 - (X0-0.5)/d2/d2/d2
 B0(2) = (Y0+1)/d1/d1/d1 - (Y0+1)/d2/d2/d2
 B0(3) = Z0/d1/d1/d1 - Z0/d2/d2/d2

 B0 = B0 * c1

 dXdx(1)=dX0dX
 dXdx(2)=dY0dX
 dXdx(3)=dZ0dX

 dXdy(1)=dX0dY
 dXdy(2)=dY0dY
 dXdy(3)=dZ0dY

 dXdz(1)=dX0dZ
 dXdz(2)=dY0dZ
 dXdz(3)=dZ0dZ

 B(1)= DOT(CROSS(dXdy,dXdz),B0)*oneuponL
 B(2)= DOT(CROSS(dXdz,dXdx),B0)*oneuponL
 B(3)= DOT(CROSS(dXdx,dXdy),B0)*oneuponL

 B = B/bscl

 omodb = 1.0/SQRT(B(1)*B(1)+B(2)*B(2)+B(3)*B(3))


 DETERMINANT = dX0dx*(dY0dy*dZ0dz-dZ0dy*dY0dz)-dX0dy*(dY0dx*dZ0dz-dZ0dx*dY0dz)+&
               dX0dz*(dY0dx*dZ0dy-dZ0dx*dY0dy)
 
 Vf(1)= (-1.0)*(dX0dt*(dY0dy*dZ0dz-dZ0dy*dY0dz) +   &
                 dY0dt*(dX0dz*dZ0dy-dX0dy*dZ0dz) +   &
                 dZ0dt*(dX0dy*dY0dz-dX0dz*dY0dy) )/DETERMINANT
 
 Vf(2)= (-1.0)*(dX0dt*(dZ0dx*dY0dz-dZ0dz*dY0dx) +   &
                 dY0dt*(dX0dx*dZ0dz-dX0dz*dZ0dx) +   &
                 dZ0dt*(dX0dz*dY0dx-dX0dx*dY0dz) )/DETERMINANT

 Vf(3)= (-1.0)*(dX0dt*(dY0dx*dZ0dy-dY0dy*dZ0dx) +   &
                 dY0dt*(dX0dy*dZ0dx-dX0dx*dZ0dy) +   &
                 dZ0dt*(dX0dx*dY0dy-dX0dy*dY0dx) )/DETERMINANT

 E = -CROSS(vf,B)
 
 ExB = CROSS(E,B)*oMODB*oMODB

END FUNCTION ExB

SUBROUTINE driftevo(R, tstart, tsteps, trange, ffs)

 ! Returns the position and field strength after the ExB drift has transported a point for a 
 ! specified time

IMPLICIT NONE

REAL(dp), DIMENSION(3), INTENT(INOUT) :: R
REAL(dp), INTENT(IN) :: tstart, trange
INTEGER, INTENT(IN) :: tsteps
REAL(dp), INTENT(OUT) :: ffs
REAL(dp), DIMENSION(3) :: finalfield, k1,k2,k3,k4
REAL(dp) :: tstep, T
INTEGER :: i

 tstep = trange/tsteps

!! print*, '*************', tstart, tstart+trange, R


DO i = 0, tsteps-1 

 T=tstart+tstep*i

   k1 = ExB(R,T)
   k2 = ExB(R+k1*tstep/2.0,T+tstep/2.0)
   k3 = ExB(R+k2*tstep/2.0,T+tstep/2.0)
   k4 = ExB(R+k3*tstep,T+tstep)
   R = R + tstep*(k1+2*k2+2*k3+k4)/6.0

!! print*, T, R

END DO

finalfield=twoDODEint(R,T+tstep)*bscl

ffs=MAGNITUDE3(finalfield)

!!print*, 'END', T+tstep, R


END SUBROUTINE driftevo

END PROGRAM fermiguessfor
