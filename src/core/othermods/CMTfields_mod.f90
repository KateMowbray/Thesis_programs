MODULE CMT_fields

! km: Adjusted to cover 3D case, not a lengthy process to comment/uncomment a few lines to revert this
!     to 2D or 2.5D but CMT_fields2half is already set up this way to avoid the extra effort
! km: The subroutines called to find the derivatives of X0, Y0 and Z0 can be changed so that this 
!     runs the braking jet case instead of the simpled CMT 
!     (this is called ABSETUP and replaces (SUB1_X0Y0Z0, SUB2_X0Y0Z0 and SUB3_X0Y0Z0)

  USE global
  USE M_products

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: CMTFIELDS
  PRIVATE :: SUB1_X0Y0Z0, dA0, SUB2_X0Y0Z0,ddA0, SUB3_X0Y0Z0, B0setup,dB0
  
  CONTAINS 
!-----------------------------------------------------------------------------!   
SUBROUTINE CMTFIELDS(R,T,E,B,DBDX,DBDY,DBDZ,DBDT,DEDX,DEDY,DEDZ,DEDT,Vf)

 REAL(num)	:: X0,Y0,Z0,dX0dX,dX0dY,dX0dZ,dY0dX,dY0dY,dY0dZ,dZ0dX,dZ0dY,dZ0dZ,dX0dt,dY0dt,dZ0dt
 REAL(num)	:: d2X0dXdY,d2X0dYdX,d2X0dXdZ,d2X0dZdX,d2X0dYdZ,d2X0dZdY,d2X0dX2,d2X0dY2,d2X0dZ2
 REAL(num)	:: d2Y0dXdY,d2Y0dYdX,d2Y0dXdZ,d2Y0dZdX,d2Y0dYdZ,d2Y0dZdY,d2Y0dX2,d2Y0dY2,d2Y0dZ2
 REAL(num)      :: d2Z0dXdY,d2Z0dYdX,d2Z0dXdZ,d2Z0dZdX,d2Z0dYdZ,d2Z0dZdY,d2Z0dX2,d2Z0dY2,d2Z0dZ2
 REAL(num)	:: d2X0dtdX, d2X0dXdt, d2X0dtdY, d2X0dYdt, d2X0dtdZ, d2X0dZdt
 REAL(num)	:: d2Y0dtdX, d2Y0dXdt, d2Y0dtdY, d2Y0dYdt, d2Y0dtdZ, d2Y0dZdt
 REAL(num)      :: d2Z0dtdX, d2Z0dXdt, d2Z0dtdY, d2Z0dYdt, d2Z0dtdZ, d2Z0dZdt
 REAL(num)	:: d2X0dt2, d2Y0dt2, d2Z0dt2
 REAL(num)	:: dA0dX0, dA0dY0, dA0dZ0
 REAL(num)	:: d2A0dx02,d2A0dy02,d2A0dx0dy0,d2A0dy0dx0
 REAL(num)	:: DETERMINANT,der_det,xder_det,yder_det,zder_det
 REAL(num), DIMENSION(3), INTENT(OUT) :: B,E
 REAL(num), DIMENSION(3), INTENT(OUT) :: DBDX,DBDY,DBDZ,DBDT,DEDX,DEDY,DEDZ,DEDT
 REAL(num), DIMENSION(3), INTENT(IN) :: R
 REAL(num), DIMENSION(3)  :: dB0dx0, dB0dy0, dB0dz0
 REAL(num), DIMENSION(3)  :: Vf,dVfdt,dVfdx,dVfdy,dVfdz,B0,X,dXdx,dXdy,dXdz,dXdt,dB0dt,dB0dx,dB0dy,dB0dz
 REAL(num), DIMENSION(3)  :: d2Xdxdt,d2Xdydt,d2Xdzdt,d2Xdx2,d2Xdydx,d2Xdzdx,d2Xdxdy,d2Xdy2,d2Xdzdy,d2Xdxdz,d2Xdydz,d2Xdz2
 REAL(num), INTENT(IN) :: T
 REAL(num)	:: oneuponL=1./Lscl, oneuponT=1./Tscl
 REAL(num)	:: LuponT=Lscl/Tscl,LuponTV=Lscl/Tscl/Vscl, oneuponTL=1./Tscl/Lscl, LoB=Lscl/Bscl, ToB=Tscl/Bscl


 !!!!!!!!!!! Magnetic Dipole !!!!!!!!!!
 ! md(1)=0.
 ! md(2)=0.
 ! md(3)=-1.
 ! RMOD2= DOT_PRODUCT(R,R)
 ! mdr  = DOT_PRODUCT(md,R)
 ! RMOD5= RMOD2**(5./2.)
 !E=0.
 !B=(3.*mdr*R-md*RMOD2)/RMOD5
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!!!! Notice that X0 and Y0 and T represent dimensionless variables !!!!

 CALL SUB1_X0Y0Z0(R,T,X0,Y0,Z0,dX0dX,dX0dY,dX0dZ,dY0dX,dY0dY,dY0dZ,dZ0dX,dZ0dY,dZ0dZ,dX0dt,dY0dt,dZ0dt)

!CALL ABSETUP(R,T,X0,Y0,Z0,DX0DX,DX0DY,DX0DZ,DX0DT,DY0DX,DY0DY,DY0DZ,DY0DT,DZ0DX,DZ0DY,DZ0DZ,DZ0DT,&
!           D2X0DX2,D2X0DXDY,D2X0DXDZ,D2X0DXDT,D2X0DY2,D2X0DYDZ,D2X0DYDT,D2X0DZ2,D2X0DZDT,D2X0DT2,&
!           D2Y0DX2,D2Y0DXDY,D2Y0DXDZ,D2Y0DXDT,D2Y0DY2,D2Y0DYDZ,D2Y0DYDT,D2Y0DZ2,D2Y0DZDT,D2Y0DT2,&
!           D2Z0DX2,D2Z0DXDY,D2Z0DXDZ,D2Z0DXDT,D2Z0DY2,D2Z0DYDZ,D2Z0DYDT,D2Z0DZ2,D2Z0DZDT,D2Z0DT2,&
!           D2X0DYDX,D2X0DZDX,D2X0DTDX,D2X0DZDY,D2X0DTDY,D2X0DTDZ,&
!           D2Y0DYDX,D2Y0DZDX,D2Y0DTDX,D2Y0DZDY,D2Y0DTDY,D2Y0DTDZ,&
!           D2Z0DYDX,D2Z0DZDX,D2Z0DTDX,D2Z0DZDY,D2Z0DTDY,D2Z0DTDZ)

! km: Find B0, final state of field 

CALL B0setup(X0,Y0,Z0,B0)
!!!!!CALL dA0(X0,Y0,dA0dX0,dA0dY0)
! km: Set up the vector X
 X(1)=X0
 X(2)=Y0
 X(3)=Z0

 dXdx(1)=dX0dx
 dXdx(2)=dY0dx
 dXdx(3)=dZ0dx

 dXdy(1)=dX0dy
 dXdy(2)=dY0dy
 dXdy(3)=dZ0dy

 dXdz(1)=dX0dz
 dXdz(2)=dY0dz
 dXdz(3)=dZ0dz

! CALL dA0(X0,Y0,dA0dX0, dA0dY0) 
 
! JT: B(1) and B(2) are equations 28 and 29 of Guiliani et al. (2005)
! B(1) = oneuponL *  ( dA0dX0 * dX0dY + dA0dY0 * dY0dY )
! B(2) = oneuponL * ( -(dA0dX0 * dX0dX + dA0dY0 * dY0dX) )
! B(3) =  dY0dY * Bzfinal					!??

! km: Calculate B field using (A.4)-(A.7) in Grady paper

 B(1)= DOT(CROSS(dXdy,dXdz),B0)*oneuponL
 B(2)= DOT(CROSS(dXdz,dXdx),B0)*oneuponL
 B(3)= DOT(CROSS(dXdx,dXdy),B0)*oneuponL + (dX0dx*DY0dY-dX0dY*dY0dX)*Bzfinal

! km: 2D set up here

!!B(1) = oneuponL * (dA0dX0 * dX0dY + dA0dY0 * dY0dY)
!!B(2) = oneuponL * (-(dA0dX0 * dX0dX + dA0dY0 * dY0dX) )
!!B(3) = 0.0 

 B=B/Bscl  !!!  B is made dimensionless

 !!!!! Velocity field is also in dimensionless units !!!!!
 !!!!! Notice the factors L/T/Vscl

 DETERMINANT = dX0dx*(dY0dy*dZ0dz-dZ0dy*dY0dz)-dX0dy*(dY0dx*dZ0dz-dZ0dx*dY0dz)+&
               dX0dz*(dY0dx*dZ0dy-dZ0dx*dY0dy)

 Vf(1)= (-LuponTV)*(dX0dt*(dY0dy*dZ0dz-dZ0dy*dY0dz) +   &
                 dY0dt*(dX0dz*dZ0dy-dX0dy*dZ0dz) +   &
                 dZ0dt*(dX0dy*dY0dz-dX0dz*dY0dy) )/DETERMINANT
 
 Vf(2)= (-LuponTV)*(dX0dt*(dZ0dx*dY0dz-dZ0dz*dY0dx) +   &
                 dY0dt*(dX0dx*dZ0dz-dX0dz*dZ0dx) +   &
                 dZ0dt*(dX0dz*dY0dx-dX0dx*dY0dz) )/DETERMINANT

 Vf(3)= (-LuponTV)*(dX0dt*(dY0dx*dZ0dy-dY0dy*dZ0dx) +   &
                 dY0dt*(dX0dy*dZ0dx-dX0dx*dZ0dy) +   &
                 dZ0dt*(dX0dx*dY0dy-dX0dy*dY0dx) )/DETERMINANT

 E = -CROSS(Vf,B)   ! This electric field is dimensionless
 ! print*,'E(3)=',E(3),'****'
 !!!! Notice that the dimensionless electric field E(3) can also be calculated 
 !!!! as
 ! E(3) = -(1.0/Tscl)*(dA0dX0*dX0dt + dA0dY0*dY0dt)/(Vscl*B0)
 !!!! see notebook for more information !!!!!!!!!!!!!!!!
 !print*,'E(3)=',E(3)

 !!!!!!! Derivatives !!!!!!!!!!!!!!!!!!

 CALL SUB2_X0Y0Z0(R,T,Y0,dY0dX,dY0dY,dY0dZ,d2X0dXdY,d2X0dYdX,d2X0dXdZ,d2X0dZdX,d2X0dYdZ,&
                  d2X0dZdY,d2X0dX2,d2X0dY2,d2X0dZ2,d2Y0dXdY,d2Y0dYdX,d2Y0dXdZ,&
                  d2Y0dZdX,d2Y0dYdZ,d2Y0dZdY,d2Y0dX2,d2Y0dY2,d2Y0dZ2,d2Z0dXdY,&
                  d2Z0dYdX,d2Z0dXdZ,d2Z0dZdX,d2Z0dYdZ,d2Z0dZdY,d2Z0dX2,&
                  d2Z0dY2,d2Z0dZ2)
 CALL SUB3_X0Y0Z0(R,T,dY0dt,d2X0dtdX,d2X0dXdt,d2X0dtdY,d2X0dYdt,d2X0dtdZ,&
                  d2X0dZdt,d2Y0dtdX,d2Y0dXdt,d2Y0dtdY,d2Y0dYdt,d2Y0dtdZ,&
                  d2Y0dZdt,d2Z0dtdX,d2Z0dXdt,d2Z0dtdY,d2Z0dYdt,d2Z0dtdZ,&
                  d2Z0dZdt,d2X0dt2,d2Y0dt2,d2Z0dt2 )
!! CALL ddA0(X0,Y0,d2A0dx02,d2A0dy02,d2A0dx0dy0,d2A0dy0dx0)

 ! Now calculate derivative of Velocity field w.r.t. time

 !!! derivative of DETERMINANT is done with respect to normalized time !!!
 der_det=d2X0dtdX*(dY0dY*dZ0dZ-dZ0dY*dY0dz) +  &
         dX0dX*(d2Y0dYdt*dZ0dZ+dY0dY*d2Z0dZdt-d2Z0dYdt*dY0dZ-dZ0dY*d2Y0dZdt)- &
         d2X0dtdY*(dY0dX*dZ0dZ-dZ0dX*dY0dZ) -  &
         dX0dY*(d2Y0dXdt*dZ0dZ+dY0dX*d2Z0dZdt-d2Z0dXdt*dY0dZ-dZ0dX*d2Y0dZdt)+ &
         d2X0dtdZ*(dY0dX*dZ0dy-dZ0dX*dY0dY) +  &
         dX0dZ*(d2Y0dXdt*dZ0dY+dY0dX*d2Z0dYdt-d2Z0dXdt*dY0dY-dZ0dX*d2Y0dYdt)

 xder_det=d2X0dX2*(dY0dY*dZ0dZ-dZ0dY*dY0dZ) +  &
          dX0dX*(d2Y0dYdX*dZ0dZ+dY0dY*d2Z0dZdX-d2Z0dYdX*dY0dZ-dZ0dY*d2Y0dZdX)-&
          d2X0dYdX*(dY0dX*dZ0dZ-dZ0dX*dY0dZ) -  &
          dX0dY*(d2Y0dX2*dZ0dZ+dY0dX*d2Z0dZdX-d2Z0dX2*dY0dZ-dZ0dX*d2Y0dZdX)+ &
          d2X0dZdX*(dY0dX*dZ0dY-dZ0dX*dY0dY) +  &
          dX0dZ*(d2Y0dX2*dZ0dY+dY0dX*d2Z0dYdX-d2Z0dX2*dY0dY-dZ0dX*d2Y0dYdX)

 yder_det=d2X0dXdY*(dY0dY*dZ0dZ-dZ0dY*dY0dZ) + &
          dX0dX*(d2Y0dY2*dZ0dZ+dY0dY*d2Z0dZdY-d2Z0dY2*dY0dZ-dZ0dY*d2Y0dZdY)- &
          d2X0dY2*(dY0dX*dZ0dZ-dZ0dX*dY0dz) -  &
          dX0dY*(d2Y0dXdY*dZ0dZ+dY0dX*d2Z0dZdY-d2Z0dXdY*dY0dZ-dZ0dX*d2Y0dZdY)+&
          d2X0dZdY*(dY0dX*dZ0dY-dZ0dx*dY0dY) +  &
          dX0dZ*(d2Y0dXdY*dZ0dY+dY0dX*d2Z0dY2-d2Z0dXdY*dY0dY-dZ0dX*d2Y0dY2)

 zder_det=d2X0dXdZ*(dY0dY*dZ0dZ-dZ0dY*dY0dZ) +  &
          dX0dX*(d2Y0dYdZ*dZ0dZ+dY0dY*d2Z0dZ2-d2Z0dYdZ*dY0dZ-dZ0dY*d2Y0dZ2)- &
          d2X0dYdZ*(dY0dX*dZ0dZ-dZ0dX*dY0dZ) -  &
          dX0dY*(d2Y0dXdZ*dZ0dZ+dY0dX*d2Z0dZ2-d2Z0dXdZ*dY0dZ-dZ0dX*d2Y0dZ2)+ &
          d2X0dZ2*(dY0dX*dZ0dY-dZ0dX*dY0dY) +  &
          dX0dZ*(d2Y0dXdZ*dZ0dY+dY0dX*d2Z0dYdZ-d2Z0dXdZ*dY0dY-dZ0dX*d2Y0dYdZ)

 !!! Time derivatives of dmensionless Vf with respect to dimensionless time!!
 dVfdt(1)=(-LuponTV)*(d2X0dt2*(dY0dY*dZ0dZ-dZ0dY*dY0dZ)/DETERMINANT+ &
          dX0dt*(d2Y0dYdt*dZ0dZ+dY0dY*d2Z0dZdt-d2Z0dYdt*dY0dZ-dZ0dY*d2Y0dZdt)/DETERMINANT + &
          d2Y0dt2*(dX0dZ*dZ0dY-dX0dY*dZ0dZ)/DETERMINANT + &
          dY0dt*(d2X0dZdt*dZ0dY+dX0dZ*d2Z0dYdt-d2X0dYdt*dZ0dZ-dX0dY*d2Z0dZdt)/DETERMINANT + &
          d2Z0dt2*(dX0dY*dY0dZ-dX0dZ*dY0dY)/DETERMINANT + &
          dZ0dt*(d2X0dYdt*dY0dZ+dX0dY*d2Y0dZdt-d2X0dZdt*dY0dY-dX0dZ*d2Y0dYdt)/DETERMINANT + &
          (dX0dt*(dY0dY*dZ0dZ-dY0dZ*dZ0dY)+dY0dt*(dX0dZ*dZ0dY-dX0dY*dZ0dZ)+&
           dZ0dt*(dX0dY*dY0dZ-dX0dZ*dY0dY))*(-der_det/DETERMINANT/DETERMINANT))!*LuponT

 dVfdt(2)=(-LuponTV)*(d2X0dt2*(dZ0dX*dY0dZ-dY0dX*dZ0dZ)/DETERMINANT+ &
          dX0dt*(d2Z0dXdt*dY0dZ+dZ0dX*d2Y0dZdt-d2Y0dXdt*dZ0dZ-dY0dX*d2Z0dZdt)/DETERMINANT+ &
          d2Y0dt2*(dX0dX*dZ0dZ-dX0dZ*dZ0dX)/DETERMINANT+ &
          dY0dt*(d2X0dXdt*dZ0dZ+dX0dX*d2Z0dZdt-d2X0dZdt*dZ0dX-dX0dZ*d2Z0dXdt)/DETERMINANT+ &
          d2Z0dt2*(dX0dZ*dY0dX-dX0dX*dY0dZ)/DETERMINANT+ &
          dZ0dt*(d2X0dZdt*dY0dX+dX0dZ*d2Y0dXdt-d2X0dXdt*dY0dZ-dX0dX*d2Y0dZdt)/DETERMINANT+ &
          (dX0dt*(dZ0dX*dY0dZ-dY0dX*dZ0dZ)+dY0dt*(dX0dX*dZ0dZ-dX0dZ*dZ0dX)+&
           dZ0dt*(dX0dZ*dY0dX-dX0dX*dY0dZ))*(-der_det/DETERMINANT/DETERMINANT))!*LuponT

 dVfdt(3)=(-LuponTV)*(d2X0dt2*(dY0dX*dZ0dY-dZ0dX*dY0dY)/DETERMINANT+ &
          dX0dt*(d2Y0dXdt*dZ0dY+dY0dX*d2Z0dYdt-d2Z0dXdt*dY0dY-dZ0dX*d2Y0dYdt)/DETERMINANT+ &
          d2Y0dt2*(dX0dY*dZ0dX-dX0dX*dZ0dY)/DETERMINANT+ &
          dY0dt*(d2X0dYdt*dZ0dX+dX0dY*d2Z0dXdt-d2X0dXdt*dZ0dY-dX0dX*d2Z0dYdt)/DETERMINANT+ &
          d2Z0dt2*(dX0dX*dY0dY-dX0dY*dY0dX)/DETERMINANT+ &
          dZ0dt*(d2X0dXdt*dY0dY+dX0dX*d2Y0dYdt-d2X0dYdt*dY0dX-dX0dY*d2Y0dXdt)/DETERMINANT+ &
          (dX0dt*(dY0dX*dZ0dY-dZ0dX*dY0dY)+dY0dt*(dX0dY*dZ0dX-dX0dX*dZ0dY)+&
           dZ0dt*(dX0dX*dY0dY-dX0dY*dY0dX))*(-der_det/DETERMINANT/DETERMINANT))!*LuponT
 
!---- SPACE DERIVATIVES: dimensionless quantities
 dVfdx(1)=(-LuponTV)*(d2X0dXdt*(dY0dY*dZ0dZ-dZ0dY*dY0dZ)/DETERMINANT+ &
          dX0dt*(d2Y0dYdX*dZ0dZ+dY0dY*d2Z0dZdX-d2Z0dYdX*dY0dZ-dZ0dY*d2Y0dZdX)/DETERMINANT+ &
          d2Y0dtdX*(dX0dz*dZ0dY-dX0dY*dZ0dZ)/DETERMINANT+ &
          dY0dt*(d2X0dZdX*dZ0dY+dX0dZ*d2Z0dYdX-d2X0dYdX*dZ0dZ-dX0dY*d2Z0dZdX)/DETERMINANT+ &
          d2Z0dtdX*(dX0dY*dY0dZ-dX0dZ*dY0dY)/DETERMINANT+ &
          dZ0dt*(d2X0dYdX*dY0dZ+dX0dY*d2Y0dZdX-d2X0dZdX*dY0dY-dX0dZ*d2Y0dYdX)/DETERMINANT+ &
          (dX0dt*(dY0dY*dZ0dZ-dZ0dY*dY0dZ)+dY0dt*(dX0dZ*dZ0dY-dX0dY*dZ0dZ)+&
           dZ0dt*(dX0dY*dY0dZ-dX0dZ*dY0dY))*(-xder_det/DETERMINANT/DETERMINANT))

 dVfdx(2)=(-LuponTV)*(d2X0dtdX*(dZ0dX*dY0dZ-dY0dX*dZ0dZ)/DETERMINANT+ &
          dX0dt*(d2Z0dX2*dY0dZ+dZ0dX*d2Y0dZdX-d2Y0dX2*dZ0dZ-dY0dX*d2Z0dZdX)/DETERMINANT+ &
          d2Y0dtdX*(dX0dX*dZ0dZ-dX0dZ*dZ0dX)/DETERMINANT+ &
          dY0dt*(d2X0dX2*dZ0dZ+dX0dX*d2Z0dZdX-d2X0dZdX*dZ0dX-dX0dZ*d2Z0dX2)/DETERMINANT+ &
          d2Z0dtdX*(dX0dZ*dY0dX-dX0dX*dY0dZ)/DETERMINANT+ &
          dZ0dt*(d2X0dZdX*dY0dX+dX0dZ*d2Y0dX2-d2X0dX2*dY0dZ-dX0dX*d2Y0dZdX)/DETERMINANT+ &
          (dX0dt*(dZ0dX*dY0dZ-dY0dX*dZ0dZ)+dY0dt*(dX0dX*dZ0dZ-dX0dZ*dZ0dX)+&
           dZ0dt*(dX0dZ*dY0dX-dX0dX*dY0dZ))*(-xder_det/DETERMINANT/DETERMINANT))

 dVfdx(3)=(-LuponTV)*(d2X0dtdX*(dY0dX*dZ0dY-dZ0dX*dY0dY)/DETERMINANT+ &
          dX0dt*(d2Y0dX2*dZ0dY+dY0dX*d2Z0dYdX-d2Z0dX2*dY0dY-dZ0dX*d2Y0dYdX)/DETERMINANT+ &
          d2Y0dtdX*(dX0dY*dZ0dX-dX0dX*dZ0dY)/DETERMINANT+ &
          dY0dt*(d2X0dYdX*dZ0dX+dX0dY*d2Z0dX2-d2X0dX2*dZ0dY-dX0dX*d2Z0dYdX)/DETERMINANT+ &
          d2Z0dtdX*(dX0dX*dY0dY-dX0dY*dY0dX)/DETERMINANT+ &
          dZ0dt*(d2X0dX2*dY0dY+dX0dX*d2Y0dYdX-d2X0dYdX*dY0dX-dX0dY*d2Y0dX2)/DETERMINANT+ &
          (dX0dt*(dY0dX*dZ0dY-dZ0dX*dY0dY)+dY0dt*(dX0dY*dZ0dX-dX0dX*dZ0dY)+&
           dZ0dt*(dX0dX*dY0dY-dX0dY*dY0dX))*(-xder_det/DETERMINANT/DETERMINANT))
!----
 dVfdy(1)=(-LuponTV)*(d2X0dtdY*(dY0dY*dZ0dZ-dZ0dY*dY0dZ)/DETERMINANT+ &
          dX0dt*(d2Y0dY2*dZ0dZ*dY0dY*d2Z0dZdY-d2Z0dY2*dY0dZ-dZ0dY*d2Y0dZdY)/DETERMINANT+ &
          d2Y0dtdY*(dX0dZ*dZ0dY-dX0dY*dZ0dZ)/DETERMINANT+ &
          dY0dt*(d2X0dZdY*dZ0dY+dX0dZ*d2Z0dY2-d2X0dY2*dZ0dZ-dX0dY*d2Z0dZdY)/DETERMINANT+ &
          d2Z0dtdY*(dX0dY*dY0dZ-dX0dZ*dY0dY)/DETERMINANT+ &
          dZ0dt*(d2X0dY2*dY0dZ+dX0dY*d2Y0dZdY-d2X0dZdY*dY0dY-dX0dZ*d2Y0dY2)/DETERMINANT+ &
          (dX0dt*(dY0dY*dZ0dZ-dZ0dY*dY0dZ)+dY0dt*(dX0dZ*dZ0dY-dX0dY*dZ0dZ)+&
           dZ0dt*(dX0dY*dY0dZ-dX0dZ*dY0dY))*(-yder_det/DETERMINANT/DETERMINANT))

 dVfdy(2)=(-LuponTV)*(d2X0dtdY*(dZ0dX*dY0dZ-dY0dX*dZ0dZ)/DETERMINANT + &
          dX0dt*(d2Z0dXdY*dY0dZ+dZ0dX*d2Y0dZdY-d2Y0dXdY*dZ0dZ-dY0dX*d2Z0dZdY)/DETERMINANT+ &
          d2Y0dtdY*(dX0dX*dZ0dZ-dX0dZ*dZ0dX)/DETERMINANT+ &
          dY0dt*(d2X0dXdY*dZ0dZ+dX0dX*d2Z0dZdY-d2X0dZdY*dZ0dX-dX0dZ*d2Z0dXdY)/DETERMINANT+ &
          d2Z0dtdY*(dX0dZ*dY0dX-dX0dX*dY0dZ)/DETERMINANT+ &
          dZ0dt*(d2X0dZdY*dY0dX+dX0dZ*d2Y0dXdY-d2X0dXdY*dY0dZ-dX0dX*d2Y0dZdY)/DETERMINANT+ &
          (dX0dt*(dZ0dX*dY0dZ-dY0dX*dZ0dZ)+dY0dt*(dX0dX*dZ0dZ-dX0dZ*dZ0dX)+&
           dZ0dt*(dX0dZ*dY0dX-dX0dX*dY0dZ))*(-yder_det/DETERMINANT/DETERMINANT))

 dVfdy(3)=(-LuponTV)*(d2X0dtdY*(dY0dX*dZ0dY-dZ0dX*dY0dY)/DETERMINANT + &
          dX0dt*(d2Y0dXdY*dZ0dY+dY0dX*d2Z0dY2-d2Z0dXdY*dY0dY-dZ0dX*d2Y0dY2)/DETERMINANT+ &
          d2Y0dtdY*(dX0dY*dZ0dX-dX0dX*dZ0dY)/DETERMINANT+ &
          dY0dt*(d2X0dY2*dZ0dX+dX0dY*d2Z0dXdY-d2X0dXdY*dZ0dY-dX0dX*d2Z0dY2)/DETERMINANT+ &
          d2Z0dtdY*(dX0dX*dY0dY-dX0dY*dY0dX)/DETERMINANT+ &
          dZ0dt*(d2X0dXdY*dY0dY+dX0dX*d2Y0dY2-d2X0dY2*dY0dX-dX0dY*d2Y0dXdY)/DETERMINANT+ &
          (dX0dt*(dY0dX*dZ0dY-dZ0dX*dY0dY)+dY0dt*(dX0dY*dZ0dX-dX0dX*dZ0dY)+&
           dZ0dt*(dX0dX*dY0dY-dX0dY*dY0dX))*(-yder_det/DETERMINANT/DETERMINANT))
!----
 dVfdz(1)=(-LuponTV)*(d2X0dtdZ*(dY0dY*dZ0dZ-dZ0dY*dY0dZ)/DETERMINANT+ &
          dX0dt*(d2Y0dYdZ*dZ0dZ+dY0dY*d2Z0dZ2-d2Z0dYdZ*dY0dZ-dZ0dY*d2Y0dZ2)/DETERMINANT+ &
          d2Y0dtdZ*(dX0dZ*dZ0dY-dX0dY*dZ0dZ)/DETERMINANT+ &
          dY0dt*(d2X0dZ2*dZ0dY+dX0dZ*d2Z0dYdZ-d2X0dYdZ*dZ0dZ-dX0dY*d2Z0dZ2)/DETERMINANT+ &
          d2Z0dtdZ*(dX0dY*dY0dZ-dX0dZ*dY0dY)/DETERMINANT+ &
          dZ0dt*(d2X0dYdZ*dY0dZ+dX0dY*d2Y0dZ2-d2X0dZ2*dY0dY-dX0dZ*d2Y0dYdZ)/DETERMINANT+ &
          (dX0dt*(dY0dY*dZ0dZ-dZ0dY*dY0dZ)+dY0dt*(dX0dZ*dZ0dY-dX0dY*dZ0dZ)+&
           dZ0dt*(dX0dY*dY0dZ-dX0dZ*dY0dY))*(-zder_det/DETERMINANT/DETERMINANT))

 dVfdz(2)=(-LuponTV)*(d2X0dtdZ*(dZ0dX*dY0dZ-dY0dX*dZ0dZ)/DETERMINANT+ &
          dX0dt*(d2Z0dXdZ*dY0dZ+dZ0dX*d2Y0dZ2-d2Y0dXdZ*dZ0dZ-dY0dX*d2Z0dZ2)/DETERMINANT+ &
          d2Y0dtdZ*(dX0dX*dZ0dZ-dX0dZ*dZ0dX)/DETERMINANT+ &
          dY0dt*(d2X0dXdZ*dZ0dZ+dX0dX*d2Z0dZ2-d2X0dZ2*dZ0dX-dX0dZ*d2Z0dXdZ)/DETERMINANT+ &
          d2Z0dtdZ*(dX0dZ*dY0dX-dX0dX*dY0dZ)/DETERMINANT+ &
          dZ0dt*(d2X0dZ2*dY0dX+dX0dZ*d2Y0dXdZ-d2X0dXdZ*dY0dZ-dX0dX*d2Y0dZ2)/DETERMINANT+ &
          (dX0dt*(dZ0dX*dY0dZ-dY0dX*dZ0dZ)+dY0dt*(dX0dX*dZ0dZ-dX0dZ*dZ0dX)+&
           dZ0dt*(dX0dZ*dY0dX-dX0dX*dY0dZ))*(-zder_det/DETERMINANT/DETERMINANT))
 
 dVfdz(3)=(-LuponTV)*(d2X0dtdZ*(dY0dX*dZ0dY-dZ0dX*dY0dY)/DETERMINANT+ &
          dX0dt*(d2Y0dXdZ*dZ0dY+dY0dX*d2Z0dYdZ-d2Z0dXdZ*dY0dY-dZ0dX*d2Y0dYdZ)/DETERMINANT+ &
          d2Y0dtdZ*(dX0dY*dZ0dX-dX0dX*dZ0dY)/DETERMINANT+ &
          dY0dt*(d2X0dYdZ*dZ0dX+dX0dY*d2Z0dXdZ-d2X0dXdZ*dZ0dY-dX0dX*d2Z0dYdZ)/DETERMINANT+ &
          d2Z0dtdZ*(dX0dX*dY0dY-dX0dY*dY0dX)/DETERMINANT+ &
          dZ0dt*(d2X0dXdZ*dY0dY+dX0dX*d2Y0dYdZ-d2X0dYdZ*dY0dX-dX0dY*d2Y0dXdZ)/DETERMINANT+ &
          (dX0dt*(dY0dX*dZ0dY-dZ0dX*dY0dY)+dY0dt*(dX0dY*dZ0dX-dX0dX*dZ0dY)+&
           dZ0dt*(dX0dX*dY0dY-dX0dY*dY0dX))*(-zder_det/DETERMINANT/DETERMINANT)) 

! Km Set up 2nd derivatives of X vector to make B calculations easier

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

!km do the same for B0 derivatives

  CALL dB0(X0,Y0,Z0,dB0dX0,dB0dY0,dB0dZ0)

!km: We can now calculate the derivatives in terms of our regular co ords x,y,z
  dB0dx = dX0dx*dB0dx0 + dY0dx*dB0dy0 + dZ0dx*dB0dz0
  dB0dy = dX0dy*dB0dx0 + dY0dy*dB0dy0 + dZ0dy*dB0dz0
  dB0dz = dX0dz*dB0dx0 + dY0dz*dB0dy0 + dZ0dz*dB0dz0
  dB0dt = dX0dt*dB0dx0 + dY0dt*dB0dy0 + dZ0dt*dB0dz0

!!!! Check: LHS and RHS should be equal !!!!!!
 !print*,'1-LHS=',(der_det+Vf(1)*xder_det+Vf(2)*yder_det+Vf(3)*zder_det)/DETERMINANT
 !print*,'2-RHS=',-(dVfdx(1)+dVfdy(2)+dVfdz(3))
 !print*,'****'

!!!! Magnetic field derivatives in dimensional form !!!
! Have to check which things are dimensional here, B0 is not
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
         DOT(CROSS(dXdx,dXdy),dB0dx)) + &
         oneuponL*(d2X0dX2*dY0dY+dX0dX*d2Y0dYdX-d2X0dYdX*dY0dX-dX0dY*d2Y0dX2)*Bzfinal   
 !dBzdy
 dBdy(3)=oneuponL*oneuponL*(DOT(CROSS(d2Xdxdy,dXdy)+CROSS(dXdx,d2Xdy2),B0)+&
         DOT(CROSS(dXdx,dXdy),dB0dy)) + &
         oneuponL*(d2X0dXdY*dY0dY+dX0dX*d2Y0dY2-d2X0dY2*dY0dX-dX0dY*d2Y0dXdY)*Bzfinal
 !dBzdz
 dBdz(3)=oneuponL*oneuponL*(DOT(CROSS(d2Xdxdz,dXdy)+CROSS(dXdx,d2Xdydz),B0)+&
         DOT(CROSS(dXdx,dXdy),dB0dz)) + &
         oneuponL*(d2X0dXdZ*dY0dY+dX0dX*d2Y0dYdZ-d2X0dYdZ*dY0dX-dX0dY*d2Y0dXdZ)*Bzfinal
 !dBxdt
 dBdt(1)=oneuponTL*(DOT(CROSS(d2Xdydt,dXdz)+CROSS(dXdy,d2Xdzdt),B0)+&
         DOT(CROSS(dXdy,dXdz),dB0dt))
 !dBydt
 dBdt(2)=oneuponTL*(DOT(CROSS(d2Xdzdt,dXdx)+CROSS(dXdz,d2Xdxdt),B0)+&
         DOT(CROSS(dXdz,dXdx),dB0dt))
 !dBzdt
 dBdt(3)=oneuponTL*(DOT(CROSS(d2Xdxdt,dXdy)+CROSS(dXdx,d2Xdydt),B0)+&
         DOT(CROSS(dXdx,dXdy),dB0dt)) + &
         oneuponL*(d2X0dXdt*dY0dY+dX0dX*d2Y0dYdt-d2X0dYdt*dY0dX-dX0dY*d2Y0dXdt)*Bzfinal

! km: 2D set up
!!dbdx(1)=oneuponl*oneuponl*((d2A0dX02*dX0dx+d2A0dY0dX0*dY0dx)*dX0dy+dA0dX0*d2X0dxdy + &
!!        (d2A0dX0dY0*dX0dx+d2A0dY02*dY0dx)*dY0dy+dA0dY0*d2Y0dxdy) 
!!dbdy(1)=oneuponl*oneuponl*((d2A0dX02*dX0dy+d2A0dY0dX0*dY0dy)*dX0dy+dA0dX0*d2X0dy2 + &
!!        (d2A0dX0dY0*dX0dy+d2A0dY02*dY0dy)*dY0dy+dA0dY0*d2Y0dy2)
!!dbdz(1)=0.0
!!dbdx(2)=-oneuponl*oneuponl*((d2A0dX02*dX0dx+d2A0dY0dX0*dY0dx)*dX0dx+dA0dX0*d2X0dx2 + &
!!        (d2A0dX0dY0*dX0dx+d2A0dY02*dY0dx)*dY0dx+dA0dY0*d2Y0dx2)
!!dbdy(2)=-oneuponl*oneuponl*((d2A0dX02*dX0dy+d2A0dY0dX0*dY0dy)*dX0dx+dA0dX0*d2X0dydx + &
!!        (d2A0dX0dY0*dX0dy+d2A0dY02*dY0dy)*dY0dx+dA0dY0*d2Y0dydx)
!!dbdz(2)=0.0
!!dbdx(3)=0.0
!!dbdy(3)=0.0
!!dbdz(3)=0.0
!!dbdt(1)=oneuponTL*(d2A0dX02*dX0dt*dX0dy+d2A0dY0dX0*dY0dt*dX0dy + &
!!                   d2A0dX0dY0*dX0dt*dY0dy+d2A0dY02*dY0dt*dY0dy + &
!!                   dA0dX0*d2X0dtdy+dA0dY0*d2Y0dtdy)
!!dbdt(2)=oneuponTL*(d2A0dX02*dX0dt*dX0dx+d2A0dY0dX0*dY0dt*dX0dx + &
!!                  d2A0dX0dY0*dX0dt*dY0dx+d2A0dY02*dY0dt*dY0dx + &
!!                   dA0dX0*d2X0dtdx+dA0dY0*d2Y0dtdx)
!!dbdt(3)=0.0
 !! Make derivatives dimensionless !!!!
 DBDX = LoB*DBDX 
 DBDY = LoB*DBDY
 DBDZ = LoB*DBDZ
 DBDT = ToB*DBDT 

 IF (T.eq.1.234) THEN
  PRINT *,DBDT
 ENDIF

 !!! The dimensionless dEdt is !!
 dEdt=-(cross(dVfdt,B)+cross(Vf,DBDT))
 !print*,'dEdt(3)=',dEdt(3),'***'

 !!! The dimensionless dEdx is !!
 dEdx=-(cross(dVfdx,B)+cross(Vf,DBDX)) 
! print*,'dEdx(3)=',dEdx(3)

 !!! The dimensionless dEdy is !!
 dEdy=-(cross(dVfdy,B)+cross(Vf,DBDY))
 !print*,'dEdy(3)=',dEdy(3)

 !!! The dimensionless dEdz is !!
 dEdz=-(cross(dVfdz,B)+cross(Vf,DBDZ))
 !print*,'dEdz(3)=',dEdz(3)


 !!!! Dimensionless derivatives of Ez can also be calculated as follows !!!
 !!!! (L/Vscl/B0) and (Tscl/Vscl/B0) reduce quantities to dimensionless.
 !!dEzdx
 !dEdx(3)=-(oneuponL/Tscl)*(L/Vscl/B0)* &
 !                 ((d2A0dx02*dx0dx+d2A0dy0dx0*dy0dx)*dx0dt+dA0dx0*d2x0dxdt+&
 !                 (d2A0dx0dy0*dx0dx+d2A0dy02*dy0dx)*dy0dt+dA0dy0*d2y0dxdt)
 !!dEzdy
 !dEdy(3)=-(oneuponL/Tscl)*(L/Vscl/B0)* &
 !                     ((d2A0dx02*dx0dy+d2A0dy0dx0*dy0dy)*dx0dt+dA0dx0*d2x0dydt+&
 !                     (d2A0dx0dy0*dx0dy+d2A0dy02*dy0dy)*dy0dt+dA0dy0*d2y0dydt)
 !!dEzdz
 !dEdz(3)=0.
 !
 !!dEdt
 !dEdt(3)=-oneuponT*oneuponT*(Tscl/Vscl/B0)* &
 !                    (d2A0dx02 * dx0dt*dx0dt + d2A0dy0dx0 *dy0dt*dx0dt + &
 !                     dA0dx0 * d2x0dt2 + &
 !                     d2A0dx0dy0 * dx0dt*dy0dt + d2A0dy02*dy0dt*dy0dt + &
 !                     dA0dy0 * d2y0dt2 )
 !print*,'DEDT(3)=',DEDT(3)
 !print*,'DEDX(3)=',DEDX(3)
 !print*,'DEDY(3)=',DEDY(3)
 !print*,'DEDZ(3)=',DEDZ(3)

 !print*,X0   ! OK
 !print*,Y0   ! OK
 !print*,dY0dY ! 0k
 !print*,dY0dX ! OK
 !print*,dY0dT ! OK
 !print*,dX0dX ! OK
 !print*,dX0dY ! OK
 !print*,dX0dt ! OK
 !print*,d2X0dTdX !OK
 !print*,d2X0dTdY !OK
 !print*,d2Y0dTdX !OK
 !print*,d2Y0dTdY !Ok
 !print*,d2X0dXdY !OK
 !print*,d2X0dYdX !OK
 !print*,d2X0dX2  !OK
 !print*,d2X0dY2  !OK
 !print*,d2X0dt2  !OK
 !print*,d2Y0dXdY  !Ok
 !print*,d2Y0dYdX  !OK
 !print*,d2Y0dX2   !OK
 !print*,d2Y0dY2   !OK
 !print*,d2Y0dt2   !OK 

 !print*,dBdy
 !print*,dBdt ! OK
 !print*,E    ! OK
 !print*,DEDT
 !print*,DEDX
 !print*,DEDY
 !print*,DEDZ

 !B(1)=0
 !B(2)=0
 !B(3)=B0*T*exp(-R(2))
 !DBDX=0; DBDY=0; DBDZ=0; 
 !DBDY(3)=(-oneuponL)*B(3)
 !DBDT(1)=0
 !DBDT(2)=0
 !DBDT(3)=B0*exp(-R(2))/Tscl
 !E=0
 !DEDX=0; DEDY=0; DEDZ=0; DEDT=0
 !Vf=0
 !dVfdx=0; dVfdy=0; dVfdz=0; dVfdt=0
 !
 !B=B/B0
 !DBDY=(L/B0)*DBDY
 !DBDT=(Tscl/B0)*DBDT

END SUBROUTINE CMTFIELDS
!-------------------------------------------------------------------------------------------!
SUBROUTINE SUB1_X0Y0Z0(R,T,X0,Y0,Z0,dX0dX,dX0dY,dX0dZ,dY0dX,dY0dY,dY0dZ,dZ0dX,dZ0dY,dZ0dZ,dX0dt,dY0dt,dZ0dt)
!JT+ reads in a given position vector at time t, 
! and calculates Lagrangian transformation, x0(=x_inf),y0(=y_inf),z0(=z_inf)
! based on eq. (36) of Guiliani et al (2005)

 REAL(num), DIMENSION(3), INTENT(IN)	:: R
 REAL(num), INTENT(IN)  		:: T
 REAL(num), INTENT(OUT)			:: X0, Y0,Z0, dX0dX, dX0dY,dX0dZ, dY0dX, dY0dY,dY0dZ,dZ0dX,dZ0dY,dZ0dZ, dX0dt, dY0dt,dZ0dt
 REAL(num)				:: a=0.9,b=0.9, delta=1.0,a3d=1.0
 REAL(num)                              :: ay=1.0, yh=0.0, oneuponden

 !!!!!!  please notice Lv/L in the expressions  !!!!!

 IF ( (1.+(R(2)/((cc*T)**esp))) .le.0.) THEN
  PRINT *, '(1.+(R(2)/((cc*T)**esp)))=',(1.+(R(2)/((cc*T)**esp)))
  PRINT *, 'R(2)=',R(2)
  PRINT *,'T=',T
 END IF

oneuponden=1.0/(a3d*a3d+R(1)*R(1)+ay*(R(2)-yh)*ay*(R(2)-yh)+R(3)*R(3))

 Y0=(cc*T)**esp*log (1.+(R(2)/((cc*T)**esp)))*  &	
      (1.+tanh((R(2)-Lv/Lscl)*a))*0.5 +           &
       (1.-tanh((R(2)-Lv/Lscl)*b))*0.5*R(2)

! X0=R(1)-delta*(Y0-R(2))*R(3)/(a3d*a3d+R(1)*R(1)+R(3)*R(3))
 X0=R(1)-delta*(Y0-R(2))*R(3)*oneuponden

 !Z0=R(3)+delta*(Y0-R(2))*R(1)/(a3d*a3d+R(1)*R(1)+R(3)*R(3))
 Z0=R(3)+delta*(Y0-R(2))*R(1)*oneuponden
 !print*,X0,Y0

!JT: are these boundary conditions?

 dY0dX=0.
 dY0dY=0.5*(1.+tanh((R(2)-Lv/Lscl)*a))/(1.+R(2)/(cc*T)**esp)+        &
      0.5*(cc*T)**esp*log(1.+R(2)/(cc*T)**esp)*               &
      (1.-tanh((R(2)-Lv/Lscl)*a)**2)*a-0.5*                         &
      (1.-tanh((R(2)-Lv/Lscl)*b)**2)*b*R(2)+0.5-0.5*tanh((R(2)-Lv/Lscl)*b) 
 dY0dZ=0.

! dX0dX=1.0+2.0*R(1)*delta*(Y0-R(2))*R(3)/(a3d*a3d+R(1)*R(1)+R(3)*R(3))**2
! dX0dY=delta*(1.0-dY0dY)*R(3)/(a3d*a3d+R(1)*R(1)+R(3)*R(3))
! dX0dZ=(-1.0)*delta*(Y0-R(2))*(a3d*a3d+R(1)*R(1)+R(3)*R(3)-2.0*R(3)*R(3))/  &
!       (a3d*a3d+R(1)*R(1)+R(3)*R(3))**2

 
 dX0dX=1.0+2.0*R(1)*delta*(Y0-R(2))*R(3)*oneuponden*oneuponden - &
       delta*dY0dX*R(3)*oneuponden
 dX0dY=delta*(1.0-dY0dY)*R(3)*oneuponden+ &
       2.0*delta*(Y0-R(2))*ay*ay*(R(2)-yh)*R(3)*oneuponden*oneuponden
 dX0dZ=(-1.0)*delta*(Y0-R(2))*(a3d*a3d+R(1)*R(1)+ay*(R(2)-yh)*ay*(R(2)-yh)-R(3)*R(3))* &
                             oneuponden*oneuponden - &
       delta*dY0dZ*R(3)*oneuponden

! dZ0dX=delta*(Y0-R(2))*(a3d*a3d+R(1)*R(1)+R(3)*R(3)-2*R(1)*R(1))/   &
!       (a3d*a3d+R(1)*R(1)+R(3)*R(3))**2
! dZ0dY=delta*(dY0dY-1.0)*R(1)/(a3d*a3d+R(1)*R(1)+R(3)*R(3))
! dZ0dZ=1.0-2*delta*(Y0-R(2))*R(1)*R(3)/(a3d*a3d+R(1)*R(1)+R(3)*R(3))**2

 dZ0dX=delta*(Y0-R(2))*(a3d*a3d+ay*(R(2)-yh)*ay*(R(2)-yh)+R(3)*R(3)-R(1)*R(1))* &
                             oneuponden*oneuponden + &
       delta*dY0dX*R(1)*oneuponden
 dZ0dY=delta*(dY0dY-1.0)*R(1)*oneuponden-2.0*delta*(Y0-R(2))*ay*ay*(R(2)-yh)*R(1)*oneuponden*oneuponden
 dZ0dZ=1.0-2*delta*(Y0-R(2))*R(1)*R(3)*oneuponden*oneuponden + &
       delta*dY0dZ*R(1)*oneuponden

! dY0dt=0.5*(cc*T)**esp*esp*log(1.+R(2)/(cc*T)**esp)*          &
!      (1.+tanh((R(2)-Lv/Lscl)*a))/T-0.5*R(2)*esp*                  &
!      (1.+tanh((R(2)-Lv/Lscl)*a))/(T*(1.+R(2)/(cc*T)**esp)) 
! dX0dt=(-1.0)*delta*dY0dt*R(3)/(a3d*a3d+R(1)*R(1)+R(3)*R(3))
! dZ0dt=delta*dY0dt*R(1)/(a3d*a3d+R(1)*R(1)+R(3)*R(3))         

 dY0dt=0.5*(cc*T)**esp*esp*log(1.+R(2)/(cc*T)**esp)*          &
      (1.+tanh((R(2)-Lv/Lscl)*a))/T-0.5*R(2)*esp*                  &
      (1.+tanh((R(2)-Lv/Lscl)*a))/(T*(1.+R(2)/(cc*T)**esp)) 
 dX0dt=(-1.0)*delta*dY0dt*R(3)*oneuponden
 dZ0dt=delta*dY0dt*R(1)*oneuponden         

 !**********************************************************
 ! X0=R(1)
 ! Y0=R(2) 

 ! dX0dX=1.
 ! dX0dY=0.
 ! dY0dX=0.
 ! dY0dY=1.

 !dX0dt=0.
 !dY0dt=0.
     
END SUBROUTINE SUB1_X0Y0Z0
!-------------------------------------------------------------------------------------------!
SUBROUTINE SUB2_X0Y0Z0(R,T,Y0,dY0dX,dY0dY,dY0dZ,d2X0dXdY,d2X0dYdX,d2X0dXdZ,d2X0dZdX,d2X0dYdZ,&
                       d2X0dZdY,d2X0dX2,d2X0dY2,d2X0dZ2,d2Y0dXdY,d2Y0dYdX, &
                       d2Y0dXdZ,d2Y0dZdX,d2Y0dYdZ,d2Y0dZdY,d2Y0dX2,d2Y0dY2,&
                       d2Y0dZ2,d2Z0dXdY,d2Z0dYdX,d2Z0dXdZ,d2Z0dZdX,d2Z0dYdZ,&
                       d2Z0dZdY,d2Z0dX2,d2Z0dY2,d2Z0dZ2)
 
 REAL(num), DIMENSION(3), INTENT(IN) :: R
 REAL(num), INTENT(IN)  	:: T, dY0dX, dY0dY, dY0dZ !included to save copying long statements from SUB1_X0Y0Z0
 REAL(num), INTENT(OUT) 	:: d2X0dXdY,d2X0dYdX,d2X0dXdZ,d2X0dZdX,d2X0dYdZ,d2X0dZdY,d2X0dX2,d2X0dY2,d2X0dZ2
 REAL(num), INTENT(OUT) 	:: d2Y0dXdY,d2Y0dYdX,d2Y0dXdZ,d2Y0dZdX,d2Y0dYdZ,d2Y0dZdY,d2Y0dX2,d2Y0dY2,d2Y0dZ2
 REAL(num), INTENT(OUT)         :: d2Z0dXdY,d2Z0dYdX,d2Z0dXdZ,d2Z0dZdX,d2Z0dYdZ,d2Z0dZdY,d2Z0dX2,d2Z0dY2,D2Z0dZ2
 REAL(num)              	:: y, den, Y0
 REAL(num)                      :: ay=1.0, yh=0.0, oneuponden
 REAL(num)              	:: a=0.9,b=0.9,delta=1.0,a3d=1.0

 y=R(2)
 !den = a3d*a3d+R(1)*R(1)+R(3)*R(3) !km: denominator from (37) and (38) of Grady(2009)
den = a3d*a3d+R(1)*R(1)+ay*(R(2)-yh)*ay*(R(2)-yh)+R(3)*R(3)
oneuponden=1.0/(a3d*a3d+R(1)*R(1)+ay*(R(2)-yh)*ay*(R(2)-yh)+R(3)*R(3))
 !!!!  please notice Lv/L in the expressions !!!!
 
 d2y0dy2= -0.1e1 / (0.1e1 + y / (cc * t) ** esp) ** 2 * (0.1e1 + tan&
     &h((y - Lv/Lscl) * a)) / (cc * t) ** esp / 0.2e1 + 0.1e1 / (0.1e1 + y /&
     & (cc * t) ** esp) * (0.1e1 - tanh((y - Lv/Lscl) * a) ** 2) * a - (cc& 
     &* t) ** esp * log(0.1e1 + y / (cc * t) ** esp) * tanh((y - Lv/Lscl) * a&
     &) * (0.1e1 - tanh((y - Lv/Lscl) * a) ** 2) * a*a + tanh((y - Lv/Lscl)& 
     &* b) * (0.1e1 - tanh((y - Lv/Lscl) * b) ** 2) * b*b * y - (0.1e1& 
     &- tanh((y - Lv/Lscl) * b) ** 2) * b
 d2Y0dXdY=0.
 d2Y0dYdX=0.
 d2Y0dXdZ=0.
 d2Y0dZdX=0.
 d2Y0dYdZ=0.
 d2Y0dZdY=0.
 d2Y0dX2=0.
 d2Y0dZ2=0.

! d2X0dXdY = 2*delta*(dY0dY-1.0)*R(1)*R(3)/den/den
! d2X0dYdX = d2X0dXdY
! d2X0dXdZ = 2*delta*(Y0-R(2))*(R(1)*den-4*R(1)*R(3)*R(3))/den/den/den
! d2X0dZdX = d2X0dXdZ
! d2X0dYdZ = (-delta)*(dY0dY-1.0)*(den-2*R(3)*R(3))/den/den
! d2X0dZdY = d2X0dYdZ
! d2X0dX2 = 2*delta*(Y0-R(2))*(R(3)*den-4*R(1)*R(1)*R(3))/den/den/den
! d2X0dY2 = (-delta)*d2Y0dY2*R(3)/den
! d2X0dZ2 = delta*(Y0-R(2))*(6*R(3)*den-8*R(3)*R(3)*R(3))/den/den/den


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

! d2Z0dXdY = delta*(dY0dY-1.0)*(den-2*R(1)*R(1))/den/den
! d2Z0dYdX = d2Z0dXdY
! d2Z0dXdZ = delta*(Y0-R(2))*(8*R(1)*R(1)*R(3)-2*R(3)*den)/den/den/den
! d2Z0dZdX = d2Z0dXdZ
! d2Z0dYdZ = (-delta)*(dY0dY-1.0)*2*R(3)*R(1)/den/den
! d2Z0dZdY = d2Z0dYdZ
! d2Z0dX2 = delta*(Y0-R(2))*(8*R(1)*R(1)*R(1)-6*R(1)*den)/den/den/den
! d2Z0dY2 = delta*d2Y0dY2*R(1)/den
! d2Z0dZ2 = (-delta)*(Y0-R(2))*(2*R(1)*den-8*R(1)*R(3)*R(3))/den/den/den

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

!*********************************************************
!d2y0dy2=0.

END SUBROUTINE SUB2_X0Y0Z0
!-------------------------------------------------------------------------------------------!
SUBROUTINE SUB3_X0Y0Z0(R,T,dY0dt,d2X0dtdX,d2X0dXdt,d2X0dtdY,d2X0dYdt,d2X0dtdZ,d2X0dZdt, &
                         d2Y0dtdX,d2Y0dXdt,d2Y0dtdY,d2Y0dYdt,d2Y0dtdZ,d2Y0dZdt, &
                         d2Z0dtdX,d2Z0dXdt,d2Z0dtdY,d2Z0dYdt,d2Z0dtdZ,d2Z0dZdt, &
                                           d2X0dt2,d2Y0dt2,d2Z0dt2 )
 REAL(num), DIMENSION(3), INTENT(IN)	:: R
 REAL(num), INTENT(IN)  		:: T, dY0dt
 REAL(num), INTENT(OUT) 		:: d2X0dtdX,d2X0dXdt,d2X0dtdY,d2X0dYdt,d2X0dtdZ,d2X0dZdt
 REAL(num), INTENT(OUT) 		:: d2Y0dtdX,d2Y0dXdt,d2Y0dtdY,d2Y0dYdt,d2Y0dtdZ,d2Y0dZdt
 REAL(num), INTENT(OUT)                 :: d2Z0dtdX,d2Z0dXdt,d2Z0dtdY,d2Z0dYdt,d2Z0dtdZ,d2Z0dZdt
 REAL(num), INTENT(OUT) 		:: d2X0dt2, d2Y0dt2, d2Z0dt2
 REAL(num)             			:: y, den
 REAL(num)				:: a=0.9, delta=1.0, a3d=1.0
 REAL(num)                              :: ay=1.0, yh=0.0, oneuponden 

y=R(2)
!den=a3d*a3d+R(1)*R(1)+R(3)*R(3)
den=a3d*a3d+R(1)*R(1)+ay*(R(2)-yh)*ay*(R(2)-yh)+R(3)*R(3)
oneuponden=1.0/(a3d*a3d+R(1)*R(1)+ay*(R(2)-yh)*ay*(R(2)-yh)+R(3)*R(3))

!!!!! Lv/L in the expressions !!!!!

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

!d2X0dtdx = 2*delta*dY0dt*R(1)*R(3)/den/den
!d2X0dxdt = d2X0dtdx
!d2X0dydt = (-delta)*d2Y0dtdy*R(3)/den
!d2X0dtdy = d2X0dydt
!d2X0dtdz = (-delta)*dY0dt*(den-2*R(3)*R(3))/den/den
!d2X0dzdt = d2X0dtdz
!d2X0dt2 = (-delta)*d2y0dt2*R(3)/den

d2X0dtdx = 2*delta*dY0dt*R(1)*R(3)*oneuponden*oneuponden-delta*d2Y0dXdt*R(3)*oneuponden
d2X0dxdt = d2X0dtdx
d2X0dydt = (-delta)*d2Y0dtdy*R(3)*oneuponden+2*delta*dY0dt*ay*ay*(R(2)-yh)*R(3)*oneuponden*oneuponden
d2X0dtdy = d2X0dydt
d2X0dtdz = (-delta)*dY0dt*(den-2*R(3)*R(3))*oneuponden*oneuponden-delta*d2Y0dzdt*R(3)*oneuponden
d2X0dzdt = d2X0dtdz
d2X0dt2 = (-delta)*d2y0dt2*R(3)*oneuponden

!d2Z0dtdx = delta*dY0dt*(den-2*R(1)*R(1))/den/den
!d2Z0dxdt = d2Z0dtdx
!d2Z0dtdy = delta*d2Y0dydt*R(1)/den
!d2Z0dydt = d2Z0dtdy
!d2Z0dtdz = (-delta)*dY0dt*2*R(1)*R(3)/den/den
!d2Z0dzdt = d2Z0dtdz
!d2Z0dt2 = delta*d2Y0dt2*R(1)/den

d2Z0dtdx = delta*dY0dt*(den-2*R(1)*R(1))*oneuponden*oneuponden+delta*d2Y0dXdt*R(1)*oneuponden
d2Z0dxdt = d2Z0dtdx
d2Z0dtdy = delta*d2Y0dydt*R(1)*oneuponden-2*delta*dY0dt*R(1)*ay*ay*(R(2)-yh)*oneuponden*oneuponden
d2Z0dydt = d2Z0dtdy
d2Z0dtdz = (-delta)*dY0dt*2*R(1)*R(3)*oneuponden*oneuponden+delta*d2Y0dZdt*R(1)*oneuponden
d2Z0dzdt = d2Z0dtdz
d2Z0dt2 = delta*d2Y0dt2*R(1)*oneuponden


!******************************************************
!d2y0dtdy=0.
!d2y0dydt=0.
!d2y0dt2=0.

END SUBROUTINE SUB3_X0Y0Z0
!-------------------------------------------------------------------------------------------!
SUBROUTINE dA0 (X0,Y0,dA0dX0, dA0dY0)
!+ returns differentials dA0/dx0 and dA0/dy0 for specific values of x0 and y0
! based on equation 31 of Guiliani et al (2005) 
!JT edited - all numbers forced to be real, and powers unrolled
 
 REAL(num), INTENT(IN)  :: X0,Y0
 REAL(num), INTENT(OUT) :: dA0dX0, dA0dY0


!!!!!!!!!!! These are the derivatives with respect to
!!!!!!!!!!! the dimensional variables, not to be used here
! dA0dX0 = 32.*(Y0+d)*c1*X0*L/  &
!         ((4.*X0**2-4.*X0*L+L**2+4*Y0**2+8.*Y0*d+4*d**2)*   &
!         (4.*X0**2+4*X0*L+L**2+4*Y0**2+8*Y0*d+4*d**2)) 
! dA0dY0 = -4.*c1*L*(4.*X0**2-L**2-4*Y0**2-8*Y0*d-4*d**2)/  & 
!         ((4*X0**2-4*X0*L+L**2+4*Y0**2+8*Y0*d+4*d**2)*     &
!         (4*X0**2+4*X0*L+L**2+4*Y0**2+8*Y0*d+4*d**2)) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!dA0dX0 = 32.*c1*(Y0*L+d)*L**3*X0/  &
! ((4*L**2*X0**2-4*L**2*X0+L**2+4*Y0**2*L**2+8*Y0*L*d+4*d**2)* &
! (4*L**2*X0**2+4*L**2*X0+L**2+4*Y0**2*L**2+8*Y0*L*d+4*d**2))
!
!dA0dY0 = 4.*c1*L**2*  &
! (-4*L**2*X0**2+L**2+4*Y0**2*L**2+8*Y0*L*d+4*d**2)/  &
! ((4*L**2*X0**2-4*L**2*X0+L**2+4*Y0**2*L**2+8*Y0*L*d+4*d**2)*  &
! (4*L**2*X0**2+4*L**2*X0+L**2+4*Y0**2*L**2+8*Y0*L*d+4*d**2))

dA0dX0 = 32.*c1*(Y0*Lscl+d)*Lscl*Lscl*Lscl*X0/  &
 ((4.*Lscl*Lscl*X0*X0-4.*Lscl*Lscl*X0+Lscl*Lscl+4.*Y0*Y0*Lscl*Lscl+8.*Y0*Lscl*d+4.*d*d)* &
 (4.*Lscl*Lscl*X0*X0+4.*Lscl*Lscl*X0+Lscl*Lscl+4.*Y0*Y0*Lscl*Lscl+8.*Y0*Lscl*d+4.*d*d))

dA0dY0 = 4.*c1*Lscl*Lscl*  &
 (-4.*Lscl*Lscl*X0*X0+Lscl*Lscl+4.*Y0*Y0*Lscl*Lscl+8.*Y0*Lscl*d+4.*d*d)/  &
 ((4.*Lscl*Lscl*X0*X0-4.*Lscl*Lscl*X0+Lscl*Lscl+4*Y0*Y0*Lscl*Lscl+8.*Y0*Lscl*d+4.*d*d)*  &
 (4.*Lscl*Lscl*X0*X0+4.*Lscl*Lscl*X0+Lscl*Lscl+4.*Y0*Y0*Lscl*Lscl+8.*Y0*Lscl*d+4.*d*d))


END SUBROUTINE dA0
!-------------------------------------------------------------------------------------------!
SUBROUTINE ddA0 (X0,Y0,d2A0dx02, d2A0dy02, d2A0dx0dy0, d2A0dy0dx0 )
! second differentials of flux function A0,
! also based on Eq. 31 of Guiliani et al (2005)

!+ edit to make all numbers real and unroll powers 
 REAL(num), INTENT(OUT) :: d2A0dx02, d2A0dy02, d2A0dx0dy0, d2A0dy0dx0
 REAL(num), INTENT(IN)  :: X0,Y0

 d2A0dx02 = 32. * c1 * (y0 * Lscl + d) * Lscl*Lscl*Lscl / (4. * Lscl*Lscl * x0*x0 - & 
       4. * Lscl*Lscl * x0 + Lscl*Lscl + 4. * y0*y0 * Lscl*Lscl + 8. * y0 * Lscl * d + &
        4. * d*d) / (4. * Lscl*Lscl * x0*x0 + 4. * Lscl*Lscl * x0 + Lscl*Lscl +  &
       4. * y0*y0 * Lscl*Lscl + 8. * y0 * Lscl * d + 4. * d*d) - 32. * c1 * (y0 &
       * Lscl + d) * Lscl*Lscl*Lscl * x0 / (4. * Lscl*Lscl * x0*x0 - 4. * Lscl*Lscl * x0 &
        + Lscl*Lscl + 4. * y0*y0 * Lscl*Lscl + 8. * y0 * Lscl * d + 4. * d*d) ** &
        2 / (4. * Lscl*Lscl * x0*x0 + 4. * Lscl*Lscl * x0 + Lscl*Lscl + 4. * y0 **  & 
       2 * Lscl*Lscl + 8. * y0 * Lscl * d + 4. * d*d) * (8. * Lscl*Lscl * x0 - 4. * &
        Lscl*Lscl) - 32. * c1 * (y0 * Lscl + d) * Lscl*Lscl*Lscl * x0 / (4. * Lscl*Lscl * x0*x0 &
         - 4. * Lscl*Lscl * x0 + Lscl*Lscl + 4. * y0*y0 * Lscl*Lscl + 8. * y0 * &
        Lscl * d + 4. * d*d) / (4. * Lscl*Lscl * x0*x0 + 4. * Lscl*Lscl * x0 + Lscl*Lscl + &
        4. * y0*y0 * Lscl*Lscl + 8. * y0 * Lscl * d + 4. * d*d) ** 2 * &
       (8. * Lscl*Lscl * x0 + 4. * Lscl*Lscl)

d2A0dy02 = 4. * c1 * Lscl*Lscl * (8. * y0 * Lscl*Lscl + 8. * Lscl * d) / (4. * Lscl*Lscl &
      * x0*x0 - 4. * Lscl*Lscl * x0 + Lscl*Lscl + 4. * y0*y0 * Lscl*Lscl + &
       8. * y0 * Lscl * d + 4. * d*d) / (4. * Lscl*Lscl * x0*x0 + 4. * Lscl*Lscl &
       * x0 + Lscl*Lscl + 4. * y0*y0 * Lscl*Lscl + 8. * y0 * Lscl * d + 4. * d*d) & 
       - 4. * c1 * Lscl*Lscl * (-4. * Lscl*Lscl * x0*x0 + Lscl*Lscl + 4. * y0*y0  &
       * Lscl*Lscl + 8. * y0 * Lscl * d + 4. * d*d) / (4. * Lscl*Lscl * x0*x0 &
       - 4. * Lscl*Lscl * x0 + Lscl*Lscl + 4. * y0*y0 * Lscl*Lscl + 8. * y0 * Lscl *  &
      d + 4. * d*d) ** 2 / (4. * Lscl*Lscl * x0*x0 + 4. * Lscl*Lscl * x0 + Lscl*Lscl &
       + 4. * y0*y0 * Lscl*Lscl + 8. * y0 * Lscl * d + 4. * d*d) * (8. * &
       y0 * Lscl*Lscl + 8. * Lscl * d) - 4. * c1 * Lscl*Lscl * (-4. * Lscl*Lscl * x0*x0 &
       + Lscl*Lscl + 4. * y0*y0 * Lscl*Lscl + 8. * y0 * Lscl * d + 4. * d*d)  &
      / (4. * Lscl*Lscl * x0*x0 - 4. * Lscl*Lscl * x0 + Lscl*Lscl + 4. * y0*y0 * &
       Lscl*Lscl + 8. * y0 * Lscl * d + 4. * d*d) / (4. * Lscl*Lscl * x0*x0 + 4. &
       * Lscl*Lscl * x0 + Lscl*Lscl + 4. * y0*y0 * Lscl*Lscl + 8. * y0 * Lscl * d +  &
      4. * d*d) ** 2 * (8. * y0 * Lscl*Lscl + 8. * Lscl * d)

d2A0dx0dy0 = 32. * c1 * Lscl ** 4. * x0 / (4. * Lscl*Lscl * x0*x0 - 4. * Lscl*Lscl &
     * x0 + Lscl*Lscl + 4. * y0*y0 * Lscl*Lscl + 8. * y0 * Lscl * d + 4. * d*d) &
      / (4. * Lscl*Lscl * x0*x0 + 4. * Lscl*Lscl * x0 + Lscl*Lscl + 4. * y0*y0  &
      * Lscl*Lscl + 8. * y0 * Lscl * d + 4. * d*d) - 32. * c1 * (y0 * Lscl + d &
     ) * Lscl*Lscl*Lscl * x0 / (4. * Lscl*Lscl * x0*x0 - 4. * Lscl*Lscl * x0 + Lscl*Lscl &
      + 4. * y0*y0 * Lscl*Lscl + 8. * y0 * Lscl * d + 4. * d*d) ** 2 / (4. * &
      Lscl*Lscl * x0*x0 + 4. * Lscl*Lscl * x0 + Lscl*Lscl + 4. * y0*y0 * Lscl*Lscl  &
     + 8. * y0 * Lscl * d + 4. * d*d) * (8. * y0 * Lscl*Lscl + 8. * Lscl * d) - &
      32. * c1 * (y0 * Lscl + d) * Lscl*Lscl*Lscl * x0 / (4. * Lscl*Lscl * x0*x0 - 4.  &
     * Lscl*Lscl * x0 + Lscl*Lscl + 4. * y0*y0 * Lscl*Lscl + 8. * y0 * Lscl * d + 4. &
      * d*d) / (4. * Lscl*Lscl * x0*x0 + 4. * Lscl*Lscl * x0 + Lscl*Lscl + 4.  &
     * y0*y0 * Lscl*Lscl + 8. * y0 * Lscl * d + 4. * d*d) ** 2 * (8. * y0 * &
      Lscl*Lscl + 8. * Lscl * d)

d2A0dy0dx0 = -32. * c1 * Lscl ** 4. * x0 / (4. * Lscl*Lscl * x0*x0 - 4. * Lscl*Lscl &
      * x0 + Lscl*Lscl + 4. * y0*y0 * Lscl*Lscl + 8. * y0 * Lscl * d + 4. * d*d) &
       / (4. * Lscl*Lscl * x0*x0 + 4. * Lscl*Lscl * x0 + Lscl*Lscl + 4. * y0*y0 &
      * Lscl*Lscl + 8. * y0 * Lscl * d + 4. * d*d) - 4. * c1 * Lscl*Lscl * (- &
     4. * Lscl*Lscl * x0*x0 + Lscl*Lscl + 4. * y0*y0 * Lscl*Lscl + 8. * y0 * Lscl  &
     * d + 4. * d*d) / (4. * Lscl*Lscl * x0*x0 - 4. * Lscl*Lscl * x0 + Lscl*Lscl &
      + 4. * y0*y0 * Lscl*Lscl + 8. * y0 * Lscl * d + 4. * d*d) ** 2 / (4. &
      * Lscl*Lscl * x0*x0 + 4. * Lscl*Lscl * x0 + Lscl*Lscl + 4. * y0*y0 * Lscl*Lscl &
      + 8. * y0 * Lscl * d + 4. * d*d) * (8. * Lscl*Lscl * x0 - 4. * Lscl*Lscl &
     ) - 4. * c1 * Lscl*Lscl * (-4. * Lscl*Lscl * x0*x0 + Lscl*Lscl + 4. * y0*y0   &
     * Lscl*Lscl + 8. * y0 * Lscl * d + 4. * d*d) / (4. * Lscl*Lscl * x0*x0  &
     - 4. * Lscl*Lscl * x0 + Lscl*Lscl + 4. * y0*y0 * Lscl*Lscl + 8. * y0 * Lscl * d &
      + 4. * d*d) / (4. * Lscl*Lscl * x0*x0 + 4. * Lscl*Lscl * x0 + Lscl*Lscl  &
     + 4. * y0*y0 * Lscl*Lscl + 8. * y0 * Lscl * d + 4. * d*d) ** 2 * (8. *  &
     Lscl*Lscl * x0 + 4. * Lscl*Lscl)


END SUBROUTINE ddA0
!-------------------------------------------------------------------------------
SUBROUTINE B0setup(X0,Y0,Z0,B0)
!km: Subroutine for setting up the final field, B0
 REAL(num),INTENT(IN) :: X0,Y0,Z0
 REAL(num), DIMENSION(3), INTENT(OUT) :: B0
 REAL(num) :: a1,a2

 a1 = SQRT((X0+0.5)*(X0+0.5)+(Y0+1)*(Y0+1)+Z0*Z0)
 a2 = SQRT((X0-0.5)*(X0-0.5)+(Y0+1)*(Y0+1)+Z0*Z0)

 B0(1) = (X0+0.5)/a1/a1/a1 - (X0-0.5)/a2/a2/a2
 B0(2) = (Y0+1)/a1/a1/a1 - (Y0+1)/a2/a2/a2
 B0(3) = Z0/a1/a1/a1 - Z0/a2/a2/a2

 B0 = B0 * c1 !km: multiply by factor c1 in (36) of Grady (2009)

!km: Important to note that even though it looks like B0 ~ 1/L^2 from its 
!    definition, the definition of c1 means that c1 ~ L so our B0 ~ 1/L

END SUBROUTINE B0setup
!-------------------------------------------------------------------------------
SUBROUTINE dB0(X0,Y0,Z0,dB0dX0,dB0dY0,dB0dZ0)
!km: Subroutine to calculate 1st derivatives of B0 wrt modified co ords X0,Y0,Z0
 REAL(num), INTENT(IN) :: X0,Y0,Z0
 REAL(num), DIMENSION(3), INTENT(OUT) :: dB0dX0,dB0dY0,dB0dZ0
 REAL(num) :: a1,a2

 a1 = SQRT((X0+0.5)*(X0+0.5)+(Y0+1)*(Y0+1)+Z0*Z0)
 a2 = SQRT((X0-0.5)*(X0-0.5)+(Y0+1)*(Y0+1)+Z0*Z0)

 dB0dX0(1)=(a1*a1-3*(X0+0.5)*(X0+0.5))/a1/a1/a1/a1/a1 - (a2*a2-3*(X0-0.5)*(X0-0.5))/a2/a2/a2/a2/a2
 dB0dX0(2)=(-Y0-1)*3*(X0+0.5)/a1/a1/a1/a1/a1 - (-Y0-1)*3*(X0-0.5)/a2/a2/a2/a2/a2
 dB0dX0(3)=(-Z0)*3*(X0+0.5)/a1/a1/a1/a1/a1 - (-Z0)*3*(X0-0.5)/a2/a2/a2/a2/a2

 dB0dY0(1)=(-X0-0.5)*3*(Y0+1)/a1/a1/a1/a1/a1 - (-X0+0.5)*3*(Y0+1)/a2/a2/a2/a2/a2
 dB0dY0(2)=(a1*a1-3*(Y0+1)*(Y0+1))/a1/a1/a1/a1/a1 - (a2*a2-3*(Y0+1)*(Y0+1))/a2/a2/a2/a2/a2
 dB0dY0(3)=(-Z0)*3*(Y0+1)/a1/a1/a1/a1/a1 - (-Z0)*3*(Y0+1)/a2/a2/a2/a2/a2

 dB0dZ0(1)=(-X0-0.5)*3*Z0/a1/a1/a1/a1/a1 - (-X0+0.5)*3*Z0/a2/a2/a2/a2/a2
 dB0dZ0(2)=(-Y0-1)*3*Z0/a1/a1/a1/a1/a1 - (-Y0-1)*3*Z0/a2/a2/a2/a2/a2
 dB0dZ0(3)=(a1*a1-3*Z0*Z0)/a1/a1/a1/a1/a1 - (a2*a2-3*Z0*Z0)/a2/a2/a2/a2/a2
 
 !km: Again, multiply by c1 factor in Grady(2009) eq (36)
 dB0dX0 = dB0dX0 * c1
 dB0dY0 = dB0dY0 * c1
 dB0dZ0 = dB0dZ0 * c1

END SUBROUTINE dB0
!-------------------------------------------------------------------------------------------!
!SUBROUTINE ddA0 (X0,Y0,d2A0dx02, d2A0dy02, d2A0dx0dy0, d2A0dy0dx0 )
!
! REAL, INTENT(OUT) :: d2A0dx02, d2A0dy02, d2A0dx0dy0, d2A0dy0dx0
! REAL, INTENT(IN)  :: X0,Y0
!
! d2A0dx02 = 32 * c1 * (y0 * L + d) * L ** 3 / (4 * L ** 2 * x0 ** 2 - & 
!     4 * L ** 2 * x0 + L ** 2 + 4 * y0 ** 2 * L ** 2 + 8 * y0 * L * d + &
!      4 * d ** 2) / (4 * L ** 2 * x0 ** 2 + 4 * L ** 2 * x0 + L ** 2 +  &
!     4 * y0 ** 2 * L ** 2 + 8 * y0 * L * d + 4 * d ** 2) - 32 * c1 * (y&
!    &0 * L + d) * L ** 3 * x0 / (4 * L ** 2 * x0 ** 2 - 4 * L ** 2 * x0 &
!      + L ** 2 + 4 * y0 ** 2 * L ** 2 + 8 * y0 * L * d + 4 * d ** 2) ** &
!      2 / (4 * L ** 2 * x0 ** 2 + 4 * L ** 2 * x0 + L ** 2 + 4 * y0 **  & 
!     2 * L ** 2 + 8 * y0 * L * d + 4 * d ** 2) * (8 * L ** 2 * x0 - 4 * &
!      L ** 2) - 32 * c1 * (y0 * L + d) * L ** 3 * x0 / (4 * L ** 2 * x0 &
!      ** 2 - 4 * L ** 2 * x0 + L ** 2 + 4 * y0 ** 2 * L ** 2 + 8 * y0 * &
!      L * d + 4 * d ** 2) / (4 * L ** 2 * x0 ** 2 + 4 * L ** 2 * x0 + L &
!      ** 2 + 4 * y0 ** 2 * L ** 2 + 8 * y0 * L * d + 4 * d ** 2) ** 2 * &
!      (8 * L ** 2 * x0 + 4 * L ** 2)
!
!d2A0dy02 = 4 * c1 * L ** 2 * (8 * y0 * L ** 2 + 8 * L * d) / (4 * L &
!     ** 2 * x0 ** 2 - 4 * L ** 2 * x0 + L ** 2 + 4 * y0 ** 2 * L ** 2 + &
!      8 * y0 * L * d + 4 * d ** 2) / (4 * L ** 2 * x0 ** 2 + 4 * L ** 2 &
!      * x0 + L ** 2 + 4 * y0 ** 2 * L ** 2 + 8 * y0 * L * d + 4 * d ** & 
!     2) - 4 * c1 * L ** 2 * (-4 * L ** 2 * x0 ** 2 + L ** 2 + 4 * y0 ** &
!      2 * L ** 2 + 8 * y0 * L * d + 4 * d ** 2) / (4 * L ** 2 * x0 ** 2 &
!      - 4 * L ** 2 * x0 + L ** 2 + 4 * y0 ** 2 * L ** 2 + 8 * y0 * L *  &
!     d + 4 * d ** 2) ** 2 / (4 * L ** 2 * x0 ** 2 + 4 * L ** 2 * x0 + L &
!      ** 2 + 4 * y0 ** 2 * L ** 2 + 8 * y0 * L * d + 4 * d ** 2) * (8 * &
!      y0 * L ** 2 + 8 * L * d) - 4 * c1 * L ** 2 * (-4 * L ** 2 * x0 ** &
!      2 + L ** 2 + 4 * y0 ** 2 * L ** 2 + 8 * y0 * L * d + 4 * d ** 2)  &
!     / (4 * L ** 2 * x0 ** 2 - 4 * L ** 2 * x0 + L ** 2 + 4 * y0 ** 2 * &
!      L ** 2 + 8 * y0 * L * d + 4 * d ** 2) / (4 * L ** 2 * x0 ** 2 + 4 &
!      * L ** 2 * x0 + L ** 2 + 4 * y0 ** 2 * L ** 2 + 8 * y0 * L * d +  &
!     4 * d ** 2) ** 2 * (8 * y0 * L ** 2 + 8 * L * d)
!
!d2A0dx0dy0 = 32 * c1 * L ** 4 * x0 / (4 * L ** 2 * x0 ** 2 - 4 * L ** &
!     2 * x0 + L ** 2 + 4 * y0 ** 2 * L ** 2 + 8 * y0 * L * d + 4 * d ** &
!      2) / (4 * L ** 2 * x0 ** 2 + 4 * L ** 2 * x0 + L ** 2 + 4 * y0 ** &
!      2 * L ** 2 + 8 * y0 * L * d + 4 * d ** 2) - 32 * c1 * (y0 * L + d &
!     ) * L ** 3 * x0 / (4 * L ** 2 * x0 ** 2 - 4 * L ** 2 * x0 + L ** 2 &
!      + 4 * y0 ** 2 * L ** 2 + 8 * y0 * L * d + 4 * d ** 2) ** 2 / (4 * &
!      L ** 2 * x0 ** 2 + 4 * L ** 2 * x0 + L ** 2 + 4 * y0 ** 2 * L **  &
!     2 + 8 * y0 * L * d + 4 * d ** 2) * (8 * y0 * L ** 2 + 8 * L * d) - &
!      32 * c1 * (y0 * L + d) * L ** 3 * x0 / (4 * L ** 2 * x0 ** 2 - 4  &
!     * L ** 2 * x0 + L ** 2 + 4 * y0 ** 2 * L ** 2 + 8 * y0 * L * d + 4 &
!      * d ** 2) / (4 * L ** 2 * x0 ** 2 + 4 * L ** 2 * x0 + L ** 2 + 4  &
!     * y0 ** 2 * L ** 2 + 8 * y0 * L * d + 4 * d ** 2) ** 2 * (8 * y0 * &
!      L ** 2 + 8 * L * d)
!
!d2A0dy0dx0 = -32 * c1 * L ** 4 * x0 / (4 * L ** 2 * x0 ** 2 - 4 * L ** &
!      2 * x0 + L ** 2 + 4 * y0 ** 2 * L ** 2 + 8 * y0 * L * d + 4 * d ** &
!      2) / (4 * L ** 2 * x0 ** 2 + 4 * L ** 2 * x0 + L ** 2 + 4 * y0 **&
!      2 * L ** 2 + 8 * y0 * L * d + 4 * d ** 2) - 4 * c1 * L ** 2 * (- &
!     4 * L ** 2 * x0 ** 2 + L ** 2 + 4 * y0 ** 2 * L ** 2 + 8 * y0 * L  &
!     * d + 4 * d ** 2) / (4 * L ** 2 * x0 ** 2 - 4 * L ** 2 * x0 + L ** &
!      2 + 4 * y0 ** 2 * L ** 2 + 8 * y0 * L * d + 4 * d ** 2) ** 2 / (4 &
!      * L ** 2 * x0 ** 2 + 4 * L ** 2 * x0 + L ** 2 + 4 * y0 ** 2 * L **&
!      2 + 8 * y0 * L * d + 4 * d ** 2) * (8 * L ** 2 * x0 - 4 * L ** 2 &
!     ) - 4 * c1 * L ** 2 * (-4 * L ** 2 * x0 ** 2 + L ** 2 + 4 * y0 **  &
!     2 * L ** 2 + 8 * y0 * L * d + 4 * d ** 2) / (4 * L ** 2 * x0 ** 2  &
!     - 4 * L ** 2 * x0 + L ** 2 + 4 * y0 ** 2 * L ** 2 + 8 * y0 * L * d &
!      + 4 * d ** 2) / (4 * L ** 2 * x0 ** 2 + 4 * L ** 2 * x0 + L ** 2  &
!     + 4 * y0 ** 2 * L ** 2 + 8 * y0 * L * d + 4 * d ** 2) ** 2 * (8 *  &
!     L ** 2 * x0 + 4 * L ** 2)
!
!
!END SUBROUTINE ddA0

FUNCTION sech(x)

REAL(num) :: sech, x

sech = 1.0/cosh(x)

return

END FUNCTION sech

FUNCTION sec(x)

REAL(num) :: sec, x

sec = 1.0/cos(x)

return

END FUNCTION sec

SUBROUTINE ABSETUP(R,T,X0,Y0,Z0,DX0DX,DX0DY,DX0DZ,DX0DT,DY0DX,DY0DY,DY0DZ,DY0DT,DZ0DX,DZ0DY,DZ0DZ,DZ0DT,&
           D2X0DX2,D2X0DXDY,D2X0DXDZ,D2X0DXDT,D2X0DY2,D2X0DYDZ,D2X0DYDT,D2X0DZ2,D2X0DZDT,D2X0DT2,&
           D2Y0DX2,D2Y0DXDY,D2Y0DXDZ,D2Y0DXDT,D2Y0DY2,D2Y0DYDZ,D2Y0DYDT,D2Y0DZ2,D2Y0DZDT,D2Y0DT2,&
           D2Z0DX2,D2Z0DXDY,D2Z0DXDZ,D2Z0DXDT,D2Z0DY2,D2Z0DYDZ,D2Z0DYDT,D2Z0DZ2,D2Z0DZDT,D2Z0DT2,&
           D2X0DYDX,D2X0DZDX,D2X0DTDX,D2X0DZDY,D2X0DTDY,D2X0DTDZ,&
           D2Y0DYDX,D2Y0DZDX,D2Y0DTDX,D2Y0DZDY,D2Y0DTDY,D2Y0DTDZ,&
           D2Z0DYDX,D2Z0DZDX,D2Z0DTDX,D2Z0DZDY,D2Z0DTDY,D2Z0DTDZ)

REAL(num), INTENT(IN) :: t
REAL(num), DIMENSION(3), INTENT(IN) :: R
REAL(num), INTENT(OUT) :: X0,Y0,Z0
REAL(num), INTENT(OUT) :: DX0DX,DX0DY,DX0DZ,DX0DT,DY0DX,DY0DY,DY0DZ,DY0DT,DZ0DX,DZ0DY,DZ0DZ,DZ0DT
REAL(num), INTENT(OUT) :: D2X0DX2,D2X0DXDY,D2X0DXDZ,D2X0DXDT,D2X0DY2,D2X0DYDZ,D2X0DYDT,D2X0DZ2,D2X0DZDT,D2X0DT2
REAL(num), INTENT(OUT) :: D2Y0DX2,D2Y0DXDY,D2Y0DXDZ,D2Y0DXDT,D2Y0DY2,D2Y0DYDZ,D2Y0DYDT,D2Y0DZ2,D2Y0DZDT,D2Y0DT2
REAL(num), INTENT(OUT) :: D2Z0DX2,D2Z0DXDY,D2Z0DXDZ,D2Z0DXDT,D2Z0DY2,D2Z0DYDZ,D2Z0DYDT,D2Z0DZ2,D2Z0DZDT,D2Z0DT2
REAL(num), INTENT(OUT) :: D2X0DYDX,D2X0DZDX,D2X0DTDX,D2X0DZDY,D2X0DTDY,D2X0DTDZ
REAL(num), INTENT(OUT) :: D2Y0DYDX,D2Y0DZDX,D2Y0DTDX,D2Y0DZDY,D2Y0DTDY,D2Y0DTDZ
REAL(num), INTENT(OUT) :: D2Z0DYDX,D2Z0DZDX,D2Z0DTDX,D2Z0DZDY,D2Z0DTDY,D2Z0DTDZ

REAL(num) :: lilf, dlilfdt, d2lilfdt2
REAL(num) :: lilg, dlilgdt, d2lilgdt2
REAL(num) :: F, dFdy, d2Fdy2
REAL(num) :: G, dGdr, d2Gdr2
REAL(num) :: J,dJdr,dJdy,d2Jdr2,d2Jdrdy,d2Jdy2
REAL(num) :: bigT,dbigTdr,dbigTdy,dbigTdt,d2bigTdr2,d2bigTdrdy,d2bigTdrdt,d2bigTdy2,d2bigTdydt,d2bigTdt2
REAL(num) :: phi,dphidr,dphidy,dphidt,d2phidr2,d2phidrdy,d2phidrdt,d2phidy2,d2phidydt,d2phidt2
REAL(num) :: ra,drdx, d2rdx2, drdz, d2rdz2, d2rdxdz
REAL(num) :: dY0dr, d2Y0dr2, d2Y0drdy, d2Y0drdt
REAL(num) :: dbigTdx, d2bigTdx2, d2bigTdxdr, d2bigTdxdy, d2bigTdxdt
REAL(num) :: dphidx, d2phidx2, d2phidxdr, d2phidxdy, d2phidxdt
REAL(num) :: partdY0dx, partd2Y0dx2, partd2Y0dxdr, partd2Y0dxdy, partd2Y0dxdt

REAL(num), PARAMETER :: lild=0.3,w=0.1,w2=2.3,beta=0.5,k=7.0,lils=0.7,szero=0.5
REAL(num), PARAMETER :: alpha=1.0,chi=1.0,zeta=0.3
REAL(num), PARAMETER :: sigma=4.0,vphi=(-2.0),yzero=7.0,xprm=100.0,yprm=1.0,tprm=5.0,tzero=100.0

REAL(num) :: den, oneuponden
REAL(num), PARAMETER :: ay=1.0/2.0, yh=0.0, delta=0.0, a3d=1.0

ra=SQRT(R(1)*R(1)+R(3)*R(3))
oneuponden=1.0/(a3d*a3d+R(1)*R(1)+ay*(R(2)-yh)*ay*(R(2)-yh)+R(3)*R(3))
den=a3d*a3d+R(1)*R(1)+ay*(R(2)-yh)*ay*(R(2)-yh)+R(3)*R(3)

drdx = R(1)/SQRT(R(1)*R(1)+R(3)*R(3))
drdz = R(3)/SQRT(R(1)*R(1)+R(3)*R(3))
d2rdx2 = R(3)*R(3)/SQRT(R(1)*R(1)+R(3)*R(3))/(R(1)*R(1)+R(3)*R(3))
d2rdxdz = (-R(1)*R(3))/SQRT(R(1)*R(1)+R(3)*R(3))/(R(1)*R(1)+R(3)*R(3))
d2rdx2 = R(1)*R(1)/SQRT(R(1)*R(1)+R(3)*R(3))/(R(1)*R(1)+R(3)*R(3))

lilf = (1.0-tanh(t-tzero))/2.0
dlilfdt = (-sech(t-tzero)*sech(t-tzero))/2.0
d2lilfdt2 = tanh(t-tzero)*sech(t-tzero)*sech(t-tzero)

lilg = 1.0-lilf
dlilgdt = (-dlilfdt)
d2lilgdt2 = (-d2lilfdt2)

F = (1.0-tanh(R(2)-yprm))/2.0
dFdy = (-sech(R(2)-yprm)*sech(R(2)-yprm))/2.0
d2Fdy2 = tanh(R(2)-yprm)*sech(R(2)-yprm)*sech(R(2)-yprm)

G =  (tanh(ra+xprm)-tanh(ra-xprm))/2.0
dGdr =  (sech(ra+xprm)*sech(ra+xprm)-sech(ra-xprm)*sech(ra-xprm))/2.0
d2Gdr2 =  tanh(ra-xprm)*sech(ra-xprm)*sech(ra-xprm)-tanh(ra+xprm)*sech(ra+xprm)*sech(ra+xprm)

J = lild*EXP(-ra*ra*R(2)/w)
dJdr = (-2.0*ra*R(2)/w)*lild*EXP(-ra*ra*R(2)/w)
dJdy = (-ra*ra/w)*lild*EXP(-ra*ra*R(2)/w)
d2Jdr2 = 2*R(2)/w*(2*ra*ra*R(2)/w-1.0)*lild*EXP(-ra*ra*R(2)/w)
d2Jdrdy = 2*ra/w*(ra*ra*R(2)/w-1.0)*lild*EXP(-ra*ra*R(2)/w)
d2Jdy2 =ra*ra*ra*ra/w/w*lild*EXP(-ra*ra*R(2)/w)

bigT = k*(1.0+(chi*sin(pi*R(2)/yzero)-1.0)*(1.0-tanh(zeta*(t-tprm)))/2.0)* &
       tan(pi/2.0*R(1)*R(1)/w2)*tanh(R(2))
dbigTdx = k*(1.0+(chi*sin(pi*R(2)/yzero)-1.0)*(1.0-tanh(zeta*(t-tprm)))/2.0)* &
          R(1)*pi/w2*sec(pi/2.0*R(1)*R(1)/w2)*sec(pi/2.0*R(1)*R(1)/w2)*tanh(R(2))
dbigTdr = 0.0
dbigTdy = k*chi*pi/yzero*cos(pi*R(2)/yzero)*(1.0-tanh(zeta*(t-tprm)))/2.0 * &
          tan(pi/2.0*R(1)*R(1)/w2)*tanh(R(2)) + &
          k*(1.0+(chi*sin(pi*R(2)/yzero)-1.0)*(1.0-tanh(zeta*(t-tprm)))/2.0)* &
          tan(pi/2.0*R(1)*R(1)/w2)*sech(R(2))*sech(R(2))
dbigTdt = k*(chi*sin(pi*R(2)/yzero)-1.0)* &
          (-zeta*sech(zeta*(t-tprm))*sech(zeta*(t-tprm)))/2.0* &
          tan(pi/2.0*R(1)*R(1)/w2)*tanh(R(2))
d2bigTdx2 = k*(1.0+(chi*sin(pi*R(2)/yzero)-1.0)*(1.0-tanh(zeta*(t-tprm)))/2.0)* &
            pi/w2*tanh(R(2))*(sec(pi/2.0*R(1)*R(1)/w2)*sec(pi/2.0*R(1)*R(1)/w2) + &
            R(1)*2*pi*R(1)/w2*sec(pi/2.0*R(1)*R(1)/w2)*sec(pi/2.0*R(1)*R(1)/w2)* &
            tan(pi/2.0*R(1)*R(1)/w2))
d2bigTdxdr = 0.0
d2bigTdxdy = k*chi*pi/yzero*cos(pi*R(2)/yzero)*(1.0-tanh(zeta*(t-tprm)))/2.0 * &
             R(1)*pi/w2*sec(pi/2.0*R(1)*R(1)/w2)*sec(pi/2.0*R(1)*R(1)/w2)*tanh(R(2)) + &
             k*(1.0+(chi*sin(pi*R(2)/yzero)-1.0)*(1.0-tanh(zeta*(t-tprm)))/2.0)* &
             R(1)*pi/w2*sec(pi/2.0*R(1)*R(1)/w2)*sec(pi/2.0*R(1)*R(1)/w2)* &
             sech(R(2))*sech(R(2))
d2bigTdxdt = k*(chi*sin(pi*R(2)/yzero)-1.0)* &
             (-zeta*sech(zeta*(t-tprm))*sech(zeta*(t-tprm)))/2.0* &
             R(1)*pi/w2*sec(pi/2.0*R(1)*R(1)/w2)*sec(pi/2.0*R(1)*R(1)/w2)*tanh(R(2))
d2bigTdr2 = 0.0
d2bigTdrdy = 0.0
d2bigTdrdt = 0.0
d2bigTdy2 = (-k)*chi*pi*pi/yzero/yzero*sin(pi*R(2)/yzero)* &
            (1.0-tanh(zeta*(t-tprm)))/2.0*tan(pi/2.0*R(1)*R(1)/w2)*tanh(R(2)) + &
            2*k*chi*pi/yzero*cos(pi*R(2)/yzero)*(1.0-tanh(zeta*(t-tprm)))/2.0 *&
            tan(pi/2.0*R(1)*R(1)/w2)*sech(R(2))*sech(R(2)) + &
            k*(1.0+(chi*sin(pi*R(2)/yzero)-1.0)*(1.0-tanh(zeta*(t-tprm)))/2.0)* &
            tan(pi/2.0*R(1)*R(1)/w2)*(-2*tanh(R(2)))*sech(R(2))*sech(R(2))
d2bigTdydt = k*chi*pi/yzero*cos(pi*R(2)/yzero)* &
             (-zeta*sech(zeta*(t-tprm))*sech(zeta*(t-tprm)))/2.0 * &
             tan(pi/2.0*R(1)*R(1)/w2)*tanh(R(2)) + &
             k*(chi*sin(pi*R(2)/yzero)-1.0)* &
             (-zeta*sech(zeta*(t-tprm))*sech(zeta*(t-tprm)))/2.0* &
             tan(pi/2.0*R(1)*R(1)/w2)*sech(R(2))*sech(R(2))
d2bigTdt2 = k*(chi*sin(pi*R(2)/yzero)-1.0)* &
            zeta*zeta*tanh(zeta*(t-tprm))*sech(zeta*(t-tprm))*sech(zeta*(t-tprm))* &
            tan(pi/2.0*R(1)*R(1)/w2)*tanh(R(2))

phi = alpha*(1.0+(chi*sin(pi*R(2)/yzero)-1.0)*(1.0-tanh(zeta*(t-tprm)))/2.0)* &
      (R(2)+1.0)**beta*(szero*ra*ra+1.0)**beta* &
      (R(2)-vphi*sigma*tanh(t/sigma)-yzero+J)+bigT
dphidx = dbigTdx
dphidr = alpha*(1.0+(chi*sin(pi*R(2)/yzero)-1.0)*(1.0-tanh(zeta*(t-tprm)))/2.0)* &
         (R(2)+1.0)**beta*beta*2*ra*szero*(szero*ra*ra+1.0)**(beta-1.0)* &
         (R(2)-vphi*sigma*tanh(t/sigma)-yzero+J) + &
         alpha*(1.0+(chi*sin(pi*R(2)/yzero)-1.0)*(1.0-tanh(zeta*(t-tprm)))/2.0)* &
         (R(2)+1.0)**beta*(szero*ra*ra+1.0)**beta*DJDr + dbigTdr
dphidy = alpha*(pi/yzero*chi*cos(pi*R(2)/yzero)*(1.0-tanh(zeta*(t-tprm)))/2.0)* &
         (R(2)+1.0)**beta*(szero*ra*ra+1.0)**beta* &
         (R(2)-vphi*sigma*tanh(t/sigma)-yzero+J) + &
         alpha*(1.0+(chi*sin(pi*R(2)/yzero)-1.0)*(1.0-tanh(zeta*(t-tprm)))/2.0)* &
         beta*(R(2)+1.0)**(beta-1.0)*(szero*ra*ra+1.0)**beta* &
         (R(2)-vphi*sigma*tanh(t/sigma)-yzero+J) + &
         alpha*(1.0+(chi*sin(pi*R(2)/yzero)-1.0)*(1.0-tanh(zeta*(t-tprm)))/2.0)* &
         (R(2)+1.0)**beta*(szero*ra*ra+1.0)**beta*(1.0+dJdy) + dbigTdy
dphidt = alpha*(chi*sin(pi*R(2)/yzero)-1.0)* &
         (-zeta*sech(zeta*(t-tprm))*sech(zeta*(t-tprm)))/2.0* &
         (R(2)+1.0)**beta*(szero*ra*ra+1.0)**beta* &
         (R(2)-vphi*sigma*tanh(t/sigma)-yzero+J) + &
         alpha*(1.0+(chi*sin(pi*R(2)/yzero)-1.0)*(1.0-tanh(zeta*(t-tprm)))/2.0)* &
         (R(2)+1.0)**beta*(szero*ra*ra+1.0)**beta* &
         (-vphi*sech(t/sigma)*sech(t/sigma))+dbigTdt
d2phidx2 = d2bigTdx2
d2phidxdr = d2bigTdxdr
d2phidxdy = d2bigTdxdy
d2phidxdt = d2bigTdxdt
d2phidr2 = (alpha*(1.0+(chi*sin(pi*R(2)/yzero)-1.0)*(1.0-tanh(zeta*(t-tprm)))/2.0)* &
           (R(2)+1.0)**beta) * &
           (2*beta*szero*(szero*ra*ra+1.0)**(beta-1.0)* &
           (R(2)-vphi*sigma*tanh(t/sigma)-yzero+J) + &
           4*beta*(beta-1)*szero*szero*ra*ra*(szero*ra*ra+1.0)**(beta-2.0)* &
           (R(2)-vphi*sigma*tanh(t/sigma)-yzero+J) + &
           4*beta*ra*szero*(szero*ra*ra+1.0)**(beta-1.0)*dJdr + &
           (szero*ra*ra+1.0)**beta*d2Jdr2) + d2bigTdr2
d2phidrdy = ((alpha*(pi/yzero*chi*cos(pi*R(2)/yzero)*(1.0-tanh(zeta*(t-tprm)))/2.0)* &
            (R(2)+1.0)**beta) + &
            (alpha*(1.0+(chi*sin(pi*R(2)/yzero)-1.0)*(1.0-tanh(zeta*(t-tprm)))/2.0)* &
            beta*(R(2)+1.0)**(beta-1.0)))* &
            (2*beta*ra*szero*(szero*ra*ra+1.0)**(beta-1.0)* &
            (R(2)-vphi*sigma*tanh(t/sigma)-yzero+J) + & 
            (szero*ra*ra+1.0)**beta*DJDr) + &
            (alpha*(1.0+(chi*sin(pi*R(2)/yzero)-1.0)*(1.0-tanh(zeta*(t-tprm)))/2.0)* &
            (R(2)+1.0)**beta)* &
            (2*beta*ra*szero*(szero*ra*ra+1.0)**(beta-1.0)*(1.0+dJdy) + &
            (szero*ra*ra+1.0)**beta*d2Jdrdy) + d2bigTdrdy
d2phidrdt = (alpha*(chi*sin(pi*R(2)/yzero)-1.0)* &
            (-zeta*sech(zeta*(t-tprm))*sech(zeta*(t-tprm)))/2.0* &
            (R(2)+1.0)**beta)* &
            (2*beta*ra*szero*(szero*ra*ra+1.0)**(beta-1.0)* &
            (R(2)-vphi*sigma*tanh(t/sigma)-yzero+J) + & 
            (szero*ra*ra+1.0)**beta*DJDr) + &
            (alpha*(1.0+(chi*sin(pi*R(2)/yzero)-1.0)*(1.0-tanh(zeta*(t-tprm)))/2.0)* &
            (R(2)+1.0)**beta)* &
            (2*beta*ra*szero*(szero*ra*ra+1.0)**(beta-1.0)* &
            (-vphi*sech(t/sigma)*sech(t/sigma))) + d2bigTdrdt
d2phidy2 = (alpha*(pi*pi/yzero/yzero*chi*(-sin(pi*R(2)/yzero))*(1.0-tanh(zeta*(t-tprm)))/2.0)* &
           (R(2)+1.0)**beta*(szero*ra*ra+1.0)**beta + &
           alpha*pi/yzero*chi*cos(pi*R(2)/yzero)*(1.0-tanh(zeta*(t-tprm)))/2.0* &
           beta*(R(2)+1.0)**(beta-1.0)*(szero*ra*ra+1.0)**beta)*&
           (R(2)-vphi*sigma*tanh(t/sigma)-yzero+J) + &
           alpha*pi/yzero*chi*cos(pi*R(2)/yzero)*(1.0-tanh(zeta*(t-tprm)))/2.0* &
           (R(2)+1.0)**beta*(szero*ra*ra+1.0)**beta*&
           (1.0+dJdy) + &
           (alpha*pi/yzero*chi*cos(pi*R(2)/yzero)*(1.0-tanh(zeta*(t-tprm)))/2.0* &
           beta*(R(2)+1.0)**(beta-1.0)*(szero*ra*ra+1.0)**beta + &
           alpha*(1.0+(chi*sin(pi*R(2)/yzero)-1.0)*(1.0-tanh(zeta*(t-tprm)))/2.0)* &
           beta*(beta-1.0)*(R(2)+1.0)**(beta-2.0)*(szero*ra*ra+1.0)**beta) * &
           (R(2)-vphi*sigma*tanh(t/sigma)-yzero+J) + &
           alpha*(1.0+(chi*sin(pi*R(2)/yzero)-1.0)*(1.0-tanh(zeta*(t-tprm)))/2.0)* &
           beta*(R(2)+1.0)**(beta-1.0)*(szero*ra*ra+1.0)**beta*(1.0+dJdy) + &
           (alpha*(pi/yzero*chi*cos(pi*R(2)/yzero)*(1.0-tanh(zeta*(t-tprm)))/2.0)* &
           (R(2)+1.0)**beta*(szero*ra*ra+1.0)**beta+ &
           alpha*(1.0+(chi*sin(pi*R(2)/yzero)-1.0)*(1.0-tanh(zeta*(t-tprm)))/2.0)* &
           beta*(R(2)+1.0)**(beta-1.0)*(szero*ra*ra+1.0)**beta)*(1.0+dJdy) + &
           alpha*(1.0+(chi*sin(pi*R(2)/yzero)-1.0)*(1.0-tanh(zeta*(t-tprm)))/2.0)* &
           (R(2)+1.0)**beta*(szero*ra*ra+1.0)**beta*d2Jdy2 + d2bigTdy2
d2phidydt = alpha*(pi/yzero*chi*cos(pi*R(2)/yzero)* &
            (-zeta*sech(zeta*(t-tprm))*sech(zeta*(t-tprm)))/2.0)* &
            (R(2)+1.0)**beta*(szero*ra*ra+1.0)**beta* &
            (R(2)-vphi*sigma*tanh(t/sigma)-yzero+J) + &
            alpha*(pi/yzero*chi*cos(pi*R(2)/yzero)*(1.0-tanh(zeta*(t-tprm)))/2.0)* &
            (R(2)+1.0)**beta*(szero*ra*ra+1.0)**beta* &
            (-vphi*sech(t/sigma)*sech(t/sigma)) + &
            alpha*(chi*sin(pi*R(2)/yzero)-1.0)* &
            (-zeta*sech(zeta*(t-tprm))*sech(zeta*(t-tprm)))/2.0* &
            beta*(R(2)+1.0)**(beta-1.0)*(szero*ra*ra+1.0)**beta* &
            (R(2)-vphi*sigma*tanh(t/sigma)-yzero+J) + &
            alpha*(1.0+(chi*sin(pi*R(2)/yzero)-1.0)*(1.0-tanh(zeta*(t-tprm)))/2.0)* &
            beta*(R(2)+1.0)**(beta-1.0)*(szero*ra*ra+1.0)**beta* &
            (-vphi*sech(t/sigma)*sech(t/sigma)) + &
            alpha*(chi*sin(pi*R(2)/yzero)-1.0)* &
            (-zeta*sech(zeta*(t-tprm))*sech(zeta*(t-tprm)))/2.0* &
            (R(2)+1.0)**beta*(szero*ra*ra+1.0)**beta* &
            (1.0+dJdy) + d2bigTdydt
d2phidt2 = alpha*(chi*sin(pi*R(2)/yzero)-1.0)* &
           zeta*zeta*tanh(zeta*(t-tprm))*sech(zeta*(t-tprm))*sech(zeta*(t-tprm))* &
           (R(2)+1.0)**beta*(szero*ra*ra+1.0)**beta* &
           (R(2)-vphi*sigma*tanh(t/sigma)-yzero+J) + &
           2*alpha*(chi*sin(pi*R(2)/yzero)-1.0)* &
           (-zeta*sech(zeta*(t-tprm))*sech(zeta*(t-tprm)))/2.0* &
           (R(2)+1.0)**beta*(szero*ra*ra+1.0)**beta* &
           (-vphi*sech(t/sigma)*sech(t/sigma)) + &
           alpha*(1.0+(chi*sin(pi*R(2)/yzero)-1.0)* &
           (1.0-tanh(zeta*(t-tprm)))/2.0)*&
           (R(2)+1.0)**beta*(szero*ra*ra+1.0)**beta* &
           (2*vphi/sigma*tanh(t/sigma)*sech(t/sigma)*sech(t/sigma))+d2bigTdt2

Y0 = (lils*LOG(1.0+R(2)/lils)*(1.0-F*G)+R(2)*F*G+R(2)*(1.0+tanh(phi))/2.0)*lilf+&
     R(2)*lilg
X0 = R(1)-delta*(Y0-R(2))*R(3)*oneuponden 
Z0 = R(3)+delta*(Y0-R(2))*R(1)*oneuponden

dY0dr = (lils*LOG(1.0+R(2)/lils)*(-F*dGdr)+R(2)*F*dGdr+ &
        R(2)*dphidr*sech(phi)*sech(phi)/2.0)*lilf
partdY0dx = R(2)*dphidx*sech(phi)*sech(phi)/2.0*lilf
partd2Y0dx2 = R(2)*(d2phidx2*sech(phi)*sech(phi)/2.0- &
                    dphidx*dphidx*tanh(phi)*sech(phi)*sech(phi))*lilf
partd2Y0dxdr = R(2)*(d2phidxdr*sech(phi)*sech(phi)/2.0- &
                    dphidx*dphidr*tanh(phi)*sech(phi)*sech(phi))*lilf
partd2Y0dxdy = (dphidx*sech(phi)*sech(phi)/2.0+ &
                R(2)*d2phidxdy*sech(phi)*sech(phi)/2.0 - &
                R(2)*dphidx*dphidy*tanh(phi)*sech(phi)*sech(phi))*lilf
partd2Y0dxdt = R(2)*(d2phidxdt*sech(phi)*sech(phi)/2.0*lilf - &
                     dphidx*dphidt*tanh(phi)*sech(phi)*sech(phi)*lilf + &
                     dphidx*sech(phi)*sech(phi)/2.0*dlilfdt)

dY0dx = drdx*dY0dr + partdY0dx
dY0dy = (1.0/(1.0+(R(2)/lils))*(1.0-F*G)+lils*LOG(1.0+R(2)/lils)*(-dFdy*G)+&
        F*G+R(2)*dFdy*G+(1+tanh(phi))/2.0+ &
        R(2)*dphidy*sech(phi)*sech(phi)/2.0)*lilf+lilg
dY0dz = drdz*dY0dr
dY0dt = (lils*LOG(1.0+R(2)/lils)*(1.0-F*G)+R(2)*F*G+ &
        R(2)*(1.0+tanh(phi))/2.0)*dlilfdt+R(2)*dphidt*sech(phi)*sech(phi)/2.0*lilf+ &
        R(2)*dlilgdt

dX0dX=1.0+2.0*R(1)*delta*(Y0-R(2))*R(3)*oneuponden*oneuponden - &
      delta*dY0dX*R(3)*oneuponden
dX0dY=delta*(1.0-dY0dY)*R(3)*oneuponden+ &
      2.0*delta*(Y0-R(2))*ay*ay*(R(2)-yh)*R(3)*oneuponden*oneuponden
dX0dZ=(-1.0)*delta*(Y0-R(2))*(a3d*a3d+R(1)*R(1)+ay*(R(2)-yh)*ay*(R(2)-yh)-R(3)*R(3))* &
                            oneuponden*oneuponden - &
      delta*dY0dZ*R(3)*oneuponden 
dX0dt = (-delta)*dY0dt*R(3)*oneuponden

dZ0dX=delta*(Y0-R(2))*(a3d*a3d+ay*(R(2)-yh)*ay*(R(2)-yh)+R(3)*R(3)-R(1)*R(1))* &
                             oneuponden*oneuponden + &
      delta*dY0dX*R(1)*oneuponden
dZ0dY=delta*(dY0dY-1.0)*R(1)*oneuponden-2.0*delta*(Y0-R(2))*ay*ay*(R(2)-yh)*R(1)*oneuponden*oneuponden
dZ0dZ=1.0-2*delta*(Y0-R(2))*R(1)*R(3)*oneuponden*oneuponden + &
       delta*dY0dZ*R(1)*oneuponden 
dZ0dt = delta*dY0dt*R(1)*oneuponden

d2Y0dr2 = (lils*LOG(1.0+R(2)/lils)*(-F*d2Gdr2)+R(2)*F*d2Gdr2+ &
        R(2)*d2phidr2*sech(phi)*sech(phi)/2.0 + &
        R(2)*dphidr*dphidr*(-tanh(phi))*sech(phi)*sech(phi))*lilf
d2Y0drdy = (1.0/(1.0+(R(2)/lils))*(-F*dGdr)+lils*LOG(1.0+R(2)/lils)*&
           (-dFdy)*dGdr+F*dGdr+R(2)*dFdy*dGdr+dphidr*sech(phi)*sech(phi)/2.0+&
           R(2)*(d2phidrdy*sech(phi)*sech(phi)/2.0+ &
           dphidr*dphidy*(-tanh(phi))*sech(phi)*sech(phi)))*lilf
d2Y0drdt = (R(2)*d2phidrdt*sech(phi)*sech(phi)/2.0+R(2)*dphidr*dphidt* &
           (-tanh(phi))*sech(phi)*sech(phi))*lilf + &
           (lils*LOG(1.0+R(2)/lils)*(-F*dGdr)+R(2)*F*dGdr+ &
           R(2)*dphidr*sech(phi)*sech(phi)/2.0)*dlilfdt

d2Y0dx2 = d2rdx2*dY0dr + drdx*drdx*d2Y0dr2 + partd2Y0dx2 + 2*drdx*partd2Y0dxdr
d2Y0dxdy = drdx*d2Y0drdy + partd2Y0dxdy
d2Y0dydx = d2Y0dxdy
d2Y0dxdz = d2rdxdz*dY0dr + drdx*drdz*d2Y0dr2 + partd2Y0dxdr*drdz
d2Y0dzdx = d2Y0dxdz
d2Y0dy2 = ((-1.0)/lils/(1+(R(2)/lils))/(1+(R(2)/lils))*(1.0-F*G)+&
          2.0/(1.0+(R(2)/lils))*(-dFdy*G)+lils*LOG(1.0+R(2)/lils)*(-d2Fdy2*G)+&
          2.0*dFdy*G+R(2)*d2Fdy2*G+dphidy*sech(phi)*sech(phi)+R(2)*&
          (d2phidy2*sech(phi)*sech(phi)/2.0-&
          dphidy*dphidy*tanh(phi)*sech(phi)*sech(phi)))*lilf
d2Y0dydz = drdz*d2Y0drdy
d2Y0dzdy = d2Y0dydz
d2Y0dz2 = d2rdz2*dY0dr + drdz*drdz*d2Y0dr2

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

d2Y0dxdt = drdx*d2Y0drdt + partd2Y0dxdt
d2Y0dtdx = d2Y0dxdt
d2Y0dydt = (dphidt*sech(phi)*sech(phi)/2.0+R(2)*d2phidydt*sech(phi)*sech(phi)/2.0-&
           R(2)*dphidy*dphidt*tanh(phi)*sech(phi)*sech(phi))*lilf + &
           (1.0/(1.0+(R(2)/lils))*(1.0-F*G)+lils*LOG(1.0+R(2)/lils)*(-dFdy*G)+&
           F*G+R(2)*dFdy*G+(1+tanh(phi))/2.0+ &
           R(2)*dphidy*sech(phi)*sech(phi)/2.0)*dlilfdt+dlilgdt
d2Y0dtdy = d2Y0dydt
d2Y0dzdt = drdz*d2Y0drdt
d2Y0dtdz = d2Y0dzdt
d2Y0dt2 = R(2)*dphidt*sech(phi)*sech(phi)*dlilfdt + &
          R(2)*sech(phi)*sech(phi)*(d2phidt2/2.0-dphidt*dphidt*tanh(phi))*lilf+&
          (lils*LOG(1.0+R(2)/lils)*(1.0-F*G)+R(2)*F*G+ &
          R(2)*(1+tanh(phi))/2.0)*d2lilfdt2+R(2)*d2lilgdt2

d2X0dxdt = 2*delta*dY0dt*R(1)*R(3)*oneuponden*oneuponden-delta*d2Y0dXdt*R(3)*oneuponden
d2X0dtdx = d2X0dxdt
d2X0dydt = (-delta)*d2Y0dydt*R(3)*oneuponden+2*delta*dY0dt*ay*ay*(R(2)-yh)*R(3)*oneuponden*oneuponden
d2X0dtdy = d2X0dydt
d2X0dzdt = (-delta)*dY0dt*(den-2*R(3)*R(3))*oneuponden*oneuponden-delta*d2Y0dzdt*R(3)*oneuponden
d2X0dtdz = d2X0dzdt
d2X0dt2 = (-delta)*d2y0dt2*R(3)*oneuponden

d2Z0dxdt = delta*dY0dt*(den-2*R(1)*R(1))*oneuponden*oneuponden+delta*d2Y0dXdt*R(1)*oneuponden
d2Z0dtdx = d2Z0dxdt
d2Z0dydt = delta*d2Y0dydt*R(1)*oneuponden-2*delta*dY0dt*R(1)*ay*ay*(R(2)-yh)*oneuponden*oneuponden
d2Z0dtdy = d2Z0dydt
d2Z0dzdt = (-delta)*dY0dt*2*R(1)*R(3)*oneuponden*oneuponden+delta*d2Y0dZdt*R(1)*oneuponden
d2Z0dtdz = d2Z0dzdt
d2Z0dt2 = delta*d2Y0dt2*R(1)*oneuponden

END SUBROUTINE ABSETUP

END MODULE CMT_fields

