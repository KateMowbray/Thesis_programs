MODULE CMT_fields2half

!km: Adjusted to cover 2.5D or 2D cases from 3D case so a lot commented out
!km: Mostly the same as CMT_fields but uses magnetic potential instead of explicitly defined field
!    to describe the final state of the field

  USE global
  USE M_products

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: CMTFIELDS2half
  PRIVATE :: SUB1_X0Y0Z0, dA0, SUB2_X0Y0Z0,ddA0, SUB3_X0Y0Z0
  
  CONTAINS 
!-----------------------------------------------------------------------------!   
SUBROUTINE CMTFIELDS2half(R,T,E,B,DBDX,DBDY,DBDZ,DBDT,DEDX,DEDY,DEDZ,DEDT,Vf)

 REAL(num)	:: X0,Y0,Z0,dX0dX,dX0dY,dX0dZ,dY0dX,dY0dY,dY0dZ,dZ0dX,dZ0dY,dZ0dZ,dX0dt,dY0dt,dZ0dt
 REAL(num)	:: d2X0dXdY,d2X0dYdX,d2X0dXdZ,d2X0dZdX,d2X0dYdZ,d2X0dZdY,d2X0dX2,d2X0dY2,d2X0dZ2
 REAL(num)	:: d2Y0dXdY,d2Y0dYdX,d2Y0dXdZ,d2Y0dZdX,d2Y0dYdZ,d2Y0dZdY,d2Y0dX2,d2Y0dY2,d2Y0dZ2
 REAL(num)      :: d2Z0dXdY,d2Z0dYdX,d2Z0dXdZ,d2Z0dZdX,d2Z0dYdZ,d2Z0dZdY,d2Z0dX2,d2Z0dY2,d2Z0dZ2
 REAL(num)	:: d2X0dtdX, d2X0dXdt, d2X0dtdY, d2X0dYdt, d2X0dtdZ, d2X0dZdt
 REAL(num)	:: d2Y0dtdX, d2Y0dXdt, d2Y0dtdY, d2Y0dYdt, d2Y0dtdZ, d2Y0dZdt
 REAL(num)      :: d2Z0dtdX, d2Z0dXdt, d2Z0dtdY, d2Z0dYdt, d2Z0dtdZ, d2Z0dZdt
 REAL(num)	:: d2X0dt2, d2Y0dt2, d2Z0dt2
 REAL(num)      :: a1,a2
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
! Find B0 from eq (36) in Grady paper
! a1 = SQRT((X0+0.5)*(X0+0.5)+(Y0+1)*(Y0+1)+Z0*Z0)
! a2 = SQRT((X0-0.5)*(X0-0.5)+(Y0+1)*(Y0+1)+Z0*Z0)

 !B0(1) = (X0+0.5)/a1/a1/a1 - (X0-0.5)/a2/a2/a2
 !B0(2) = (Y0+1)/a1/a1/a1 - (Y0+1)/a2/a2/a2
 !B0(3) = Z0/a1/a1/a1 - Z0/a2/a2/a2

 !B0 = B0 * c1  ! km: Multiply by factor c1 in (36) of Grady (2009)
! km: Set up the vector X
 !X(1)=X0
 !X(2)=Y0
 !X(3)=Z0

 !dXdx(1)=dX0dx
 !dXdx(2)=dY0dx
 !dXdx(3)=dZ0dx

 !dXdy(1)=dX0dy
 !dXdy(2)=dY0dy
 !dXdy(3)=dZ0dy

 !dXdz(1)=dX0dz
 !dXdz(2)=dY0dz
 !dXdz(3)=dZ0dz

 CALL dA0(X0,Y0,dA0dX0, dA0dY0) 
 
! JT: B(1) and B(2) are equations 28 and 29 of Guiliani et al. (2005)
 B(1) = oneuponL *  ( dA0dX0 * dX0dY + dA0dY0 * dY0dY )
 B(2) = oneuponL * ( -(dA0dX0 * dX0dX + dA0dY0 * dY0dX) )
 B(3) = oneuponL*((dY0dx*dZ0dy-dZ0dX*dY0dY)*(dA0dY0) -&
                  (dZ0dx*dX0dy-dX0dx*dZ0dy)*(dA0dX0))+&
                  (dX0dx*dY0dY-dY0dx*dX0dy)*Bzfinal

 !B(1)= DOT(CROSS(dXdy,dXdz),B0)*oneuponL
 !B(2)= DOT(CROSS(dXdz,dXdx),B0)*oneuponL
 !B(3)= DOT(CROSS(dXdx,dXdy),B0)*oneuponL
 
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

 CALL SUB2_X0Y0Z0(R,T,Y0,dY0dY,d2X0dXdY,d2X0dYdX,d2X0dXdZ,d2X0dZdX,d2X0dYdZ,&
                  d2X0dZdY,d2X0dX2,d2X0dY2,d2X0dZ2,d2Y0dXdY,d2Y0dYdX,d2Y0dXdZ,&
                  d2Y0dZdX,d2Y0dYdZ,d2Y0dZdY,d2Y0dX2,d2Y0dY2,d2Y0dZ2,d2Z0dXdY,&
                  d2Z0dYdX,d2Z0dXdZ,d2Z0dZdX,d2Z0dYdZ,d2Z0dZdY,d2Z0dX2,&
                  d2Z0dY2,d2Z0dZ2)
 CALL SUB3_X0Y0Z0(R,T,dY0dt,d2X0dtdX,d2X0dXdt,d2X0dtdY,d2X0dYdt,d2X0dtdZ,&
                  d2X0dZdt,d2Y0dtdX,d2Y0dXdt,d2Y0dtdY,d2Y0dYdt,d2Y0dtdZ,&
                  d2Y0dZdt,d2Z0dtdX,d2Z0dXdt,d2Z0dtdY,d2Z0dYdt,d2Z0dtdZ,&
                  d2Z0dZdt,d2X0dt2,d2Y0dt2,d2Z0dt2 )
 CALL ddA0(X0,Y0,d2A0dx02,d2A0dy02,d2A0dx0dy0,d2A0dy0dx0)

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

!!!! Check: LHS and RHS should be equal !!!!!!
 !print*,'1-LHS=',(der_det+Vf(1)*xder_det+Vf(2)*yder_det+Vf(3)*zder_det)/DETERMINANT
 !print*,'2-RHS=',-(dVfdx(1)+dVfdy(2)+dVfdz(3))
 !print*,'****'

!!!! Magnetic field derivatives in dimensional form !!!
! Have to check which things are dimensional here, B0 is not
 !dBxdx
 dBdx(1)=oneuponL*oneuponL*((d2A0dX02*dX0dx+d2A0dY0dX0*dY0dx)*dX0dy+dA0dx0*d2X0dxdy+&
            (d2A0dX0dY0*dX0dx+d2A0dY02*dY0dx)*dY0dy + dA0dY0*d2Y0dxdy)
 !dBxdy
 dBdy(1)=oneuponL*oneuponL*((d2A0dX02*dX0dy+d2A0dY0dX0*dY0dy)*dX0dy+dA0dX0*d2X0dy2+&
            (d2A0dX0dY0*dX0dy + d2A0dY02*dY0dy)*dY0dy + dA0dY0*d2Y0dy2)
 !dBxdz
 dBdz(1)=0.
 !dBydx
 dBdx(2)=-oneuponL*oneuponL*((d2A0dX02*dX0dx+d2A0dY0dX0*dY0dx)*dX0dx+dA0dX0*d2X0dx2+&
            (d2A0dX0dY0*dX0dx+d2A0dY02*dY0dx)*dY0dx + dA0dY0*d2Y0dx2)
 !dBydy
 dBdy(2)=-oneuponL*oneuponL*((d2A0dX02*dX0dy+d2A0dY0dX0*dY0dy)*dX0dx+dA0dX0*d2X0dydx +&
            (d2A0dX0dY0*dX0dy+d2A0dY02*dY0dy)*dY0dx + dA0dY0*d2Y0dydx)
 !dBydz
 dBdz(2)=0.
 !dBzdx
 dBdx(3)=oneuponL*oneuponL*((d2Y0dx2*dZ0dy+dY0dx*d2Z0dydx-d2Z0dx2*dY0dy-dZ0dx*d2Y0dydx)*&
                            (dA0dY0) +&
         (dY0dx*dZ0dy-dZ0dx*dY0dy)*(d2A0dY02*dY0dx+d2A0dY0dX0*dX0dx) -&
         (d2Z0dx2*dX0dy+dZ0dx*d2X0dydx-d2X0dx2*dZ0dy-dX0dx*d2Z0dydx)*(dA0dX0) -&
         (dZ0dx*dX0dy-dX0dx*dZ0dy)*(d2A0dX02*dX0dx+d2A0dX0dY0*dY0dx))+&
         oneuponL*(d2X0dx2*dY0dy+dX0dx*d2Y0dydx-d2Y0dx2*dX0dy-dY0dx*d2X0dydx)*Bzfinal

 !dBzdy
 dBdy(3)=oneuponL*oneuponL*&
         ((d2Y0dxdy*dZ0dy+dY0dx*d2Z0dy2-d2Z0dxdy*dY0dy-dZ0dx*d2Y0dy2)*(dA0dY0)+&
         (dY0dx*dZ0dy-dZ0dx*dY0dy)*(d2A0dY0dX0*dX0dY+d2A0dY02*dY0dY)-&
         (d2Z0dxdy*dX0dy+dZ0dx*d2X0dy2-d2X0dxdy*dZ0dy-dX0dx*d2Z0dy2)*(dA0dX0)-&
         (dZ0dx*dX0dy-dX0dx*dZ0dy)*(d2A0dX02*dX0dY+d2A0dX0dY0*dY0dY))+&
         oneuponL*(d2X0dxdy*dY0dy+dX0dx*d2Y0dy2-d2Y0dxdy*dX0dy-dY0dx*d2X0dy2)*Bzfinal
 !dBzdz
 dBdz(3)=0.0
 !dBxdt
 dBdt(1)=oneuponTL*(d2A0dX02*dX0dt*dX0dy + d2A0dY0dX0*dY0dt*dX0dy+&
                    d2A0dX0dY0*dX0dt*dY0dy + d2A0dY02*dY0dt*dY0dy+&
                    dA0dX0*d2X0dtdy + dA0dY0*d2Y0dtdy)
 !dBydt
 dBdt(2)=oneuponTL*(d2A0dX02*dX0dt*dX0dx + d2A0dY0dX0*dY0dt*dX0dx+&
                    d2A0dX0dY0*dX0dt*dY0dx + d2A0dY02*dY0dt*dY0dx+&
                    dA0dX0*d2X0dtdx + dA0dY0*d2Y0dtdx)
 !dBzdt
 dBdt(3)=oneuponTL*&
         ((d2Y0dxdt*dZ0dy+dY0dx*d2Z0dydt-d2Z0dxdt*dY0dy-dZ0dx*d2Y0dydt)*(dA0dY0*dY0dy+dA0dX0*dX0dy)+&
         (dY0dx*dZ0dy-dZ0dx*dY0dy)*(dA0dY0*d2Y0dydt+dY0dy*(d2A0dX0dY0*dX0dt+d2A0dY02*dY0dt)+&
         dA0dX0*d2X0dydt+dX0dy*(d2A0dX02*dX0dt+d2A0dX0dY0*dY0dt))-&
         (d2Z0dxdt*dX0dy+dZ0dx*d2X0dydt-d2X0dxdt*dZ0dy-dX0dx*d2Z0dydt)*(dA0dX0*dX0dx+dA0dY0*dY0dx)-&
         (dZ0dx*dX0dy-dX0dx*dZ0dy)*(dA0dX0*d2X0dxdt+dX0dx*(d2A0dX02*dX0dt+d2A0dX0dY0*dY0dt)+&
         dA0dY0*d2Y0dxdt+dY0dx*(d2A0dX0dY0*dX0dt+d2A0dY02*dY0dt)))+&
         oneuponT*(d2X0dxdt*dY0dy+dX0dx*d2Y0dydt-d2Y0dxdt*dX0dy-dY0dx*d2X0dydt)*Bzfinal

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

END SUBROUTINE CMTFIELDS2half
!-------------------------------------------------------------------------------------------!
SUBROUTINE SUB1_X0Y0Z0(R,T,X0,Y0,Z0,dX0dX,dX0dY,dX0dZ,dY0dX,dY0dY,dY0dZ,dZ0dX,dZ0dY,dZ0dZ,dX0dt,dY0dt,dZ0dt)
!JT+ reads in a given position vector at time t, 
! and calculates Lagrangian transformation, x0(=x_inf),y0(=y_inf),z0(=z_inf)
! based on eq. (36) of Guiliani et al (2005)

 REAL(num), DIMENSION(3), INTENT(IN)	:: R
 REAL(num), INTENT(IN)  		:: T
 REAL(num), INTENT(OUT)			:: X0, Y0,Z0, dX0dX, dX0dY,dX0dZ, dY0dX, dY0dY,dY0dZ,dZ0dX,dZ0dY,dZ0dZ, dX0dt, dY0dt,dZ0dt
 REAL(num)				:: a=0.9,b=0.9, delta=0.0,a3d=1.0

 !!!!!!  please notice Lv/L in the expressions  !!!!!

 IF ( (1.+(R(2)/((cc*T)**esp))) .le.0.) THEN
  PRINT *, '(1.+(R(2)/((cc*T)**esp)))=',(1.+(R(2)/((cc*T)**esp)))
  PRINT *, 'R(2)=',R(2)
  PRINT *,'T=',T
 END IF

 Y0=(cc*T)**esp*log (1.+(R(2)/((cc*T)**esp)))*  &	
      (1.+tanh((R(2)-Lv/Lscl)*a))*0.5 +           &
       (1.-tanh((R(2)-Lv/Lscl)*b))*0.5*R(2)

 X0=R(1)

 Z0=R(3)+delta*(Y0-R(2))*R(1)/(a3d*a3d+R(1)*R(1))

!print*,X0,Y0

!JT: are these boundary conditions?

 dY0dX=0.
 dY0dY=0.5*(1.+tanh((R(2)-Lv/Lscl)*a))/(1.+R(2)/(cc*T)**esp)+        &
      0.5*(cc*T)**esp*log(1.+R(2)/(cc*T)**esp)*               &
      (1.-tanh((R(2)-Lv/Lscl)*a)**2)*a-0.5*                         &
      (1.-tanh((R(2)-Lv/Lscl)*b)**2)*b*R(2)+0.5-0.5*tanh((R(2)-Lv/Lscl)*b) 
 dY0dZ=0.

 dX0dX=1.0
 dX0dY=0.0
 dX0dZ=0.0

 dZ0dX=delta*(Y0-R(2))*(a3d*a3d+R(1)*R(1)+R(3)*R(3)-2*R(1)*R(1))/   &
       (a3d*a3d+R(1)*R(1)+R(3)*R(3))**2
 dZ0dY=delta*(dY0dY-1.0)*R(1)/(a3d*a3d+R(1)*R(1))
 dZ0dZ=1.0

 dY0dt=0.5*(cc*T)**esp*esp*log(1.+R(2)/(cc*T)**esp)*          &
      (1.+tanh((R(2)-Lv/Lscl)*a))/T-0.5*R(2)*esp*                  &
      (1.+tanh((R(2)-Lv/Lscl)*a))/(T*(1.+R(2)/(cc*T)**esp)) 
 dX0dt=0.0
 dZ0dt=delta*dY0dt*R(1)/(a3d*a3d+R(1)*R(1))         
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
SUBROUTINE SUB2_X0Y0Z0(R,T,Y0,dY0dY,d2X0dXdY,d2X0dYdX,d2X0dXdZ,d2X0dZdX,d2X0dYdZ,&
                       d2X0dZdY,d2X0dX2,d2X0dY2,d2X0dZ2,d2Y0dXdY,d2Y0dYdX, &
                       d2Y0dXdZ,d2Y0dZdX,d2Y0dYdZ,d2Y0dZdY,d2Y0dX2,d2Y0dY2,&
                       d2Y0dZ2,d2Z0dXdY,d2Z0dYdX,d2Z0dXdZ,d2Z0dZdX,d2Z0dYdZ,&
                       d2Z0dZdY,d2Z0dX2,d2Z0dY2,d2Z0dZ2)
 
 REAL(num), DIMENSION(3), INTENT(IN) :: R
 REAL(num), INTENT(IN)  	:: T, DY0dY !included to save copying long statement from SUB1_X0Y0Z0
 REAL(num), INTENT(OUT) 	:: d2X0dXdY,d2X0dYdX,d2X0dXdZ,d2X0dZdX,d2X0dYdZ,d2X0dZdY,d2X0dX2,d2X0dY2,d2X0dZ2
 REAL(num), INTENT(OUT) 	:: d2Y0dXdY,d2Y0dYdX,d2Y0dXdZ,d2Y0dZdX,d2Y0dYdZ,d2Y0dZdY,d2Y0dX2,d2Y0dY2,d2Y0dZ2
 REAL(num), INTENT(OUT)         :: d2Z0dXdY,d2Z0dYdX,d2Z0dXdZ,d2Z0dZdX,d2Z0dYdZ,d2Z0dZdY,d2Z0dX2,d2Z0dY2,D2Z0dZ2
 REAL(num)              	:: y, den, Y0
 REAL(num)              	:: a=0.9,b=0.9,delta=0.0,a3d=1.0

 y=R(2)
 den = a3d*a3d+R(1)*R(1) !km: denominator from (37) and (38) of Grady(2009)
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

 d2X0dXdY = 0.
 d2X0dYdX = 0.
 d2X0dXdZ = 0.
 d2X0dZdX = 0.
 d2X0dYdZ = 0.
 d2X0dZdY = 0.
 d2X0dX2 = 0.
 d2X0dY2 = 0.
 d2X0dZ2 = 0.

 d2Z0dXdY = delta*(dY0dY-1.0)*(den-2*R(1)*R(1))/den/den
 d2Z0dYdX = d2Z0dXdY
 d2Z0dXdZ = 0.
 d2Z0dZdX = d2Z0dXdZ
 d2Z0dYdZ = 0.
 d2Z0dZdY = d2Z0dYdZ
 d2Z0dX2 = delta*(Y0-R(2))*(8*R(1)*R(1)*R(1)-6*R(1)*den)/den/den/den
 d2Z0dY2 = delta*d2Y0dY2*R(1)/den
 d2Z0dZ2 = 0.

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
 REAL(num)				:: a=0.9, delta=0.0, a3d=1.0
 
y=R(2)
den=a3d*a3d+R(1)*R(1)

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

d2X0dtdx = 0.
d2X0dxdt = 0.
d2X0dydt = 0.
d2X0dtdy = 0.
d2X0dtdz = 0.
d2X0dzdt = 0.
d2X0dt2 = 0.

d2Z0dtdx = delta*dY0dt*(den-2*R(1)*R(1))/den/den
d2Z0dxdt = d2Z0dtdx
d2Z0dtdy = delta*d2Y0dtdy*R(1)/den
d2Z0dydt = d2Z0dtdy
d2Z0dtdz = 0.
d2Z0dzdt = 0.
d2Z0dt2 = delta*d2Y0dt2*R(1)/den


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

END MODULE CMT_fields2half

